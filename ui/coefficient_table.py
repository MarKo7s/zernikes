"""Grouped editable Zernike coefficient tables by aberration family."""

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

try:
    from zernike_catalog import ABERRATION_GROUPS, get_mode_name, iter_modes_by_group
except ImportError:
    from ..zernike_catalog import ABERRATION_GROUPS, get_mode_name, iter_modes_by_group

_COL_INDEX = 0
_COL_N = 1
_COL_M = 2
_COL_NAME = 3
_COL_H = 4
_COL_V = 5
_HEADERS = ("index", "n", "m", "name", "H", "V")
_ROW_HEIGHT = 22


def _read_only_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(Qt.ItemFlag.ItemIsEnabled)
    return item


def _editable_item(value: float) -> QTableWidgetItem:
    item = QTableWidgetItem(f"{value:.6g}")
    item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsEditable)
    return item


def _configure_table_columns(table: QTableWidget) -> None:
    """Fix coefficient column widths; name column absorbs extra horizontal space."""
    table.resizeColumnsToContents()
    header = table.horizontalHeader()
    header.setStretchLastSection(False)
    coef_width = max(table.columnWidth(_COL_H), table.columnWidth(_COL_V), 72)
    for col in range(table.columnCount()):
        if col == _COL_NAME:
            header.setSectionResizeMode(col, QHeaderView.ResizeMode.Stretch)
        else:
            header.setSectionResizeMode(col, QHeaderView.ResizeMode.Fixed)
            if col in (_COL_H, _COL_V):
                table.setColumnWidth(col, coef_width)
    table.verticalHeader().setDefaultSectionSize(_ROW_HEIGHT)
    table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)


class _ModeGroupTable(QGroupBox):
    """One aberration family table inside a group box."""

    coefficientsChanged = Signal()

    def __init__(self, title: str, modes: list[tuple[int, int, int]], zernikes, parent=None):
        super().__init__(title, parent)
        self._z = zernikes
        self._row_indices = [idx for idx, _, _ in modes]

        self._table = QTableWidget()
        self._table.setColumnCount(len(_HEADERS))
        self._table.setHorizontalHeaderLabels(list(_HEADERS))
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        self._table.setFont(QFont("Consolas", 9))
        self._table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(self._table)

        self._build_rows(modes)
        _configure_table_columns(self._table)
        self._table.setFixedHeight(
            self._table.horizontalHeader().height() + _ROW_HEIGHT * len(modes) + 2
        )
        self._table.cellChanged.connect(self._on_cell_changed)

    def _build_rows(self, modes: list[tuple[int, int, int]]) -> None:
        self._table.setRowCount(len(modes))
        for row, (idx, n, m) in enumerate(modes):
            name = get_mode_name(n, m)
            self._table.setItem(row, _COL_INDEX, _read_only_item(str(idx)))
            self._table.setItem(row, _COL_N, _read_only_item(str(n)))
            self._table.setItem(row, _COL_M, _read_only_item(str(m)))
            self._table.setItem(row, _COL_NAME, _read_only_item(name))
            self._table.setItem(row, _COL_H, _editable_item(float(self._z.w_h[idx])))
            self._table.setItem(row, _COL_V, _editable_item(float(self._z.w_v[idx])))

    def refresh_from_zernikes(self) -> None:
        self._table.blockSignals(True)
        try:
            for row, idx in enumerate(self._row_indices):
                self._table.item(row, _COL_H).setText(f"{float(self._z.w_h[idx]):.6g}")
                self._table.item(row, _COL_V).setText(f"{float(self._z.w_v[idx]):.6g}")
        finally:
            self._table.blockSignals(False)

    def _on_cell_changed(self, row: int, col: int) -> None:
        if col not in (_COL_H, _COL_V):
            return

        item = self._table.item(row, col)
        if item is None:
            return

        try:
            value = float(item.text())
        except ValueError:
            self.refresh_from_zernikes()
            return

        idx = self._row_indices[row]
        value_f32 = np.float32(value)
        if col == _COL_H:
            self._z.w_h[idx] = value_f32
        else:
            self._z.w_v[idx] = value_f32

        item.setText(f"{value:.6g}")
        self.coefficientsChanged.emit()


class CoefficientTablesPanel(QWidget):
    """Scrollable stack of aberration-group coefficient tables."""

    coefficientsChanged = Signal()

    def __init__(self, zernikes, parent=None):
        super().__init__(parent)
        self._z = zernikes
        self._group_tables: list[_ModeGroupTable] = []
        self._viewport = None

        save_backup_btn = QPushButton("Save to backup")
        save_backup_btn.clicked.connect(self._save_to_backup)
        reload_backup_btn = QPushButton("Reload from backup")
        reload_backup_btn.clicked.connect(self._reload_from_backup)

        self._backup_status = QLabel()
        backup_box = QGroupBox("Backup")
        backup_row = QHBoxLayout(backup_box)
        backup_row.addWidget(save_backup_btn)
        backup_row.addWidget(reload_backup_btn)
        backup_row.addWidget(self._backup_status)

        viewport_btn = QPushButton("Launch phase mask viewport")
        viewport_btn.clicked.connect(self._launch_viewport)

        actions_row = QHBoxLayout()
        actions_row.addWidget(backup_box)
        actions_row.addStretch(1)
        actions_row.addWidget(viewport_btn)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self._scroll_content = QWidget()
        self._scroll_layout = QVBoxLayout(self._scroll_content)
        self._scroll_layout.setContentsMargins(0, 0, 0, 0)
        self._scroll.setWidget(self._scroll_content)

        self._build_group_tables()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(actions_row)
        layout.addWidget(self._scroll, stretch=1)

    def _build_group_tables(self) -> None:
        while self._scroll_layout.count():
            item = self._scroll_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._group_tables.clear()

        modes_by_group = iter_modes_by_group(self._z.max_zernike_radial_number)
        title_by_key = dict(ABERRATION_GROUPS)

        for key, _label in ABERRATION_GROUPS:
            modes = modes_by_group.get(key)
            if not modes:
                continue
            group_table = _ModeGroupTable(title_by_key[key], modes, self._z)
            group_table.coefficientsChanged.connect(self.coefficientsChanged.emit)
            self._group_tables.append(group_table)
            self._scroll_layout.addWidget(group_table)

        self._scroll_layout.addStretch(1)

    def refresh_from_zernikes(self) -> None:
        for group_table in self._group_tables:
            group_table.refresh_from_zernikes()

    def refresh_phase_viewport(self) -> None:
        if self._viewport is not None:
            self._viewport.refresh()

    def _launch_viewport(self) -> None:
        from .phase_viewport import PhaseViewportWindow

        if self._viewport is None:
            self._viewport = PhaseViewportWindow(self._z)
            self._viewport.destroyed.connect(self._on_viewport_destroyed)
        self._viewport.refresh()
        self._viewport.show()
        self._viewport.raise_()
        self._viewport.activateWindow()

    def _on_viewport_destroyed(self, *_args) -> None:
        self._viewport = None

    def apply_backup_status(self, text: str, color: str) -> None:
        """Show the backup timestamp line painted by the parent widget."""
        self._backup_status.setText(text)
        self._backup_status.setStyleSheet(f"color: {color};")

    def _save_to_backup(self) -> None:
        self._z.write_backup_of_zernike_coefs()
        self.coefficientsChanged.emit()

    def _reload_from_backup(self) -> None:
        self._z.read_backup_of_zernike_coefs()
        self.refresh_from_zernikes()
        self.coefficientsChanged.emit()


# Backward-compatible alias
CoefficientTable = CoefficientTablesPanel
