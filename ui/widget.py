"""Standalone Zernike coefficient editor widget."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .coefficient_table import CoefficientTablesPanel

_DEFAULT_WIDTH = 960
_DEFAULT_HEIGHT = 820
_MIN_WIDTH = 720
_MIN_HEIGHT = 600
_SYNC_STYLE = {
    True: ("up to date", "#1b7f3a"),
    False: ("not up to date", "#c0392b"),
    None: ("unresolved", "#b8860b"),
}


def _sync_line(prefix, entry):
    """Status text and color for one JSON or backup entry."""
    word, color = _SYNC_STYLE[entry["sync"]]
    when = entry["timestamp"]
    if when is None:
        return f"{prefix} — {word}", color
    return f"{prefix} {when.strftime('%d.%m.%Y %H:%M:%S')} — {word}", color


class ZernikesWidget(QWidget):
    """Edit Zernikes coefficients, load/save JSON, and regenerate phase masks."""

    phaseMaskChanged = Signal()

    def __init__(self, zernikes, parent=None):
        super().__init__(parent)
        self._z = zernikes
        self.setWindowTitle("Zernikes")
        self.resize(_DEFAULT_WIDTH, _DEFAULT_HEIGHT)
        self.setMinimumSize(_MIN_WIDTH, _MIN_HEIGHT)

        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse_load)
        reload_btn = QPushButton("Reload from file")
        reload_btn.clicked.connect(self._reload_from_file)
        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self._save)
        save_as_btn = QPushButton("Save As…")
        save_as_btn.clicked.connect(self._save_as)

        self._reload_btn = reload_btn

        buttons_row = QHBoxLayout()
        buttons_row.addWidget(browse_btn)
        buttons_row.addWidget(reload_btn)
        buttons_row.addWidget(save_btn)
        buttons_row.addWidget(save_as_btn)
        self._json_status = QLabel()
        buttons_row.addWidget(self._json_status)
        buttons_row.addStretch(1)

        self._path_label = QLabel()
        self._path_label.setWordWrap(True)

        params_box = QGroupBox("Parameters")
        params_form = QFormLayout(params_box)
        z = self._z
        params_form.addRow("max radial order", QLabel(str(z.max_zernike_radial_number)))
        params_form.addRow(
            "mask size (px)",
            QLabel(f"{z.mask_size_x_in_pixels} × {z.mask_size_y_in_pixels}"),
        )
        params_form.addRow(
            "aperture diameter (m)",
            QLabel(f"{z.aperture_diameter_in_m:g}"),
        )
        params_form.addRow(
            "pixel size (m)",
            QLabel(f"{z.SLM_pixel_size_in_um:g}"),
        )

        self._coef_tables = CoefficientTablesPanel(z)
        self._coef_tables.coefficientsChanged.connect(self._regenerate_phase)

        layout = QVBoxLayout(self)
        layout.addLayout(buttons_row)
        layout.addWidget(self._path_label)
        layout.addWidget(params_box)
        layout.addWidget(self._coef_tables, stretch=1)

        self.refresh_from_zernikes()

    
    def refresh_from_zernikes(self):
        """Sync all controls from the Zernikes instance. Call this when Zenike Instance was modified outside of the widget."""
        #!This is necessary to ensure that the widget is in sync with the Zenike Instance.
        #!For example, if the Zenike Instance was modified outside of the widget, the widget will not be aware of the change.
        #!This function is called when the Zenike Instance is loaded from a file, or when the Zenike Instance is modified.

        z = self._z

        if z.coefs_path:
            self._path_label.setText(f"Loaded: {z.coefs_path}")
        else:
            self._path_label.setText("No coefficients file loaded")
        self._reload_btn.setEnabled(z.coefs_path is not None)
        self._coef_tables.refresh_from_zernikes()
        self._apply_sync_status()

    def _apply_sync_status(self):
        """Paint JSON and backup status from the engine check."""
        res = self._z.coeffs_status_sync()
        text, color = _sync_line("last saved", res["json"])
        self._json_status.setText(text)
        self._json_status.setStyleSheet(f"color: {color};")
        text, color = _sync_line("last backup", res["backup"])
        self._coef_tables.apply_backup_status(text, color)

    def _regenerate_phase(self):
        """Regenerate the phase mask from the current coefficients and update the display."""

        self._z.calculate_znk_phase_mask(restore_coefs=False)
        self._coef_tables.refresh_phase_viewport()
        self._apply_sync_status()
        self.phaseMaskChanged.emit()

    def _after_load(self):
        self._z.write_backup_of_zernike_coefs()
        self._z.calculate_znk_phase_mask(restore_coefs=False)
        self.refresh_from_zernikes()
        self._coef_tables.refresh_phase_viewport()
        self.phaseMaskChanged.emit()

    def _browse_load(self):
        start = self._z.coefs_path or ""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Zernike coefficients",
            start,
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        self._z.load_coefs(path)
        self._after_load()

    def _reload_from_file(self):
        if not self._z.coefs_path:
            return
        self._z.load_coefs(self._z.coefs_path)
        self._after_load()

    def _save(self):
        if self._z.coefs_path:
            self._z.save_coefs(self._z.coefs_path)
            self.refresh_from_zernikes()
            return
        self._save_as()

    def _save_as(self):
        start = self._z.coefs_path or ""
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Zernike coefficients",
            start,
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        if not path.lower().endswith(".json"):
            path = f"{path}.json"
        self._z.save_coefs(path)
        self.refresh_from_zernikes()
