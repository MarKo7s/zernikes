"""Standalone Zernike coefficient editor widget."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from .coefficient_table import CoefficientTablesPanel

_DEFAULT_WIDTH = 960
_DEFAULT_HEIGHT = 820
_MIN_WIDTH = 720
_MIN_HEIGHT = 600


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

        self._h_center_x = QSpinBox()
        self._h_center_y = QSpinBox()
        self._v_center_x = QSpinBox()
        self._v_center_y = QSpinBox()
        for spin in (self._h_center_x, self._h_center_y, self._v_center_x, self._v_center_y):
            spin.setRange(0, z.N - 1)
            spin.setSuffix(" px")
            spin.valueChanged.connect(self._on_center_changed)

        h_center_box = QGroupBox("Mask center H")
        h_center_form = QFormLayout(h_center_box)
        h_center_form.addRow("x", self._h_center_x)
        h_center_form.addRow("y", self._h_center_y)

        v_center_box = QGroupBox("Mask center V")
        v_center_form = QFormLayout(v_center_box)
        v_center_form.addRow("x", self._v_center_x)
        v_center_form.addRow("y", self._v_center_y)

        centers_row = QHBoxLayout()
        centers_row.addWidget(h_center_box, stretch=1)
        centers_row.addWidget(v_center_box, stretch=1)

        self._coef_tables = CoefficientTablesPanel(z)
        self._coef_tables.coefficientsChanged.connect(self._regenerate_phase)

        layout = QVBoxLayout(self)
        layout.addLayout(buttons_row)
        layout.addWidget(self._path_label)
        layout.addWidget(params_box)
        layout.addLayout(centers_row)
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

        for spin, attr in (
            (self._h_center_x, "h_mask_center_x"),
            (self._h_center_y, "h_mask_center_y"),
            (self._v_center_x, "v_mask_center_x"),
            (self._v_center_y, "v_mask_center_y"),
        ):
            spin.blockSignals(True)
            spin.setValue(int(getattr(z, attr, 0)))
            spin.blockSignals(False)

        self._coef_tables.refresh_from_zernikes()

    def _regenerate_phase(self):
        """Regenerate the phase mask from the current coefficients and update the display."""

        self._z.calculate_znk_phase_mask(restore_coefs=False)
        self._coef_tables.refresh_phase_viewport()
        self.phaseMaskChanged.emit()

    def _on_center_changed(self):
        self._z.h_mask_center_x = self._h_center_x.value()
        self._z.h_mask_center_y = self._h_center_y.value()
        self._z.v_mask_center_x = self._v_center_x.value()
        self._z.v_mask_center_y = self._v_center_y.value()

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
