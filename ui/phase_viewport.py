"""Optional pyqtgraph viewport for Zernike H/V phase masks."""

import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QVBoxLayout, QWidget

PI = math.pi


def _phase_colormap() -> pg.ColorMap:
    """Grayscale: -π black → +π white."""
    pos = np.array([0.0, 1.0])
    color = np.array([[0, 0, 0, 255], [255, 255, 255, 255]], dtype=np.ubyte)
    return pg.ColorMap(pos, color)


class PhaseViewportWindow(QWidget):
    """Side-by-side H/V phase maps with a shared phase color bar."""

    def __init__(self, zernikes, parent=None):
        super().__init__(parent)
        self._z = zernikes
        self._did_initial_fit = False

        self.setWindowTitle("Zernike phase masks")
        self.resize(920, 480)

        pg.setConfigOptions(useOpenGL=False, antialias=False)

        self._graphics = pg.GraphicsLayoutWidget()
        self._plot_h = self._graphics.addPlot(row=0, col=0, title="H phase")
        self._plot_v = self._graphics.addPlot(row=0, col=1, title="V phase")

        self._img_h = pg.ImageItem(axisOrder="row-major")
        self._img_v = pg.ImageItem(axisOrder="row-major")
        self._plot_h.addItem(self._img_h)
        self._plot_v.addItem(self._img_v)

        cmap = _phase_colormap()
        for image in (self._img_h, self._img_v):
            image.setColorMap(cmap)
            image.setLevels([-PI, PI])

        for plot in (self._plot_h, self._plot_v):
            plot.setAspectLocked(True)
            plot.invertY(True)
            plot.showGrid(x=False, y=False)

        self._colorbar = pg.ColorBarItem(
            values=(-PI, PI),
            colorMap=cmap,
            label="phase (rad)",
            width=24,
        )
        self._graphics.addItem(self._colorbar, row=0, col=2)
        self._colorbar.setImageItem(self._img_h)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(self._graphics)

        self.refresh()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        QTimer.singleShot(0, self._fit_views)

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        QTimer.singleShot(0, self._fit_views)

    def _fit_views(self) -> None:
        n = self._z.N
        for plot in (self._plot_h, self._plot_v):
            plot.setLimits(
                xMin=0,
                xMax=n,
                yMin=0,
                yMax=n,
                minXRange=4,
                minYRange=4,
                maxXRange=n,
                maxYRange=n,
            )
            plot.autoRange(padding=0)

    def refresh(self) -> None:
        """Update images from the current Zernikes phase buffers."""
        phase_h = np.asarray(self._z.znk_phase_h, dtype=np.float32)
        phase_v = np.asarray(self._z.znk_phase_v, dtype=np.float32)
        self._img_h.setImage(phase_h, autoLevels=False)
        self._img_v.setImage(phase_v, autoLevels=False)
        if not self._did_initial_fit:
            self._did_initial_fit = True
            QTimer.singleShot(0, self._fit_views)
