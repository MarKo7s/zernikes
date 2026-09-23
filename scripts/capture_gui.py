#!/usr/bin/env python3
"""Grab a live Zernikes GUI screenshot for README (docs/images/gui.png).

Uses ZernikesWidget on an in-memory Zernikes instance (no SLM, no saved file).
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# This script lives in `scripts/`; repo root is its parent directory.
REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DEFAULT = REPO_ROOT / "docs" / "images" / "gui.png"


def build_window():
    """Return (window, teardown) for a clean coefficient editor."""
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))

    from zernikes import Zernikes
    from zernikes.ui import ZernikesWidget

    z = Zernikes(
        max_zernike_radial_number=4,
        mask_size_x_in_pixels=64,
        mask_size_y_in_pixels=64,
        coefs_path=None,
    )
    gui = ZernikesWidget(z)
    gui.show()

    def teardown() -> None:
        gui.close()

    return gui, teardown


def main() -> None:
    parser = argparse.ArgumentParser(description="Capture GUI screenshot for README.")
    parser.add_argument("-o", "--output", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--wait-ms", type=int, default=4000)
    parser.add_argument(
        "--offscreen",
        action="store_true",
        help="Use QT_QPA_PLATFORM=offscreen (may blank pyqtgraph on Windows).",
    )
    args = parser.parse_args()

    if args.offscreen:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    import pyqtgraph as pg
    from PySide6.QtWidgets import QApplication

    pg.setConfigOptions(useOpenGL=False, antialias=True)

    app = QApplication.instance() or QApplication(sys.argv)
    window, teardown = build_window()
    window.resize(960, 820)
    window.show()
    window.raise_()
    window.activateWindow()

    deadline = time.perf_counter() + max(args.wait_ms, 200) / 1000.0
    while time.perf_counter() < deadline:
        app.processEvents()
        time.sleep(0.02)
    app.processEvents()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pix = window.grab()
    if pix.isNull() or pix.width() < 10:
        teardown()
        raise SystemExit("grab() returned an empty pixmap")
    if not pix.save(str(args.output), "PNG"):
        teardown()
        raise SystemExit(f"Failed to write {args.output}")

    teardown()
    print(f"Wrote {args.output}")
    app.quit()


if __name__ == "__main__":
    main()
