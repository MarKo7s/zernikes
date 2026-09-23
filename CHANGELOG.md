# Changelog

All notable changes to zernikes are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.1.0] - 2026-09-23

### Changed

- Coefficient JSON now starts with a preamble (`schema`, `version`, `saved_at`). A bare `Zernikes()` loads `config/zernike_coeffs.json` when that file exists; `coefs_path=None` stays a clean start.
- Mask centers are no longer stored on `Zernikes` or in the coefficient file.
- `coeffs_status_sync()` reports whether current weights match the last JSON save and the RAM backup. The GUI shows that next to Save and in the backup box.

## [1.0.0] - 2026-09-16

### Added

- Pip-installable package layout (`pyproject.toml`, git-tag release via `scripts/release.py`).
- `Zernikes` basis and H/V phase-mask synthesis with JSON coefficient I/O, in-memory backup, and tilt/defocus scans.
- Optional PySide6 coefficient editor (`zernikes.ui`) with grouped aberration tables and pyqtgraph phase viewport.
