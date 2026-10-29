# Changelog

All notable changes to zernikes are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.2.0] - 2026-10-29

### Fixed

- Piston (basis row 0) is 1 inside the pupil, so its weight adds a uniform phase of `w·π` radians. The Modelab center subtraction had left that row at 0.

### Added

- In the coefficient editor, Up and Down step the selected H or V cell by `COEFFICIENT_KEY_STEP` (0.1).
- `calculate_znk_phase_mask` emits `onPhaseMaskChanged` after the masks are written.
- `ZernikesWidget.refresh_all_from_zernikes` (and its async form) rebuilds the phase mask from the current weights. `enable_user_interface` locks the editor during a remote run.
- `modify_zernike_coefs` takes `pol` (`"H"`, `"V"`, or `"HV"`) and applies tilt and defocus offsets to that polarization only.

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
