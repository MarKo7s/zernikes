# Changelog

All notable changes to zernikes are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.0.0] - 2026-09-16

### Added

- Pip-installable package layout (`pyproject.toml`, git-tag release via `scripts/release.py`).
- `Zernikes` basis and H/V phase-mask synthesis with JSON coefficient I/O, in-memory backup, and tilt/defocus scans.
- Optional PySide6 coefficient editor (`zernikes.ui`) with grouped aberration tables and pyqtgraph phase viewport.
