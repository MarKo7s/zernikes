# Zernikes

Modelab Zernike basis and phase-mask synthesis for SLM holography.

Add the repo folder to `sys.path`, then import the module file directly:

```python
from pathlib import Path
import sys

zernikes_dir = Path("/path/to/algorithms/zernikes")
sys.path.insert(0, str(zernikes_dir))

from zernikes import Zernikes
```

## Coefficient model

Weights live in unified vectors `w_h` and `w_v`, one entry per basis index:

| Index | (n, m) | Name (n ≤ 4) |
|------:|--------|--------------|
| 0 | (0, 0) | Piston |
| 1 | (1, −1) | Tilt X |
| 2 | (1, +1) | Tilt Y |
| 4 | (2, 0) | Defocus |
| … | … | … |

Named attributes are aliases into the array:

```python
z.piston_h = 0.1      # same as z.w_h[0] = 0.1
z.tilt_x_h = 0.05     # same as z.w_h[1] = 0.05
z.w_h[Zernikes.nm_to_idx(2, 0)] = 0.2  # defocus
```

## Quick start

```python
from zernikes import Zernikes

z = Zernikes(
    max_zernike_radial_number=6,
    mask_size_x_in_pixels=64,
    mask_size_y_in_pixels=64,
)

df = z.generate_table()       # pandas DataFrame if pandas is installed
z.calculate_znk_phase_mask()

phase_h = z.znk_phase_h       # real wrapped phase in [-π, π] radians
phase_v = z.znk_phase_v
```

## Aberration groups

Modes are classified in [`zernike_catalog.py`](zernike_catalog.py) for the UI and for group-wise optimization:

| Key | Rule |
|-----|------|
| `piston` | n = 0 |
| `tilt` | n = 1 |
| `defocus` | m = 0, n = 2 |
| `spherical` | m = 0, n ≥ 4 |
| `astigmatism` | \|m\| = 2 |
| `coma` | \|m\| = 1, n ≥ 3 |
| `foil` | \|m\| ≥ 3 |

```python
from zernike_catalog import iter_modes_by_group, ABERRATION_GROUPS

groups = iter_modes_by_group(z.max_zernike_radial_number)
# groups["defocus"] → [(4, 2, 0), ...]   each tuple is (index, n, m)

coefs, idxs = z.filter_coefficients_by_aberration_keys("H", "defocus", "foil")
# coefs = z.w_h[idxs]   idxs are sorted, unique basis indices
```

## In-memory backup

Snapshot and restore coefficient vectors without touching disk:

```python
z.write_backup_of_zernike_coefs()   # copy w_h/w_v → b_w_h/b_w_v
z.read_backup_of_zernike_coefs()  # restore from backup
```

File load in the widget also refreshes the backup after each load.

## Init vs scan workflow

**Normal init** (load or clean → backup → synthesize):

```python
z = Zernikes(..., coefs_path="coefs.json")  # sets z.coefs_path; or z.clean()
# __init__ already calls write_backup_of_zernike_coefs()
# and calculate_znk_phase_mask(restore_coefs=False)
```

**Tilt/defocus scan** (temporary offsets, then restore):

```python
z.modify_zernike_coefs(scan_tilt=True, tilt_x=0.01, tilt_y=0.0)
z.calculate_znk_phase_mask(restore_coefs=True)

z.modify_zernike_coefs(scan_defocus=True, defocus=0.05)
z.calculate_znk_phase_mask(restore_coefs=True)
```

## JSON coefficients

Save / load uses compact JSON (piston, tilts, and higher-order `w` for n ≥ 2).
If the file's `max_zernike_radial_number` differs from the instance, the instance wins;
extra file modes are ignored and missing modes default to zero (with warnings).

```python
z.save_coefs("coefs.json")    # updates z.coefs_path
z.load_coefs("coefs.json")
z.calculate_znk_phase_mask()
```

## SLM integration

Pass real phase masks to `pyLCOS` (not complex fields):

```python
slm.ModMuxMask.H.zernike = z.znk_phase_h
slm.ModMuxMask.V.zernike = z.znk_phase_v
slm.setmask(pol="HV")
```

Phase wrapping matches `pyLCOS._addMasks`: coefficients are summed in the basis,
then wrapped to `[-π, π]` via `((π·φ + π) mod 2π) − π`.

## UI widget

PySide6 editor for an existing `Zernikes` instance. Requires **PySide6**; the optional phase viewport also needs **pyqtgraph**.

Features:

- Load / save / reload JSON coefficients
- Coefficient tables grouped by aberration family (piston, tilt, defocus, …)
- In-memory **Save to backup** / **Reload from backup**
- Optional **Launch phase mask viewport** (H and V phase maps, grayscale −π…π)

```python
from zernikes import Zernikes
from ui.widget import ZernikesWidget

z = Zernikes(max_zernike_radial_number=6, mask_size_x_in_pixels=256, mask_size_y_in_pixels=256)
gui = ZernikesWidget(z)
gui.show()
```

When coefficients are changed **outside** the widget (e.g. from another script or optimizer), sync the GUI and phase:

```python
z.tilt_x_h = 0.05
gui.refresh_from_zernikes()   # updates tables and path label
gui._regenerate_phase()       # recalculates masks and open viewport
```

See [`tests/zernikes_widget.ipynb`](tests/zernikes_widget.ipynb) for a notebook launcher (`%gui qt`).

## Module layout

| File | Role |
|------|------|
| `zernikes.py` | `Zernikes` class: basis, phase synthesis, scans, `filter_coefficients_by_aberration_keys`, `generate_table` |
| `zernike_catalog.py` | Index ↔ (n, m) ↔ names; `aberration_group`, `iter_modes_by_group` |
| `utils.py` | JSON load/save and coefficient assembly |
| `ui/widget.py` | `ZernikesWidget` main editor |
| `ui/coefficient_table.py` | Grouped editable coefficient tables |
| `ui/phase_viewport.py` | Optional pyqtgraph H/V phase viewer |

## Static helpers

```python
Zernikes.nm_to_idx(2, 0)   # → 4 (defocus)
Zernikes.idx_to_nm(4)      # → (2, 0)
```
