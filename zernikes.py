"""Zernike phase masks for SLM holography (Modelab basis convention).

Weights live in unified vectors ``w_h`` and ``w_v`` indexed by basis row.
``piston_h``, ``tilt_x_h``, and ``tilt_y_h`` are property aliases into ``w_h[0:3]``.
"""

from math import factorial
from pathlib import Path

import numexpr as ne
import numpy as np

#Double import so it can be used as module and in the notebook
try:
    from .utils import (
        default_coefficients_dict,
        dict_to_zernike_coefficients,
        load_coefficients_json,
        save_coefficients_json,
        zernike_coefficients_to_dict,
    )
    from .zernike_catalog import (
        IDX_DEFOCUS,
        expected_zernike_count,
        get_mode_name,
        idx_to_nm as catalog_idx_to_nm,
        iter_modes,
        nm_to_idx as catalog_nm_to_idx,
        iter_modes_by_group,
        ABERRATION_GROUPS
    )
except ImportError:
    from utils import (
        default_coefficients_dict,
        dict_to_zernike_coefficients,
        load_coefficients_json,
        save_coefficients_json,
        zernike_coefficients_to_dict,
    )
    from zernike_catalog import (
        IDX_DEFOCUS,
        expected_zernike_count,
        get_mode_name,
        idx_to_nm as catalog_idx_to_nm,
        iter_modes,
        nm_to_idx as catalog_nm_to_idx,
        iter_modes_by_group,
        ABERRATION_GROUPS
    )


def cart2pol(x, y):
    """Cartesian to polar (rho, phi)."""
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return rho, phi


class Zernikes:
    """Build Modelab Zernike basis on an SLM grid and synthesize H/V phase masks.

    Coefficients are stored in ``w_h`` and ``w_v`` (one entry per basis index).
    Named attributes ``piston_h``, ``tilt_x_h``, ``tilt_y_h`` (and V counterparts)
    are views into indices 0, 1, and 2.
    """

    def __init__(
        self,
        max_zernike_radial_number=6,
        mask_size_x_in_pixels=1024,
        mask_size_y_in_pixels=1024,
        aperture_diameter_in_m=7.5e-3,
        SLM_pixel_size_in_um=9.2e-6,
        coefs_path=None,
    ):
        """Create Zernike basis, load coefficients, and compute initial phase masks.

        Args:
            max_zernike_radial_number: Maximum radial order n for the basis.
            mask_size_x_in_pixels: Mask width in pixels.
            mask_size_y_in_pixels: Mask height in pixels.
            aperture_diameter_in_m: Aperture diameter in metres.
            SLM_pixel_size_in_um: Pixel size in metres (9.2 µm stored as 9.2e-6).
            coefs_path: Optional JSON coefficients file; zeros if omitted.
        """
        self.max_zernike_radial_number = max_zernike_radial_number
        self.mask_size_x_in_pixels = mask_size_x_in_pixels
        self.mask_size_y_in_pixels = mask_size_y_in_pixels
        self.aperture_diameter_in_m = aperture_diameter_in_m
        self.SLM_pixel_size_in_um = SLM_pixel_size_in_um

        half_width = mask_size_x_in_pixels / 2.0
        half_height = mask_size_y_in_pixels / 2.0

        x, y = np.meshgrid(
            np.linspace(0, mask_size_x_in_pixels - 1, mask_size_x_in_pixels, dtype=np.float32),
            np.linspace(0, mask_size_y_in_pixels - 1, mask_size_y_in_pixels, dtype=np.float32),
        )
        rho, phi = cart2pol(x - half_width, y - half_height)

        self.N = rho.shape[0]
        r_norm = aperture_diameter_in_m / SLM_pixel_size_in_um / 2
        rho = rho / r_norm

        self.zernike_basis_array = self.generate_zernike_basis(
            rho, phi, max_zernike_radial_number
        )
        self._init_buffers()
        self.coefs_path = None

        if coefs_path is not None:
            self.load_coefs(coefs_path)
        else:
            self.clean()

        self.write_backup_of_zernike_coefs()
        self.calculate_znk_phase_mask(restore_coefs=False)

    @staticmethod
    def nm_to_idx(n, m):
        """Basis index for mode (n, m). Weight lives at ``w_h[nm_to_idx(n, m)]``."""
        return catalog_nm_to_idx(n, m)

    @staticmethod
    def idx_to_nm(idx):
        """Return (n, m) for basis row ``idx`` (e.g. ``idx_to_nm(4)`` -> defocus)."""
        return catalog_idx_to_nm(idx)

    def _init_buffers(self):
        """Preallocate coefficient, scratch, and output phase buffers."""
        count = expected_zernike_count(self.max_zernike_radial_number)
        self.w_h = np.zeros(count, dtype=np.float32)
        self.w_v = np.zeros(count, dtype=np.float32)

        self.znk_phase_h = np.zeros((self.N, self.N), dtype=np.float32)
        self.znk_phase_v = np.zeros((self.N, self.N), dtype=np.float32)
        self._zern_phase_flat_h = np.zeros(self.N * self.N, dtype=np.float32)
        self._zern_phase_flat_v = np.zeros(self.N * self.N, dtype=np.float32)
        self._basis_flat = self.zernike_basis_array.reshape(self.zernike_count, -1)

    @property
    def piston_h(self):
        """H-pol piston (basis index 0). Alias into ``w_h[0]``."""
        return self.w_h[0]

    @piston_h.setter
    def piston_h(self, value):
        # piston_h = x  <=>  w_h[0] = x
        self.w_h[0] = np.float32(value)

    @property
    def tilt_x_h(self):
        """H-pol tilt X (basis index 1, Z_1^{-1}). Alias into ``w_h[1]``."""
        return self.w_h[1]

    @tilt_x_h.setter
    def tilt_x_h(self, value):
        # tilt_x_h = x  <=>  w_h[1] = x
        self.w_h[1] = np.float32(value)

    @property
    def tilt_y_h(self):
        """H-pol tilt Y (basis index 2, Z_1^{+1}). Alias into ``w_h[2]``."""
        return self.w_h[2]

    @tilt_y_h.setter
    def tilt_y_h(self, value):
        # tilt_y_h = x  <=>  w_h[2] = x
        self.w_h[2] = np.float32(value)

    @property
    def piston_v(self):
        """V-pol piston (basis index 0). Alias into ``w_v[0]``."""
        return self.w_v[0]

    @piston_v.setter
    def piston_v(self, value):
        self.w_v[0] = np.float32(value)

    @property
    def tilt_x_v(self):
        """V-pol tilt X (basis index 1). Alias into ``w_v[1]``."""
        return self.w_v[1]

    @tilt_x_v.setter
    def tilt_x_v(self, value):
        self.w_v[1] = np.float32(value)

    @property
    def tilt_y_v(self):
        """V-pol tilt Y (basis index 2). Alias into ``w_v[2]``."""
        return self.w_v[2]

    @tilt_y_v.setter
    def tilt_y_v(self, value):
        self.w_v[2] = np.float32(value)

    @staticmethod
    def _radial_coefficients(n, m_abs):
        coeffs = []
        for k in range((n - m_abs) // 2 + 1):
            c = ((-1) ** k) * factorial(n - k) / (
                factorial(k)
                * factorial((n + m_abs) // 2 - k)
                * factorial((n - m_abs) // 2 - k)
            )
            coeffs.append((np.float32(c), n - 2 * k))
        return coeffs

    def _build_rho_powers(self, rho, max_radial_idx):
        """Cache rho**k powers up to 2*max_radial_idx."""
        max_power = 2 * max_radial_idx
        rho_powers = [np.ones_like(rho, dtype=np.float32)]
        for power in range(1, max_power + 1):
            rho_powers.append(rho ** power)
        return rho_powers

    def _compute_radial_nm(self, rho, rho_powers, n, m_abs, radial_coeffs_cache):
        """Radial polynomial R_n^|m|(rho) on the instance grid."""
        r_nm = np.zeros((self.N, self.N), dtype=np.float32)
        key = (n, m_abs)
        if key not in radial_coeffs_cache:
            radial_coeffs_cache[key] = self._radial_coefficients(n, m_abs)
        for coeff, power in radial_coeffs_cache[key]:
            if power == 0:
                r_nm += coeff
            else:
                r_nm += coeff * rho_powers[power]
        r_nm = np.where(rho > 1, 0, r_nm)
        r_nm = np.where(rho == 1, 1, r_nm)
        return r_nm

    def _calculate_modelab_angular_convention(self, r_nm, phi, n, m, odd=True):
        """Angular part with Modelab sign convention."""
        if odd:
            if n % 2 == 0:
                z_nm = (-1) ** (m // 2 + 1) * r_nm * np.sin(m * phi)
            else:
                z_nm = (-1) ** ((m + 1) // 2) * r_nm * np.cos(m * phi)
        else:
            if (n % 2 == 0) and (m == 0):
                z_nm = r_nm * np.cos(m * phi) + (-1) ** (n // 2 + 1)
            elif (n % 2 == 0) and (m != 0):
                z_nm = (-1) ** (m // 2) * r_nm * np.cos(m * phi)
            else:
                z_nm = (-1) ** ((m + 1) // 2) * r_nm * np.sin(m * phi)
        return z_nm

    def generate_zernike_basis(self, rho, phi, max_radial_idx):
        """Precompute Modelab Zernike modes on the polar grid."""
        self.zernike_count = expected_zernike_count(max_radial_idx)
        zernike_basis_array = np.zeros((self.zernike_count, self.N, self.N), dtype=np.float32)

        rho_powers = self._build_rho_powers(rho, max_radial_idx)
        radial_coeffs_cache = {}
        radial_cache = {}

        for n in range(max_radial_idx + 1):
            for m in range(-n, n + 1, 2):
                m_abs = abs(m)
                radial_key = (n, m_abs)
                if radial_key not in radial_cache:
                    radial_cache[radial_key] = self._compute_radial_nm(
                        rho, rho_powers, n, m_abs, radial_coeffs_cache
                    )
                z_nm = self._calculate_modelab_angular_convention(
                    radial_cache[radial_key], phi, n, m_abs, odd=(m < 0)
                )
                idx = (n * (n + 2) + m) // 2
                zernike_basis_array[idx] = z_nm

        return zernike_basis_array

    def _apply_coefficients(self, coefficients):
        """Load mask centers and unified weight vectors from a coefficient dict."""
        self.h_mask_center_x = coefficients["h_mask_center_x"]
        self.h_mask_center_y = coefficients["h_mask_center_y"]
        self.v_mask_center_x = coefficients["v_mask_center_x"]
        self.v_mask_center_y = coefficients["v_mask_center_y"]
        self.w_h = np.asarray(coefficients["w_h"], dtype=np.float32)
        self.w_v = np.asarray(coefficients["w_v"], dtype=np.float32)

    def clean(self):
        """Reset all coefficients and mask centers to zero."""
        self._apply_coefficients(
            dict_to_zernike_coefficients(
                default_coefficients_dict(self.max_zernike_radial_number),
                max_zernike_radial_number=self.max_zernike_radial_number,
            )
        )

    def load_coefs(self, path):
        """Load coefficients from JSON (instance radial order wins on mismatch)."""
        data = load_coefficients_json(path)
        coefficients = dict_to_zernike_coefficients(
            data,
            max_zernike_radial_number=self.max_zernike_radial_number,
        )
        self._apply_coefficients(coefficients)
        self.coefs_path = str(Path(path).resolve())

    def save_coefs(self, path):
        """Save current coefficients to JSON."""
        data = zernike_coefficients_to_dict(self)
        save_coefficients_json(path, data)
        self.coefs_path = str(Path(path).resolve())

    def write_backup_of_zernike_coefs(self):
        """Snapshot current ``w_h`` / ``w_v`` as the base coefficient state."""
        self.b_w_h = np.copy(self.w_h)
        self.b_w_v = np.copy(self.w_v)

    def read_backup_of_zernike_coefs(self):
        """Restore ``w_h`` / ``w_v`` from the last backup snapshot."""
        self.w_h = np.copy(self.b_w_h)
        self.w_v = np.copy(self.b_w_v)

    def modify_zernike_coefs(
        self,
        apply_tilt_offset=False,
        apply_defocus_offset=False,
        scan_tilt=False,
        scan_defocus=False,
        **kwargs,
    ):
        """Apply temporary tilt/defocus offsets or scans (not used at init).

        Snapshots backup first. Follow with ``calculate_znk_phase_mask(restore_coefs=True)``.
        """
        self.write_backup_of_zernike_coefs()

        if apply_tilt_offset or apply_defocus_offset:
            if "MASK_corrections" not in kwargs:
                raise ValueError(
                    '"apply_tilt_offset" or "apply_defocus_offset" requires '
                    '"MASK_corrections" dictionary.'
                )
            mask_corrections = kwargs["MASK_corrections"]

            if apply_tilt_offset:
                self.tilt_x_h += mask_corrections["tilt_x_H"]
                self.tilt_y_h += mask_corrections["tilt_y_H"]
                self.tilt_x_v += mask_corrections["tilt_x_V"]
                self.tilt_y_v += mask_corrections["tilt_y_V"]

            if apply_defocus_offset:
                self.w_h[IDX_DEFOCUS] += mask_corrections["defocus_H"]
                self.w_v[IDX_DEFOCUS] += mask_corrections["defocus_V"]

        if scan_tilt:
            if not all(key in kwargs for key in ("tilt_x", "tilt_y")):
                raise ValueError('scan_tilt=True requires "tilt_x" and "tilt_y" keyword arguments.')
            self.tilt_x_h += kwargs["tilt_x"]
            self.tilt_y_h += kwargs["tilt_y"]
            self.tilt_x_v += kwargs["tilt_x"]
            self.tilt_y_v += kwargs["tilt_y"]

        if scan_defocus:
            if "defocus" not in kwargs:
                raise ValueError('scan_defocus=True requires "defocus" keyword argument.')
            self.w_h[IDX_DEFOCUS] += kwargs["defocus"]
            self.w_v[IDX_DEFOCUS] += kwargs["defocus"]

    def calculate_znk_phase_mask(self, restore_coefs=False):
        """Synthesize wrapped phase masks ``znk_phase_h`` and ``znk_phase_v``.

        Computes the phase as ``w @ basis`` (1D weights × 2D basis), then wraps
        to [-pi, pi] radians (same convention as ``pyLCOS._addMasks``).

        Args:
            restore_coefs: If True, restore ``w_h``/``w_v`` from backup after computing
                (use after ``modify_zernike_coefs`` scans).
        """
        np.matmul(self.w_h, self._basis_flat, out=self._zern_phase_flat_h)
        np.matmul(self.w_v, self._basis_flat, out=self._zern_phase_flat_v)

        pi = np.pi
        z_h = self._zern_phase_flat_h.reshape(self.N, self.N)
        z_v = self._zern_phase_flat_v.reshape(self.N, self.N)
        ne.evaluate(
            "(((pi * z + pi) % (2*pi)) - pi)",
            out=self.znk_phase_h,
            local_dict={"z": z_h, "pi": pi},
        )
        ne.evaluate(
            "(((pi * z + pi) % (2*pi)) - pi)",
            out=self.znk_phase_v,
            local_dict={"z": z_v, "pi": pi},
        )

        if restore_coefs:
            self.read_backup_of_zernike_coefs()

    def make_zernike_fields(self, restore_coefs=False):
        """Deprecated alias for ``calculate_znk_phase_mask``."""
        self.calculate_znk_phase_mask(restore_coefs=restore_coefs)

    def filter_coefficients_by_aberration_keys(self, pol: str, *aberration_keys: str ):
        """Filter coefficients by polarization and aberration keys. Return the coefficients and the indices of the filtered coefficients."""

        groups = iter_modes_by_group(self.max_zernike_radial_number)
        idx_list = []
        for aberration_key in aberration_keys:
            if aberration_key not in groups:
                raise KeyError(f"Unknown or empty aberration group: {aberration_key!r}, supported keys are: {ABERRATION_GROUPS}")
            idx_list.extend(idx for idx, _, _ in groups[aberration_key])

        idx_list = sorted(set(idx_list)) #eliminate duplicated and sorted them ascending

        if pol == "H":
            return self.w_h[idx_list], idx_list
        elif pol == "V":
            return self.w_v[idx_list], idx_list
        else:
            raise ValueError(f"Invalid polarization: {pol}")


    def generate_table(self, dataframe = True):
        """Print and return coefficient table: index, n, m, H, V, name.

        Uses pandas if available; otherwise prints aligned text and returns a list of dicts.
        """
        rows = []
        for idx, n, m in iter_modes(self.max_zernike_radial_number):
            rows.append(
                {
                    "index": idx,
                    "n": n,
                    "m": m,
                    "H": float(self.w_h[idx]),
                    "V": float(self.w_v[idx]),
                    "name": get_mode_name(n, m),
                }
            )

        try:
            if dataframe:
                import pandas as pd

                df = pd.DataFrame(rows)
                self.table_output = df.to_string(index=False)

                return df

        except ImportError:
            header = f"{'index':>5}  {'n':>2}  {'m':>3}  {'H':>10}  {'V':>10}  name"

            return header, rows
