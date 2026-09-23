"""JSON coefficient I/O and assembly for the Zernikes module."""

import json
import warnings
from pathlib import Path

import numpy as np

try:
    from .zernike_catalog import (
        basis_index,
        expected_zernike_count,
        iter_modes,
    )
except ImportError:
    from zernike_catalog import (
        basis_index,
        expected_zernike_count,
        iter_modes,
    )

COEFFICIENTS_SCHEMA = "zernikes.coefficients"
COEFFICIENTS_VERSION = 1


def zernike_coeffs_default_path():
    """Last-saved coefficient file inside the package config folder."""
    return Path(__file__).resolve().parent / "config" / "zernike_coeffs.json"


def _preamble(saved_at):
    """File identity written with every coefficient JSON."""
    return {
        "schema": COEFFICIENTS_SCHEMA,
        "version": COEFFICIENTS_VERSION,
        "saved_at": saved_at,
    }


def build_nm_lists(max_zernike_radial_number):
    """Build parallel n, m lists for higher-order modes (n >= 2)."""
    n_list = []
    m_list = []
    for n_idx in range(2, max_zernike_radial_number + 1):
        for m_idx in range(-n_idx, n_idx + 1, 2):
            n_list.append(n_idx)
            m_list.append(m_idx)
    return n_list, m_list


def expected_w_length(max_zernike_radial_number):
    """Length of compact higher-order w vector in JSON (modes with n >= 2)."""
    return sum(n_idx + 1 for n_idx in range(2, max_zernike_radial_number + 1))


def default_coefficients_dict(max_zernike_radial_number):
    """Default zero coefficient dictionary for JSON export."""
    w_len = expected_w_length(max_zernike_radial_number)
    n_list, m_list = build_nm_lists(max_zernike_radial_number)
    zero_w = [0.0] * w_len

    return {
        "preamble": _preamble(None),
        "max_zernike_radial_number": max_zernike_radial_number,
        "H": {
            "piston": 0.0,
            "tilt_x": 0.0,
            "tilt_y": 0.0,
            "w": zero_w,
        },
        "V": {
            "piston": 0.0,
            "tilt_x": 0.0,
            "tilt_y": 0.0,
            "w": zero_w,
        },
        "nm": [[n, m] for n, m in zip(n_list, m_list)],
    }


def _validate_coefficients_data(data):
    """Validate JSON structure against the file's own max radial order."""
    preamble = data.get("preamble")
    schema = preamble.get("schema") if isinstance(preamble, dict) else None
    if schema != COEFFICIENTS_SCHEMA:
        raise ValueError(
            f"Coefficients file schema is {schema!r}. Expected {COEFFICIENTS_SCHEMA!r}."
        )
    version = preamble.get("version")
    if version != COEFFICIENTS_VERSION:
        raise ValueError(
            f"Unsupported coefficients version: {version!r}. Expected {COEFFICIENTS_VERSION}."
        )

    file_max_n = data.get("max_zernike_radial_number")
    if file_max_n is None:
        raise ValueError("Coefficients file is missing 'max_zernike_radial_number'.")

    expected_len = expected_w_length(file_max_n)
    for pol in ("H", "V"):
        pol_data = data.get(pol)
        if pol_data is None:
            raise ValueError(f"Coefficients file is missing '{pol}' section.")
        for key in ("piston", "tilt_x", "tilt_y", "w"):
            if key not in pol_data:
                raise ValueError(f"Coefficients file is missing '{pol}.{key}'.")
        if len(pol_data["w"]) != expected_len:
            raise ValueError(
                f"Expected {expected_len} higher-order coefficients for '{pol}.w', "
                f"got {len(pol_data['w'])}."
            )


def _file_pol_lookup(data, pol):
    """Build (n, m) -> weight lookup from one polarization section of a JSON file."""
    pol_data = data[pol]
    file_max_n = data["max_zernike_radial_number"]
    n_list, m_list = build_nm_lists(file_max_n)

    if "nm" in data:
        file_nm = [tuple(pair) for pair in data["nm"]]
        expected_nm = list(zip(n_list, m_list))
        if file_nm != expected_nm:
            raise ValueError("Coefficients 'nm' list does not match expected Zernike ordering.")

    lookup = {
        (0, 0): float(pol_data["piston"]),
        (1, -1): float(pol_data["tilt_x"]),
        (1, 1): float(pol_data["tilt_y"]),
    }
    for (n, m), weight in zip(zip(n_list, m_list), pol_data["w"]):
        lookup[(n, m)] = float(weight)
    return lookup


def _warn_order_mismatch(file_max_n, instance_max_n):
    """Warn when file radial order differs from the Zernikes instance."""
    if file_max_n == instance_max_n:
        return
    if file_max_n > instance_max_n:
        warnings.warn(
            f"Coefficients file max_zernike_radial_number={file_max_n} is greater than "
            f"instance max_zernike_radial_number={instance_max_n}; "
            f"higher-order file coefficients are ignored.",
            stacklevel=3,
        )
    else:
        warnings.warn(
            f"Coefficients file max_zernike_radial_number={file_max_n} is less than "
            f"instance max_zernike_radial_number={instance_max_n}; "
            f"missing modes default to zero.",
            stacklevel=3,
        )


def _assemble_weight_vector(pol_lookup, instance_max_n, file_max_n):
    """Map file coefficients onto a unified instance-sized weight vector."""
    count = expected_zernike_count(instance_max_n)
    weights = np.zeros(count, dtype=np.float32)
    ignored = 0

    for idx, n, m in iter_modes(instance_max_n):
        if (n, m) in pol_lookup:
            weights[idx] = np.float32(pol_lookup[(n, m)])
        elif n <= file_max_n:
            pass
        else:
            ignored += 1

    file_modes = set(pol_lookup.keys()) - {(0, 0), (1, -1), (1, 1)}
    for n, m in file_modes:
        if n > instance_max_n:
            ignored += 1

    if ignored > 0 and file_max_n > instance_max_n:
        warnings.warn(
            f"Ignored {ignored} coefficient(s) from file with n > {instance_max_n}.",
            stacklevel=3,
        )

    return weights


def dict_to_zernike_coefficients(data, max_zernike_radial_number=None):
    """Convert JSON dict to internal coefficient arrays for a Zernikes instance."""
    _validate_coefficients_data(data)

    if max_zernike_radial_number is None:
        raise ValueError("max_zernike_radial_number is required when loading coefficients.")

    file_max_n = data["max_zernike_radial_number"]
    _warn_order_mismatch(file_max_n, max_zernike_radial_number)

    lookup_h = _file_pol_lookup(data, "H")
    lookup_v = _file_pol_lookup(data, "V")

    return {
        "w_h": _assemble_weight_vector(lookup_h, max_zernike_radial_number, file_max_n),
        "w_v": _assemble_weight_vector(lookup_v, max_zernike_radial_number, file_max_n),
    }


def zernike_coefficients_to_dict(zernikes, saved_at):
    """Serialize Zernikes instance coefficients to a JSON-ready dictionary."""
    max_n = zernikes.max_zernike_radial_number
    n_list, m_list = build_nm_lists(max_n)
    compact_w_h = []
    compact_w_v = []
    for n, m in zip(n_list, m_list):
        idx = basis_index(n, m)
        compact_w_h.append(float(zernikes.w_h[idx]))
        compact_w_v.append(float(zernikes.w_v[idx]))

    return {
        "preamble": _preamble(saved_at),
        "max_zernike_radial_number": max_n,
        "H": {
            "piston": float(zernikes.w_h[0]),
            "tilt_x": float(zernikes.w_h[1]),
            "tilt_y": float(zernikes.w_h[2]),
            "w": compact_w_h,
        },
        "V": {
            "piston": float(zernikes.w_v[0]),
            "tilt_x": float(zernikes.w_v[1]),
            "tilt_y": float(zernikes.w_v[2]),
            "w": compact_w_v,
        },
        "nm": [[n, m] for n, m in zip(n_list, m_list)],
    }


def save_coefficients_json(path, data):
    """Write coefficient dictionary to a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def load_coefficients_json(path):
    """Load coefficient dictionary from a JSON file."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)
