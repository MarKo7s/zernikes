import json
from pathlib import Path

import numpy as np

COEFFICIENTS_VERSION = 1


def build_nm_lists(max_zernike_radial_number):
    n_list = []
    m_list = []
    for n_idx in range(2, max_zernike_radial_number + 1):
        for m_idx in range(-n_idx, n_idx + 1, 2):
            n_list.append(n_idx)
            m_list.append(m_idx)
    return n_list, m_list


def expected_w_length(max_zernike_radial_number):
    return sum(2 * n_idx + 1 for n_idx in range(2, max_zernike_radial_number + 1))


def default_coefficients_dict(max_zernike_radial_number):
    w_len = expected_w_length(max_zernike_radial_number)
    n_list, m_list = build_nm_lists(max_zernike_radial_number)
    zero_w = [0.0] * w_len

    return {
        "version": COEFFICIENTS_VERSION,
        "max_zernike_radial_number": max_zernike_radial_number,
        "mask_centers": {
            "H": {"x": 0, "y": 0},
            "V": {"x": 0, "y": 0},
        },
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


def zernike_coefficients_to_dict(zernikes):
    w_h = np.asarray(zernikes.w_h, dtype=np.float32)
    w_v = np.asarray(zernikes.w_v, dtype=np.float32)
    n_list, m_list = build_nm_lists(zernikes.max_zernike_radial_number)

    return {
        "version": COEFFICIENTS_VERSION,
        "max_zernike_radial_number": zernikes.max_zernike_radial_number,
        "mask_centers": {
            "H": {
                "x": int(zernikes.h_mask_center_x),
                "y": int(zernikes.h_mask_center_y),
            },
            "V": {
                "x": int(zernikes.v_mask_center_x),
                "y": int(zernikes.v_mask_center_y),
            },
        },
        "H": {
            "piston": float(zernikes.piston_h),
            "tilt_x": float(zernikes.tilt_x_h),
            "tilt_y": float(zernikes.tilt_y_h),
            "w": w_h.tolist(),
        },
        "V": {
            "piston": float(zernikes.piston_v),
            "tilt_x": float(zernikes.tilt_x_v),
            "tilt_y": float(zernikes.tilt_y_v),
            "w": w_v.tolist(),
        },
        "nm": [[n, m] for n, m in zip(n_list, m_list)],
    }


def _validate_coefficients_data(data, max_zernike_radial_number=None):
    if data.get("version") != COEFFICIENTS_VERSION:
        raise ValueError(
            f"Unsupported coefficients version: {data.get('version')!r}. "
            f"Expected {COEFFICIENTS_VERSION}."
        )

    file_max_n = data.get("max_zernike_radial_number")
    if file_max_n is None:
        raise ValueError("Coefficients file is missing 'max_zernike_radial_number'.")

    if max_zernike_radial_number is not None and file_max_n != max_zernike_radial_number:
        raise ValueError(
            f"Coefficients max_zernike_radial_number={file_max_n} does not match "
            f"Zernikes instance max_zernike_radial_number={max_zernike_radial_number}."
        )

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

    mask_centers = data.get("mask_centers")
    if mask_centers is None:
        raise ValueError("Coefficients file is missing 'mask_centers'.")
    for pol in ("H", "V"):
        center = mask_centers.get(pol)
        if center is None or "x" not in center or "y" not in center:
            raise ValueError(f"Coefficients file is missing 'mask_centers.{pol}'.")


def dict_to_zernike_coefficients(data, max_zernike_radial_number=None):
    _validate_coefficients_data(data, max_zernike_radial_number)

    file_max_n = data["max_zernike_radial_number"]
    n_list, m_list = build_nm_lists(file_max_n)

    if "nm" in data:
        file_nm = [tuple(pair) for pair in data["nm"]]
        expected_nm = list(zip(n_list, m_list))
        if file_nm != expected_nm:
            raise ValueError("Coefficients 'nm' list does not match expected Zernike ordering.")

    return {
        "h_mask_center_x": int(data["mask_centers"]["H"]["x"]),
        "h_mask_center_y": int(data["mask_centers"]["H"]["y"]),
        "v_mask_center_x": int(data["mask_centers"]["V"]["x"]),
        "v_mask_center_y": int(data["mask_centers"]["V"]["y"]),
        "piston_h": np.float32(data["H"]["piston"]),
        "tilt_x_h": np.float32(data["H"]["tilt_x"]),
        "tilt_y_h": np.float32(data["H"]["tilt_y"]),
        "piston_v": np.float32(data["V"]["piston"]),
        "tilt_x_v": np.float32(data["V"]["tilt_x"]),
        "tilt_y_v": np.float32(data["V"]["tilt_y"]),
        "w_h": np.asarray(data["H"]["w"], dtype=np.float32),
        "w_v": np.asarray(data["V"]["w"], dtype=np.float32),
        "n": n_list,
        "m": m_list,
        "zern_idxs": np.asarray(
            [(n * (n + 2) + m) // 2 for n, m in zip(n_list, m_list)],
            dtype=np.int64,
        ),
    }


def save_coefficients_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def load_coefficients_json(path):
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)
