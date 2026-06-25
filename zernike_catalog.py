"""Static Zernike mode names and index mapping for the SLM Zernikes module."""


#! Some Constants:
MAX_NAMED_ORDER = 4

IDX_PISTON = 0
IDX_TILT_X = 1
IDX_TILT_Y = 2

MODE_NAMES = {
    (0, 0): "Piston",
    (1, -1): "Tilt X",
    (1, 1): "Tilt Y",
    (2, -2): "Astigmatism 45°",
    (2, 0): "Defocus",
    (2, 2): "Astigmatism 0°/90°",
    (3, -3): "Trefoil Y",
    (3, -1): "Coma Y",
    (3, 1): "Coma X",
    (3, 3): "Trefoil X",
    (4, -4): "Tetrafoil Y",
    (4, -2): "Secondary astigmatism 45°",
    (4, 0): "Spherical",
    (4, 2): "Secondary astigmatism 0°/90°",
    (4, 4): "Tetrafoil X",
}

ABERRATION_GROUPS = (
    ("piston", "Piston"),
    ("tilt", "Tilt"),
    ("defocus", "Defocus"),
    ("spherical", "Spherical aberration"),
    ("astigmatism", "Astigmatism"),
    ("coma", "Coma"),
    ("foil", "Foils"),
)


#! Functions:
def _validate_nm(n: int, m: int) -> None:
    if n < 0:
        raise ValueError(f"Radial order n must be >= 0, got {n}.")
    if abs(m) > n:
        raise ValueError(f"|m| must be <= n for Zernike mode (n={n}, m={m}).")
    if (n - m) % 2 != 0:
        raise ValueError(f"m must have the same parity as n (n={n}, m={m}).")


def basis_index(n: int, m: int) -> int:
    """Internal basis row index for Zernike mode (n, m)."""
    _validate_nm(n, m)
    return (n * (n + 2) + m) // 2


def nm_to_idx(n: int, m: int) -> int:
    """Alias for basis_index(n, m)."""
    return basis_index(n, m)


IDX_DEFOCUS = basis_index(2, 0)


def aberration_group(n: int, m: int) -> str:
    """Classify a Zernike mode into an aberration family for UI grouping."""
    if n == 0:
        return "piston"
    if n == 1:
        return "tilt"
    m_abs = abs(m)
    if m_abs == 0:
        return "defocus" if n == 2 else "spherical"
    if m_abs == 1:
        return "coma"
    if m_abs == 2:
        return "astigmatism"
    return "foil"


def iter_modes_by_group(max_zernike_radial_number: int) -> dict[str, list[tuple[int, int, int]]]:
    """Map aberration group key to ``(index, n, m)`` lists (non-empty groups only)."""
    groups: dict[str, list[tuple[int, int, int]]] = {key: [] for key, _ in ABERRATION_GROUPS}
    for idx, n, m in iter_modes(max_zernike_radial_number):
        groups[aberration_group(n, m)].append((idx, n, m))
    return {key: modes for key, modes in groups.items() if modes}


def idx_to_nm(idx: int) -> tuple[int, int]:
    """Return (n, m) for a basis row index."""
    if idx < 0:
        raise ValueError(f"Index must be >= 0, got {idx}.")
    n = 0
    while True:
        first = (n * (n + 2) - n) // 2
        last = (n * (n + 2) + n) // 2
        if first <= idx <= last:
            m = 2 * idx - n * (n + 2)
            return n, m
        if n > idx:
            break
        n += 1
    raise ValueError(f"Index {idx} is not a valid Zernike basis index.")


def get_mode_name(n: int, m: int) -> str:
    """Human-readable mode name through order 4; empty string for n > 4."""
    if n > MAX_NAMED_ORDER:
        return ""
    return MODE_NAMES.get((n, m), "")


def iter_modes(max_zernike_radial_number: int):
    """Yield (index, n, m) in basis generation order for radial orders 0..max."""
    for n in range(max_zernike_radial_number + 1):
        for m in range(-n, n + 1, 2):
            yield basis_index(n, m), n, m


def expected_zernike_count(max_zernike_radial_number: int) -> int:
    """Number of basis rows allocated for a given max radial order."""
    return (max_zernike_radial_number + 1) ** 2
