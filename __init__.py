"""ModeLab Zernike basis and phase-mask synthesis for SLM holography."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("zernikes")
except PackageNotFoundError:
    __version__ = "1.0.0"

from .pyZernikes import Zernikes

__all__ = [
    "Zernikes",
    "__version__",
]
