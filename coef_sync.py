"""On-demand comparison of current Zernike weights against saved copies."""

from dataclasses import dataclass
from datetime import datetime

import numpy as np


@dataclass
class CoefSync:
    """Timestamps for the JSON mirror and the RAM backup.

    Weight arrays stay on the Zernikes instance. This object only remembers
    when each copy was taken and compares them when asked.
    """

    json_saved_at: datetime | None = None
    backup_at: datetime | None = None

    def status(self, w_h, w_v, json_w_h, json_w_v, b_w_h, b_w_v):
        """Return JSON and backup sync against the current weight vectors."""
        return {
            "json": {
                "timestamp": self.json_saved_at,
                "sync": _vectors_match(w_h, w_v, json_w_h, json_w_v),
            },
            "backup": {
                "timestamp": self.backup_at,
                "sync": _vectors_match(w_h, w_v, b_w_h, b_w_v),
            },
        }


def _vectors_match(w_h, w_v, other_h, other_v):
    """True if both polarizations match, False if not, None if the copy is missing."""
    if other_h is None or other_v is None:
        return None
    return bool(np.array_equal(w_h, other_h) and np.array_equal(w_v, other_v))
