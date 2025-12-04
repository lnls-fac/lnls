"""."""
import os as _os

from lnls.timer import Timer, TimerError

from . import dialog, ids, rotating_coil, utils

__all__ = ['utils', 'rotating_coil', 'ids', 'dialog']

with open(_os.path.join(__path__[0], 'VERSION'), 'r') as _f:
    __version__ = _f.read().strip()
