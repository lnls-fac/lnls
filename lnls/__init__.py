import os as _os
from lnls.system import *
from lnls.timer import *
from . import ids
from . import utils
from . import rotating_coil
from . import dialog

__all__ = ['utils', 'rotating_coil', 'ids', 'dialog']

with open(_os.path.join(__path__[0], 'VERSION'), 'r') as _f:
    __version__ = _f.read().strip()
