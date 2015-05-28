import os as _os
from lnls.system import *
from lnls.timer import *
from . import ids

with open(_os.path.join(__path__[0], 'VERSION'), 'r') as _f:
    __version__ = _f.read().strip()
