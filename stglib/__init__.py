from . import eco
from . import exo
from . import iq
from . import hobo
from . import aqd
from . import rsk
from . import indexvel
from . import core
from . import troll
from . import argonaut
from .core import utils, cmd, waves
from .core.utils import read_globalatts

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
