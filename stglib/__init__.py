from . import aqd, argonaut, core, eco, exo, hobo, indexvel, iq, rdi, rsk, troll
from ._version import get_versions
from .core import cmd, utils, waves
from .core.utils import read_globalatts
from .aqd import qaqc

__version__ = get_versions()["version"]
del get_versions
