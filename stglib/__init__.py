from . import (
    _version,
    aqd,
    argonaut,
    core,
    eco,
    eofe,
    exo,
    hobo,
    indexvel,
    iq,
    lisst,
    rdi,
    rsk,
    sig,
    troll,
    vec,
    wxt,
)
from .aqd import aqdutils
from .core import cmd, utils, waves
from .core.utils import read_globalatts

__version__ = _version.get_versions()["version"]
