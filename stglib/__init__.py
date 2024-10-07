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
    mc,
    rdi,
    rsk,
    sg,
    sig,
    tcm,
    troll,
    vec,
    wxt,
)
from .aqd import aqdutils
from .core import cmd, utils, waves
from .core.utils import read_globalatts
from .sg import sgutils

__version__ = _version.get_versions()["version"]
