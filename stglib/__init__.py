from . import (
    _version,
    abss,
    aqd,
    argonaut,
    core,
    eco,
    eofe,
    exo,
    glx,
    hobo,
    indexvel,
    iq,
    lisst,
    mc,
    rdi,
    rsk,
    sg,
    sig,
    tb,
    tcm,
    troll,
    vec,
    wxt,
)
from .aqd import aqdutils
from .core import cmd, filter, qaqc, utils, waves
from .core.utils import read_globalatts
from .sg import sgutils

__version__ = _version.get_versions()["version"]
