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
    son,
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
from .son import sonutils

__version__ = _version.get_versions()["version"]

__cfmax__ = "CF-1.11"
