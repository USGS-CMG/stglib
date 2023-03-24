import dolfyn
import pandas as pd
from dolfyn.adv import api

from ..aqd import aqdutils
from ..core import utils


def raw_to_cdf(metadata):
    basefile = metadata["basefile"]

    if "prefix" in metadata:
        prefix = metadata["prefix"]
    else:
        prefix = ""
    basefile = prefix + basefile
    print(basefile)

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata, inst_type="SIG")

    ds = dolfyn.read(f"{basefile}.ad2cp")

    ds = utils.write_metadata(ds, metadata)
    ds = utils.ensure_cf(ds)

    cdf_filename = prefix + ds.attrs["filename"] + "-raw.nc"

    dolfyn.save(ds, cdf_filename, compression=True)
    print(f"Finished writing data to {cdf_filename}")
