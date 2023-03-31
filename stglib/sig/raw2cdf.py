import time

import dolfyn

from ..aqd import aqdutils
from ..core import utils


def raw_to_cdf(metadata):
    basefile = metadata["basefile"]

    if "prefix" in metadata:
        prefix = metadata["prefix"]
    else:
        prefix = ""
    basefile = prefix + basefile

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata, inst_type="SIG")

    start_time = time.time()
    ds = dolfyn.read(f"{basefile}.ad2cp")
    end_time = time.time()
    print(f"Finished loading {basefile}.ad2cp in {end_time-start_time:.1f} seconds")

    ds = utils.write_metadata(ds, metadata)
    ds = utils.ensure_cf(ds)

    cdf_filename = prefix + ds.attrs["filename"] + "-raw.nc"

    dolfyn.save(ds, cdf_filename, compression=True)
    print(f"Finished writing data to {cdf_filename}")
