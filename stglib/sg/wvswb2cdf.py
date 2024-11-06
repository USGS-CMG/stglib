from ..core import utils
from . import sgutils


def wb_to_cdf(metadata):
    """
    Load a raw .wb and .hex file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    # Get metadata from .hex file
    hexmeta = sgutils.read_hex(basefile + ".hex")

    # Append to metadata variable
    metadata.update(hexmeta)

    # Read in data
    ds = sgutils.read_wb(basefile + ".wb")

    # Convert pressure from psia to dbar
    ds["P_1"] = ds.P_1 / 14.503773800722 * 10

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    cdf_filename = ds.attrs["filename"] + "-waves-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds
