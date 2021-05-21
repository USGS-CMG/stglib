import pandas as pd

from ..core import utils

def csv_to_cdf(metadata):

    basefile = metadata["basefile"]

    df = pd.read_csv(basefile + '.txt', infer_datetime_format=True)

    df = df.rename(columns={'Time': 'time'}).set_index('time')

    ds = df.to_xarray()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds['time'] = pd.DatetimeIndex(ds["time"])

    if "Turbidity" in ds:
        ds = ds.rename({'Turbidity': 'Turb'})

        ds["Turb"].attrs.update(
            {"units": "Nephelometric turbidity units (NTU)", "long_name": "Turbidity"}
        )

    ds = utils.shift_time(ds, 0)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds
