import unittest
import stglib
import xarray as xr
import numpy as np


class TestTimes(unittest.TestCase):

    def test_epic_time_conversion(self):
        bbvcf = xr.open_dataset(
            ('http://geoport.whoi.edu/thredds/dodsC/silt/usgs/Projects/'
             'stellwagen/CF-1.6/CHINCOTEAGUE/10191Aaqd-a.nc'))
        bbvepic = xr.open_dataset(
            ('http://stellwagen.er.usgs.gov/thredds/dodsC/TSdata/'
             'CHINCOTEAGUE/10191Aaqd-a.nc'), decode_times=False)

        difftime = bbvcf['time'].values - stglib.utils.epic_to_datetime(
           bbvepic['time'], bbvepic['time2']).values

        # values should be equal to at least 1 ms
        assert np.all(difftime <= np.timedelta64(1, 'ms'))


if __name__ == '__main__':
    unittest.main()
