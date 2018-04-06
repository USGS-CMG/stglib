import unittest
import xarray as xr
import numpy as np

class TestStg(unittest.TestCase):

    def test_xarray_values(self):
        a = np.array([0, 1, 2, 3, 4])
        b = xr.DataArray([0, 1, 2, 3, 4]).values

        np.testing.assert_equal(a, b)

if __name__ == '__main__':
    unittest.main()
