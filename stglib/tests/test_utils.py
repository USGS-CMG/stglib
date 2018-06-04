import unittest
import stglib
import xarray as xr
import numpy as np
import pandas as pd


class TestUtils(unittest.TestCase):

    def test_rotate(self):
        expected = (-np.sqrt(2)/2, np.sqrt(2)/2)
        result = stglib.aqd.qaqc.rotate(0, 1, -45)
        np.testing.assert_almost_equal(expected, result)

        expected = (np.sqrt(2)/2, np.sqrt(2)/2)
        result = stglib.aqd.qaqc.rotate(0, 1, 45)
        np.testing.assert_almost_equal(expected, result)
        print(expected, result)

        expected = (1, 0)
        result = stglib.aqd.qaqc.rotate(0, 1, 90)
        np.testing.assert_almost_equal(expected, result)

        expected = (np.sqrt(2)/2, -np.sqrt(2)/2)
        result = stglib.aqd.qaqc.rotate(0, 1, 135)
        np.testing.assert_almost_equal(expected, result)

        expected = (0, -1)
        result = stglib.aqd.qaqc.rotate(0, 1, 180)
        np.testing.assert_almost_equal(expected, result)

        expected = (-1, 0)
        result = stglib.aqd.qaqc.rotate(0, 1, 270)
        np.testing.assert_almost_equal(expected, result)
