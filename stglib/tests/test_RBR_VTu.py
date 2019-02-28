import unittest
import stglib
import yaml

class Test_RBR_VTu(unittest.TestCase):
    
    def setUp(self):
        filepath = '../../examples/RBRTu_glob_att.txt'
        gatts = stglib.utils.read_globalatts(filepath)
        # these are the most important metadata
        self.assertEqual(gatts['MOORING'],'1110')
        self.assertEqual(gatts['WATER_DEPTH'],8.6)
        self.assertEqual(gatts['latitude'],29.7114)
        self.assertEqual(gatts['longitude'],-81.218633)
        
        filepath = '../../examples/RBRTu_config.yaml'
        with open(filepath) as f:
            config = yaml.safe_load(f)
        self.assertEqual(config['basefile'],'RBRTu_rawdata')
        self.assertEqual(config['instrument_type'],'Virtuoso Tu')
        
        for k in config:
            gatts[k] = config[k]
        
        self.metadata = gatts
        
    def test_rsk_to_cdf(self):
            
        # need to adjust metadata for the fact that the data file is stored
        # over in examples, like the other examples
        self.metadata['basefile'] = '../../examples/'+self.metadata['basefile']
        self.ds = stglib.rsk.rsk2cdf.rsk_to_cdf(self.metadata)
        self.assertEqual(float(self.ds['Turb'][int(len(self.ds['Turb'])/2)]),21.066973271097595)
        
        # note that this warning occurs even when the tests are passed,
        # ~\AppData\Local\Continuum\miniconda3\envs\IOOS\lib\importlib\_bootstrap.py:219: 
        # ImportWarning: can't resolve package from __spec__ or __package__, falling back on __name__ and __path__

    def test_cdf_to_nc(self):

        cdffilename = self.metadata['filename']+'-raw.cdf'
        
        self.ds = stglib.rsk.cdf2nc.cdf_to_nc(cdffilename)
        self.assertEqual(float(self.ds['Turb'][int(len(self.ds['Turb'])/2)]),72.64784091035335)        
        
if __name__ == '__main__':
    unittest.main()
