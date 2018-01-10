from __future__ import print_function, division
import pandas as pd
import numpy as np
import scipy.stats

def parse_qrev_xml(doc, negateq=False):
    adcp = {}
    lendoc = len(doc['Channel']['Transect'])
    adcp['starttime'] = pd.to_datetime([doc['Channel']['Transect'][n]['StartDateTime']['#text'] for n in range(lendoc)])
    adcp['endtime'] = pd.to_datetime([doc['Channel']['Transect'][n]['EndDateTime']['#text'] for n in range(lendoc)])
    adcp['time'] = pd.to_datetime(np.mean([adcp['starttime'].view('i8'), adcp['endtime'].view('i8')], axis=0).astype('datetime64[ns]'))
    adcp['q'] = np.asarray([float(doc['Channel']['Transect'][n]['Discharge']['Total']['#text']) for n in range(lendoc)])
    if negateq:
        adcp['q'] = -adcp['q']
    adcp['AreaQrev'] = np.asarray([float(doc['Channel']['Transect'][n]['Other']['Area']['#text']) for n in range(lendoc)])
    adcp['filename'] = np.asarray([doc['Channel']['Transect'][n]['Filename']['#text'] for n in range(lendoc)])

    return adcp

def linregress(adcp):
    adcp['slope'], adcp['intercept'], adcp['r_value'], adcp['p_value'], adcp['std_err'] = scipy.stats.linregress(adcp['veli'], adcp['Vca'])

    return adcp
