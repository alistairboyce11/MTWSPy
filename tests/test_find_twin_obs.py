# tests/test_find_twin_obs.py

import MTWSPy.find_twin_obs as find_twin_obs
import MTWSPy.v01_phasenames as v01_phasenames
import MTWSPy.toolkit as toolkit
from obspy import read
import os
import numpy as np

def test_get_tt_times():

    pf = 'params_in.yaml'
    params = toolkit.get_params(pf)

    lf_name = 'temp_lf.log'
    of_name = 'temp_of.log'

    logfile = open(lf_name,'w')
    outfile = open(of_name,'w')

    input_dict={}
    input_dict['params_in'] = params
    input_dict['event_id'] = '20080101XXXXXX'
    input_dict['phases'] = v01_phasenames.phases()
    filename = './MTWSPy/data/obs/e2008/py_formatted/20080101063232/20080101063232_II_AAK.00.LHT'
    seis = read(filename)

    fail = 0

    input_dict, filename, seis, logfile, outfile, fail = find_twin_obs.get_tt_times(input_dict, filename, seis, logfile, outfile, fail)

    assert 'ttt' in input_dict
    assert seis[0].stats['station'] == 'AAK'
    assert len(seis[0].stats['traveltimes']) == 15
    assert np.round(seis[0].stats['traveltimes']['S'],1) == 72.4

    logfile.close()
    outfile.close()
    os.remove(lf_name)
    os.remove(of_name)




