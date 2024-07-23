# tests/test_find_twin_obs.py

from MTWSPy.find_twin import FindTwinObs
find_twin_obs = FindTwinObs()

from MTWSPy.v01_phasenames import Phases

from MTWSPy.toolkit import Toolkit
tk = Toolkit()

from obspy import read
import os
import numpy as np

def test_get_tt_times():

    pf = 'params_in.yaml'
    params = tk.get_params(pf)

    lf_name = 'test_output/temp_lf.log'
    of_name = 'test_output/temp_of.log'

    logfile = open(lf_name,'w')
    outfile = open(of_name,'w')

    input_dict={}
    input_dict['params_in'] = params
    input_dict['event_id'] = '20080101XXXXXX'
    input_dict['phases'] = Phases().get_phase_dictionary()
    filename = './MTWSPy/data/obs/e2008/py_formatted/20080101063232/20080101063232_IU_AFI.00.LHT'
    seis = read(filename)

    fail = 0

    input_dict, filename, seis, logfile, outfile, fail = find_twin_obs.get_tt_times(input_dict, filename, seis, logfile, outfile, fail)

    assert 'ttt' in input_dict
    assert seis[0].stats['station'] == 'AFI'
    assert len(seis[0].stats['traveltimes']) == 15
    assert np.round(seis[0].stats['traveltimes']['SS'],1) == 2176.7

    logfile.close()
    outfile.close()
    os.remove(lf_name)
    os.remove(of_name)



def main():
    test_get_tt_times()

    return


if __name__ == '__main__':
    main()


