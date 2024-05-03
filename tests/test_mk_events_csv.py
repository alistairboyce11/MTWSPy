# tests/test_mk_events_csv.py

import MTWSPy.mk_events_csv as mk_events_csv
import MTWSPy.toolkit as toolkit
import os

def test_get_catalogs():

    pf = 'params_in.yaml'
    params = toolkit.get_params(pf)

    logfile = open('temp.log','w')

    cmt_cat = mk_events_csv.get_cmt_catalog(params, logfile)
    assert len(cmt_cat) > 0

    dmt_cat = mk_events_csv.get_dmt_catalog_year(params, logfile)
    logfile.close()
    os.remove('temp.log')

    assert len(dmt_cat) > 0
