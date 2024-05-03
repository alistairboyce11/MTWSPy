# tests/test_match_catalog.py

from MTWSPy.match_catalog import MatchCatalog
from MTWSPy.toolkit import Toolkit
tk = Toolkit()
import os

def test_get_catalogs():

    pf = 'params_in.yaml'
    params = tk.get_params(pf)

    mc = MatchCatalog(params)

    logfile = open('temp.log','w')

    cmt_cat = mc.get_cmt_catalog()
    assert len(cmt_cat) > 0

    data_cat = mc.get_data_catalog_year()
    logfile.close()
    os.remove('temp.log')

    assert len(data_cat) > 0
