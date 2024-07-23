# tests/test_toolkit.py

from MTWSPy.toolkit import Toolkit
tk = Toolkit()

def test_check_files():

    loc = './MTWSPy/data/obs/e2008/py_formatted'
    year = '2008'
    component = 'T'
    numfiles = tk.check_files(loc, year, component)

    assert numfiles == 2109


def test_get_params():

    pf = 'params_in.yaml'

    params = tk.get_params(pf)

    assert 'obs' in params['obs_loc']
    assert 'syn' in params['syn_loc']


def main():
    test_check_files()
    test_get_params()

    return


if __name__ == '__main__':
    main()


