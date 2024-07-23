# tests/test_match_twin_files.py

def test_find_common_events():


    from MTWSPy.toolkit import Toolkit
    tk = Toolkit()

    pf = 'params_in.yaml'
    params = tk.get_params(pf)

    from MTWSPy.match_twin_files import MatchTwinFiles
    match_twin_files = MatchTwinFiles(params)




def main():
    test_find_common_events()

    return


if __name__ == '__main__':
    main()


