# tests/test_v01_phasenames.py

from MTWSPy.v01_phasenames import Phases

def test_phases():

    phases = Phases().get_phase_dictionary()

    assert 'Master' in phases


