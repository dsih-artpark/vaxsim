import numpy as np
import pytest
from vaxsim.model import sirsv_model_with_weibull_random_vaccination

@pytest.fixture
def baseline_params():
    """Fixture providing baseline model parameters."""
    return {
        'beta': 0.125,
        'gamma': 0.07,
        'vax_rate': 0.00833,
        'weibull_shape_vax': 3,
        'weibull_scale_vax': 220,
        'weibull_shape_rec': 3,
        'weibull_scale_rec': 1380,
        'days': 100,
        'S0': 639996,
        'I0': 4,
        'R0': 180000,
        'V0': 180000,
        'seed_rate': 0,
        'vax_period': 180,
        'vax_duration': 30,
        'start_vax_day': 30,
    }

def test_population_conservation(baseline_params):
    """Test if total population remains constant throughout simulation."""
    S, I, R, V = sirsv_model_with_weibull_random_vaccination(baseline_params, 'test')
    N0 = sum([baseline_params[k] for k in ['S0', 'I0', 'R0', 'V0']])
    
    for t in range(len(S)):
        assert np.isclose(S[t] + I[t] + R[t] + V[t], N0)

def test_initial_conditions(baseline_params):
    """Test if model correctly sets initial conditions."""
    S, I, R, V = sirsv_model_with_weibull_random_vaccination(baseline_params, 'test')
    
    assert S[0] == baseline_params['S0']
    assert I[0] == baseline_params['I0']
    assert R[0] == baseline_params['R0']
    assert V[0] == baseline_params['V0']

def test_negative_populations(baseline_params):
    """Test that populations never become negative."""
    S, I, R, V = sirsv_model_with_weibull_random_vaccination(baseline_params, 'test')
    
    assert all(S >= 0)
    assert all(I >= 0)
    assert all(R >= 0)
    assert all(V >= 0)