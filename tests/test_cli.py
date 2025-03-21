import pytest
from pathlib import Path
import tempfile
import yaml
import sys
import pandas as pd
import argparse
from vaxsim.cli import load_params, parse_seed_infection, main

@pytest.fixture
def mock_argv(monkeypatch):
    """Fixture to mock sys.argv."""
    def _mock_argv(args):
        monkeypatch.setattr(sys.argv, args)
    return _mock_argv

@pytest.fixture
def mock_data_file(tmp_path):
    """Create a mock data file for testing."""
    data = pd.DataFrame({
        'date': pd.date_range(start='2020-01-01', periods=10),
        'sero_eff': [0.4] * 10,
        'diva': [0.2] * 10
    })
    data_path = tmp_path / 'data copy.csv'
    data.to_csv(data_path, index=False)
    return data_path

@pytest.fixture
def sample_params():
    """Create a temporary params.yaml file for testing."""
    params = {
        'baseline': {
            'beta': 0.125,
            'gamma': 0.07,
            'vax_rate': 0.00833,
            'days': 10,  # Reduced for faster testing
            'S0': 639996,
            'I0': 4,
            'R0': 180000,
            'V0': 180000,
            'seed_rate': 0,
            'vax_period': 180,
            'vax_duration': 30,
            'start_vax_day': 30,
            'weibull_shape_vax': 3,
            'weibull_scale_vax': 220,
            'weibull_shape_rec': 3,
            'weibull_scale_rec': 1380,
        }
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        params_path = Path(tmpdir) / 'params.yaml'
        with params_path.open('w') as f:
            yaml.dump(params, f)
        yield params_path

def test_load_params(sample_params, monkeypatch):
    """Test parameter loading from different locations."""
    monkeypatch.chdir(sample_params.parent)
    params = load_params()
    assert 'baseline' in params
    assert params['baseline']['beta'] == 0.125

def test_parse_seed_infection():
    """Test parsing of seed infection arguments."""
    # Test valid inputs
    assert parse_seed_infection("random:10") == ("random", 10)
    assert parse_seed_infection("none:0") == ("none", 0)
    
    # Test invalid inputs with expected error format
    expected_error = "Seed infection format must be 'method:rate', e.g., 'random:10'"
    
    # Test invalid method
    with pytest.raises(argparse.ArgumentTypeError) as exc_info:
        parse_seed_infection("invalid:10")
    assert str(exc_info.value) == expected_error
    
    # Test invalid format
    with pytest.raises(argparse.ArgumentTypeError) as exc_info:
        parse_seed_infection("random:abc")
    assert str(exc_info.value) == expected_error

def test_cli_arguments(monkeypatch, capsys, mock_data_file, sample_params):
    """Test CLI argument parsing."""
    test_args = ['vaxsim', '--scenario', 'baseline', '--model_type', 'random']
    monkeypatch.setattr(sys, 'argv', test_args)
    monkeypatch.chdir(mock_data_file.parent)
    
    try:
        main()
    except SystemExit as e:
        assert e.code == 0

def test_cli_output_creation(tmp_path, monkeypatch, mock_data_file, sample_params):
    """Test if CLI creates necessary output directories."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    monkeypatch.chdir(mock_data_file.parent)
    
    test_args = ['vaxsim', '--scenario', 'baseline', '--model_type', 'random']
    monkeypatch.setattr(sys, 'argv', test_args)
    
    try:
        main()
        assert (output_dir / "logs").exists()
    except SystemExit:
        pass

def test_seed_infection_integration(monkeypatch, mock_data_file, sample_params):
    """Test seed infection integration."""
    test_cases = [
        ['vaxsim', '--scenario', 'baseline', '--seed_infection', 'random:5'],
        ['vaxsim', '--scenario', 'baseline', '--seed_infection', 'none:0']
    ]
    
    for args in test_cases:
        monkeypatch.setattr(sys, 'argv', args)
        monkeypatch.chdir(mock_data_file.parent)
        try:
            main()
        except SystemExit:
            pass