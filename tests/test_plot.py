import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from vaxsim.plot import (
    plot_model,
    plot_histogram,
    plot_waning,
    plot_parameter_sweep
)

@pytest.fixture
def sample_time_series():
    """Create sample time series data for testing."""
    days = 100
    return {
        'S': np.ones(days) * 600000,
        'I': np.zeros(days),
        'R': np.ones(days) * 200000,
        'V': np.ones(days) * 200000,
        'days': days
    }

@pytest.fixture
def mock_data_file(tmp_path):
    """Create a mock data CSV file."""
    data = pd.DataFrame({
        'date': pd.date_range(start='2020-01-01', periods=10, freq='ME'),
        'sero_eff': [0.4] * 10,
        'diva': [0.2] * 10
    })
    data_path = tmp_path / 'data copy.csv'
    data.to_csv(data_path)
    return data_path

@pytest.fixture
def output_dir(tmp_path):
    """Create temporary output directory."""
    out_dir = tmp_path / 'output' / 'plots'
    out_dir.mkdir(parents=True)
    return out_dir

def test_plot_model(sample_time_series, mock_data_file, output_dir, monkeypatch):
    """Test model plotting function."""
    # Change working directory to where mock data is
    monkeypatch.chdir(mock_data_file.parent)
    
    plot_model(
        sample_time_series['S'],
        sample_time_series['I'],
        sample_time_series['R'],
        sample_time_series['V'],
        sample_time_series['days'],
        'test_scenario',
        'random',
        output_dir=output_dir
    )
    
    expected_file = output_dir / 'test_scenario_plot_random.png'
    assert expected_file.exists()

def test_plot_histogram(output_dir):
    """Test histogram plotting function."""
    decay_times_vax = np.random.weibull(3, 1000) * 220
    decay_times_rec = np.random.weibull(3, 1000) * 1380
    
    plot_histogram(
        decay_times_vax,
        decay_times_rec,
        'test_scenario',
        round_counter=1,
        start=True
    )
    
    expected_path = Path('output/diagnosis/test_scenario')
    expected_file = expected_path / 'decay_times_test_scenario_round_1_begin.png'
    assert expected_file.exists()

def test_plot_waning(sample_time_series, output_dir):
    """Test waning immunity plot function."""
    plot_waning(
        sample_time_series['S'],
        sample_time_series['I'],
        sample_time_series['R'],
        sample_time_series['V'],
        sample_time_series['days'],
        'test_scenario',
        'random',
        output_dir=output_dir
    )
    
    expected_file = output_dir / 'test_scenario_waning_random.png'
    assert expected_file.exists()

def test_parameter_sweep(tmp_path, monkeypatch):
    """Test parameter sweep plotting function."""
    # Create output directory structure
    sweep_dir = tmp_path / "output" / "sweep"
    sweep_dir.mkdir(parents=True, exist_ok=True)
    
    # Create sample results
    results = [
        {
            'param1': p1, 
            'param2': p2, 
            'protected': np.random.random()
        }
        for p1 in [1, 2, 3]
        for p2 in [4, 5, 6]
    ]
    
    # Temporarily change working directory for test
    with monkeypatch.context() as m:
        m.chdir(tmp_path)
        
        plot_parameter_sweep(
            results,
            'param1',
            'param2',
            output_variable='protected',
            model_type='random'
        )
        
        # Check if file was created
        expected_file = sweep_dir / "parameter_sweep_param1_param2_protected_1_random.png"
        assert expected_file.exists()

@pytest.fixture(autouse=True)
def cleanup():
    """Clean up created files after tests."""
    plt.close('all')  # Close all figures
    yield
    # Clean up output directories after tests
    import shutil
    for path in ['output/plots', 'output/diagnosis', 'output/sweep']:
        if Path(path).exists():
            shutil.rmtree(path)