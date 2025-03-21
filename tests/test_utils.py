import numpy as np
import pytest
import pandas as pd
from vaxsim.utils import (
    model_loss, 
    auc_below_threshold, 
    generate_seed_schedule
)

@pytest.fixture
def sample_time_series():
    """Fixture providing sample time series data."""
    days = 100
    # Create sample data with pandas DatetimeIndex using 'ME' for month end
    date_range = pd.date_range(start='2020-01-01', periods=3, freq='ME')
    return {
        'S': np.ones(days) * 600000,
        'I': np.zeros(days),
        'R': np.ones(days) * 200000,
        'V': np.ones(days) * 200000,
        'data': pd.DataFrame({
            'sero_eff': [0.4, 0.5, 0.6],
            'diva': [0.2, 0.25, 0.3]
        }, index=date_range)
    }

def test_model_loss(sample_time_series):
    """Test model loss calculation."""
    loss = model_loss(
        sample_time_series['S'],
        sample_time_series['I'],
        sample_time_series['R'],
        sample_time_series['V'],
        sample_time_series['data']
    )
    assert isinstance(loss, float)
    assert loss >= 0

def test_auc_below_threshold():
    """Test area under curve calculation."""
    days = 100
    S = np.ones(days) * 600000
    I = np.zeros(days)
    R = np.ones(days) * 200000
    V = np.ones(days) * 200000
    
    auc = auc_below_threshold(S, I, R, V, days)
    
    assert isinstance(auc, float)
    assert auc >= 0
    assert auc <= 100  # Percentage should be between 0 and 100

def test_seed_schedule_generation():
    """Test infection seeding schedule generation."""
    days = 100
    schedule = generate_seed_schedule(
        method='random',
        min_day=0,
        max_day=days,
        days=days,
        num_seeds=5
    )
    
    assert isinstance(schedule, list)
    assert len(schedule) == days
    assert sum(schedule) == 5  # Check number of seeding events
    assert all(x in [0, 1] for x in schedule)  # Check binary values