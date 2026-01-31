"""
Unit tests for the ClumpIsotope model functions.

These tests verify the basic functionality of the model implementation
and ensure that changes haven't broken core functionality.
"""

import sys
import os
import numpy as np
import pytest

# Add the source directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'clump_history', 'src'))

from clump_history.model import compute_history


def test_compute_history_basic():
    """Test basic functionality of compute_history function."""
    # Simple thermal history: constant temperature over time
    time_myr = np.array([10, 5, 0])  # Going from past to present
    T_celsius = np.array([50, 50, 50])  # Constant temperature
    T_k = T_celsius + 273.15  # Convert to Kelvin
    ed = 120  # Typical activation energy
    d0_std = 0.05  # Typical initial std
    
    D, Dstd, Deq = compute_history(time_myr, T_k, ed, d0_std)
    
    # Check that all outputs have the same length as inputs
    assert len(D) == len(time_myr)
    assert len(Dstd) == len(time_myr)
    assert len(Deq) == len(time_myr)
    
    # Check that all outputs are numpy arrays
    assert isinstance(D, np.ndarray)
    assert isinstance(Dstd, np.ndarray)
    assert isinstance(Deq, np.ndarray)


def test_compute_history_different_temps():
    """Test compute_history with varying temperatures."""
    # Varying temperature history
    time_myr = np.array([10, 5, 0])
    T_celsius = np.array([100, 75, 50])  # Cooling over time
    T_k = T_celsius + 273.15
    ed = 120
    d0_std = 0.05
    
    D, Dstd, Deq = compute_history(time_myr, T_k, ed, d0_std)
    
    # Verify dimensions match
    assert len(D) == len(time_myr)
    assert len(Dstd) == len(time_myr)
    assert len(Deq) == len(time_myr)
    
    # With cooling history, equilibrium values should decrease toward present
    # (higher temperatures yield lower Delta47 values)
    assert Deq[0] < Deq[-1]  # Past (hotter) should have lower eq value than present (cooler)


def test_compute_history_edge_cases():
    """Test edge cases for compute_history function."""
    # Single time point
    time_myr = np.array([0])
    T_k = np.array([373.15])  # 100°C
    ed = 120
    d0_std = 0.05
    
    D, Dstd, Deq = compute_history(time_myr, T_k, ed, d0_std)
    
    assert len(D) == 1
    assert len(Dstd) == 1
    assert len(Deq) == 1


if __name__ == "__main__":
    # Run the tests
    test_compute_history_basic()
    test_compute_history_different_temps()
    test_compute_history_edge_cases()
    print("All tests passed!")