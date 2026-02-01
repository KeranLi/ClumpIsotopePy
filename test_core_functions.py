#!/usr/bin/env python
"""
Test core functions without launching GUI
"""
import sys
sys.path.insert(0, 'clump_history/src')

import pandas as pd
import numpy as np
import isotopylog as ipl
from clump_history.io import load_thermal_history, load_test_data
from clump_history.model import compute_history
from clump_history.fit import constrained_u_fit

print("=" * 60)
print("Testing Clump History Core Functions")
print("=" * 60)

# Test 1: Load data
print("\n[1/5] Testing data loading...")
try:
    thermal_df = pd.read_csv('datasets/Thermal_History_Hu.csv')
    time_myr = thermal_df['Time/Myr'].values
    avg_temp_c = thermal_df['Avg_T/Celsius'].values
    avg_temp_k = avg_temp_c + 273.15
    print(f"  - Loaded thermal history: {len(time_myr)} points")
    print(f"  - Time range: {time_myr.min()} - {time_myr.max()} Myr")
    print(f"  - Temp range: {avg_temp_c.min()} - {avg_temp_c.max()} C")
    
    test_df = pd.read_csv('datasets/acutal_test_Hu.csv')
    delta47 = test_df['Delta47'].values
    delta47_err = test_df['SD'].values
    print(f"  - Loaded test data: {len(delta47)} samples")
    print("  PASSED!")
except Exception as e:
    print(f"  FAILED: {e}")
    sys.exit(1)

# Test 2: EDistribution
print("\n[2/5] Testing EDistribution loading...")
try:
    ed = ipl.EDistribution.from_literature(
        mineral='calcite',
        reference='HH21'
    )
    print("  - Loaded HH21 calcite EDistribution")
    print("  PASSED!")
except Exception as e:
    print(f"  FAILED: {e}")
    sys.exit(1)

# Test 3: U-Fit function
print("\n[3/5] Testing U-Fit adjustment...")
try:
    adjusted_temp = constrained_u_fit(
        time_myr, avg_temp_k, 
        start_x=550, end_x=600, 
        new_max_temp_k=150+273.15, 
        plot=False
    )
    print(f"  - Adjusted temperature range: {adjusted_temp.min()-273.15:.1f} - {adjusted_temp.max()-273.15:.1f} C")
    print("  PASSED!")
except Exception as e:
    print(f"  FAILED: {e}")
    sys.exit(1)

# Test 4: Compute history (single scenario)
print("\n[4/5] Testing forward modeling...")
try:
    D, Dstd, Deq = compute_history(time_myr, avg_temp_k, ed, d0_std=0.02)
    print(f"  - Computed D47 range: {D.min():.3f} - {D.max():.3f}")
    print(f"  - Final D47: {D[-1]:.3f} +/- {Dstd[-1]:.3f}")
    print("  PASSED!")
except Exception as e:
    print(f"  FAILED: {e}")
    sys.exit(1)

# Test 5: Multiple scenarios
print("\n[5/5] Testing multiple peak temperatures...")
try:
    peak_temps = [150, 200, 250]
    scenarios = []
    
    # Initial scenario
    D, Dstd, Deq = compute_history(time_myr, avg_temp_k, ed, d0_std=0.02)
    scenarios.append(('initial', D, Dstd, Deq))
    print(f"  - Initial scenario: D47 final = {D[-1]:.3f}")
    
    # Peak scenarios
    for temp_c in peak_temps:
        temp_k = temp_c + 273.15
        adjusted = constrained_u_fit(time_myr, avg_temp_k, 550, 600, temp_k, plot=False)
        D, Dstd, Deq = compute_history(time_myr, adjusted, ed, d0_std=0.02)
        scenarios.append((f'{temp_c}C', D, Dstd, Deq))
        print(f"  - Peak {temp_c}C scenario: D47 final = {D[-1]:.3f}")
    
    print("  PASSED!")
except Exception as e:
    print(f"  FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 60)
print("All tests passed! You can now run the GUI.")
print("=" * 60)
