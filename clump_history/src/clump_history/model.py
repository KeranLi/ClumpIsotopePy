"""
Carbonate Clumped Isotope Reordering Model

This module implements theoretical models for the temperature-dependent reordering 
of clumped isotopes in geological samples over time. The models simulate how 
carbonate clumped isotope compositions (Δ47) change due to solid-state isotope 
exchange reactions during burial and thermal maturation.

The implementation is based on:
1. Stolper, D. A., Eiler, J. M. (2015). The kinetics of solid-state 
   isotope-exchange reactions for clumped isotopes: A study of inorganic 
   calcites and apatites from natural and experimental samples. 
   American Journal of Science, 315(5), 363-401.

2. Hemingway, J. D., & Henkes, G. A. (2021). A disordered kinetic model 
   for clumped isotope bond reordering in carbonates. 
   Earth and Planetary Science Letters, 575, 117177.

Authors:
    - Original implementation: Keran Li (Nanjing University)
    - Theoretical basis: Clumped isotope geochemistry community
"""

import numpy as np
import isotopylog as ipl

# Conversion factor from million years to seconds
# 1 Myr = 1e6 years * 365 days/year * 24 hours/day * 60 minutes/hour * 60 seconds/minute
MYR_TO_SEC = 1e6 * 365 * 24 * 60 * 60


def compute_history(time_myr, T_k, ed, d0_std):
    """
    Compute the evolution of clumped isotope composition (Δ47) over a given thermal history.
    
    This function simulates the temperature-dependent reordering of clumped isotopes 
    in a geological sample subjected to a specific thermal history. It uses isotopylog 
    to perform the geologic history integration based on solid-state isotope exchange kinetics.
    
    Parameters
    ----------
    time_myr : array-like
        Time values in million years (Ma). Conventionally, time increases into the past,
        so the oldest time point comes first in the series.
    T_k : array-like
        Temperature values in Kelvin corresponding to each time segment. 
        Should match the time points in `time_myr`.
    ed : float
        Activation energy parameter for the isotope exchange reaction [kJ/mol].
        This represents the energy barrier for solid-state reordering.
    d0_std : float
        Initial standard deviation of the clumped isotope signal, representing 
        the initial state of the sample before significant reordering occurred.
        
    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing three arrays:
        - D : Array of reordered Δ47 values through time (reverse chronological order)
        - Dstd : Array of Δ47 standard deviations through time (reverse chronological order)
        - Deq : Array of equilibrium Δ47 values calculated from temperature history
        
    Notes
    -----
    The algorithm:
    1. Converts time from Myr to seconds
    2. Reverses the time and temperature series to work forward in time
    3. Computes initial equilibrium Δ47 value from the first temperature
    4. Integrates the reordering equations using isotopylog
    5. Reverses the results back to the original time direction
    
    The reordering process depends on both temperature and duration of exposure,
    with higher temperatures causing more rapid equilibration toward the 
    temperature-dependent equilibrium Δ47 value.
    
    References
    ----------
    .. [1] Stolper, D. A., & Eiler, J. M. (2015). American Journal of Science, 315(5), 363-401.
    .. [2] Hemingway, J. D., & Henkes, G. A. (2021). Earth and Planetary Science Letters, 575, 117177.
    """
    time_myr = np.asarray(time_myr, dtype=float)
    T_k = np.asarray(T_k, dtype=float)

    # Reverse the time and temperature arrays to process from present to past
    # This is necessary because the integration proceeds forward in time
    time_rev = time_myr[::-1]
    T_rev = T_k[::-1]

    # Convert time from Myr to seconds
    t_sec = time_rev * MYR_TO_SEC
    # Adjust time so that t=0 corresponds to the first time point
    t_sec = t_sec[0] - t_sec

    # Calculate the initial equilibrium Δ47 value from the first temperature
    D0 = ipl.Deq_from_T(T_rev[0])
    # Initialize the starting conditions [Δ47 value, d47/dt, other params]
    d0 = [D0, 0, 0]

    # Perform the geologic history integration using isotopylog
    # This computes how Δ47 evolves over the thermal history
    D_rev, Dstd_rev = ipl.geologic_history(t_sec, T_rev, ed, d0, d0_std=[d0_std])

    # Reverse the results back to the original time direction (past to present)
    D = D_rev[::-1]
    Dstd = Dstd_rev[::-1]
    # Calculate equilibrium Δ47 values for comparison with actual values
    Deq = ipl.Deq_from_T(T_k)

    return D, Dstd, Deq