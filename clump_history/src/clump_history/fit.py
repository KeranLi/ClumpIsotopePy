"""
Constrained U-Fit Adjustment Module for Thermal History Modeling

This module implements a constrained polynomial fitting algorithm to adjust 
thermal histories by introducing temperature peaks within specified time windows.
The U-fit approach maintains geological constraints while allowing for targeted
adjustments to thermal models based on independent temperature indicators.

The method fits a 4th-order polynomial to a specific time window to create a 
U-shaped temperature profile that reaches a specified maximum temperature 
while maintaining continuity with surrounding temperature values.

Authors:
    - Original implementation: Keran Li (Nanjing University)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def constrained_u_fit(time_myr, values, start_x, end_x, new_max_temp_k, plot=False):
    """
    Apply a constrained U-shaped adjustment to a thermal history within a specified time window.
    
    This function modifies a thermal history by fitting a 4th-order polynomial within 
    a specified time window to achieve a desired maximum temperature while maintaining 
    smooth transitions with adjacent segments. The resulting temperature profile has 
    a U-shape characteristic of localized heating events.
    
    The algorithm enforces five constraints:
    1. Polynomial passes through the starting temperature value
    2. Polynomial derivative matches the initial slope at the start
    3. Polynomial passes through the ending temperature value
    4. Polynomial reaches the specified maximum temperature at the midpoint
    5. Polynomial derivative is zero at the midpoint (ensuring a maximum)
    
    Parameters
    ----------
    time_myr : array-like
        Time values in million years (Ma). Conventionally, time increases into the past.
    values : array-like
        Temperature values corresponding to each time point in `time_myr`.
    start_x : float
        Start of the time window for adjustment [Myr].
    end_x : float
        End of the time window for adjustment [Myr].
    new_max_temp_k : float
        New maximum temperature to achieve at the midpoint of the window [Kelvin].
    plot : bool, optional
        Whether to display a diagnostic plot comparing original and adjusted profiles.
        
    Returns
    -------
    numpy.ndarray
        Adjusted temperature array with the same shape as the input `values`.
        
    Raises
    ------
    ValueError
        If the specified time window contains fewer than 3 points.
        
    Notes
    -----
    The U-fit adjustment is particularly useful for modeling localized heating events 
    such as igneous intrusions, hydrothermal circulation, or contact metamorphism 
    that affected the sample's thermal history and Δ47 composition.
    
    The method ensures continuity and smoothness at the boundaries of the adjustment
    window while allowing for substantial temperature changes within the window.
    
    Examples
    --------
    >>> import numpy as np
    >>> time = np.linspace(0, 100, 100)  # 0 to 100 Ma
    >>> temp = 10 + 0.1 * time  # Baseline thermal history
    >>> adjusted = constrained_u_fit(time, temp, 30, 50, 400, plot=True)
    """
    time_myr = np.asarray(time_myr, dtype=float)
    values = np.asarray(values, dtype=float)

    # Create a mask for the specified time window
    mask = (time_myr >= start_x) & (time_myr <= end_x)
    local_time = time_myr[mask]
    local_values = values[mask]

    # Validate that the window contains sufficient data points
    if local_time.size < 3:
        raise ValueError(
            f"Peak window [{start_x}, {end_x}] has too few points ({local_time.size}). "
            "At least 3 points are required for the polynomial fit."
        )

    # Extract boundary conditions for the polynomial
    start_value = local_values[0]
    end_value = local_values[-1]
    mid_x = (start_x + end_x) / 2.0
    # Estimate the initial slope from the first two points
    start_slope = (local_values[1] - local_values[0]) / (local_time[1] - local_time[0])

    # Define the 4th-order polynomial and its derivative
    def poly(x, a, b, c, d, e):
        """4th-order polynomial: ax^4 + bx^3 + cx^2 + dx + e"""
        return a * x**4 + b * x**3 + c * x**2 + d * x + e

    def dpoly(x, a, b, c, d):
        """First derivative of the polynomial: 4ax^3 + 3bx^2 + 2cx + d"""
        return 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d

    # Define constraint equations for the polynomial coefficients
    def constraints(params):
        """
        System of equations defining the constraints for the polynomial coefficients.
        
        The solution finds coefficients [a, b, c, d, e] that satisfy:
        1. Polynomial passes through start point
        2. Polynomial derivative matches initial slope at start
        3. Polynomial passes through end point
        4. Polynomial achieves specified maximum at midpoint
        5. Polynomial derivative is zero at midpoint (ensuring extremum)
        """
        a, b, c, d, e = params
        return [
            # Constraint 1: Polynomial passes through starting temperature
            poly(start_x, a, b, c, d, e) - start_value,
            # Constraint 2: Matches initial slope at start
            dpoly(start_x, a, b, c, d) - start_slope,
            # Constraint 3: Polynomial passes through ending temperature
            poly(end_x, a, b, c, d, e) - end_value,
            # Constraint 4: Achieves specified maximum temperature at midpoint
            poly(mid_x, a, b, c, d, e) - new_max_temp_k,
            # Constraint 5: Derivative is zero at midpoint (ensures extremum)
            dpoly(mid_x, a, b, c, d),
        ]

    # Initial guess for polynomial coefficients [a, b, c, d, e]
    # Based on linear approximation with the known slope and value at start
    initial_guess = [0.0, 0.0, 0.0, float(start_slope), float(start_value)]
    
    # Solve for the polynomial coefficients that satisfy all constraints
    a, b, c, d, e = fsolve(constraints, initial_guess)

    # Apply the polynomial to adjust temperatures within the window
    adjusted_local = poly(local_time, a, b, c, d, e)
    
    # Create a copy of the original values and apply adjustments only within the window
    adjusted = values.copy()
    adjusted[mask] = adjusted_local

    # Optionally display a diagnostic plot
    if plot:
        plt.figure(figsize=(10, 5))
        plt.plot(time_myr, values, lw=1.5, label="Original Thermal History")
        plt.plot(local_time, adjusted_local, lw=2, label="Adjusted Profile")
        plt.axvline(start_x, ls="--", lw=1, color='gray', label="Window Boundaries")
        plt.axvline(end_x, ls="--", lw=1, color='gray')
        plt.scatter([mid_x], [new_max_temp_k], s=80, label="Target Peak", zorder=5)
        plt.xlabel("Time (Myr)")
        plt.ylabel("Temperature (K)")
        plt.title("Constrained U-Fit Thermal History Adjustment")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return adjusted