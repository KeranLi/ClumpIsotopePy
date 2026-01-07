import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def constrained_u_fit(time_myr, values, start_x, end_x, new_max_temp_k, plot=False):
    time_myr = np.asarray(time_myr, dtype=float)
    values = np.asarray(values, dtype=float)

    mask = (time_myr >= start_x) & (time_myr <= end_x)
    local_time = time_myr[mask]
    local_values = values[mask]

    if local_time.size < 3:
        raise ValueError(
            f"Peak window [{start_x}, {end_x}] has too few points ({local_time.size})."
        )

    start_value = local_values[0]
    end_value = local_values[-1]
    mid_x = (start_x + end_x) / 2.0
    start_slope = (local_values[1] - local_values[0]) / (local_time[1] - local_time[0])

    def poly(x, a, b, c, d, e):
        return a * x**4 + b * x**3 + c * x**2 + d * x + e

    def dpoly(x, a, b, c, d):
        return 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d

    def constraints(params):
        a, b, c, d, e = params
        return [
            poly(start_x, a, b, c, d, e) - start_value,
            dpoly(start_x, a, b, c, d) - start_slope,
            poly(end_x, a, b, c, d, e) - end_value,
            poly(mid_x, a, b, c, d, e) - new_max_temp_k,
            dpoly(mid_x, a, b, c, d),
        ]

    initial_guess = [0.0, 0.0, 0.0, float(start_slope), float(start_value)]
    a, b, c, d, e = fsolve(constraints, initial_guess)

    adjusted_local = poly(local_time, a, b, c, d, e)
    adjusted = values.copy()
    adjusted[mask] = adjusted_local

    if plot:
        plt.figure(figsize=(10, 5))
        plt.plot(time_myr, values, lw=1.5, label="Original")
        plt.plot(local_time, adjusted_local, lw=2, label="Adjusted")
        plt.axvline(start_x, ls="--", lw=1)
        plt.axvline(end_x, ls="--", lw=1)
        plt.scatter([mid_x], [new_max_temp_k], s=80, label="Peak")
        plt.xlabel("Time (Myr)")
        plt.ylabel("Temperature (K)")
        plt.legend()
        plt.tight_layout()
        plt.show()

    return adjusted
