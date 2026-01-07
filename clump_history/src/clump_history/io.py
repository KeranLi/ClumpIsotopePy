from pathlib import Path
import numpy as np
import pandas as pd


def load_thermal_history(csv_path: Path, time_col: str, avg_col: str):
    df = pd.read_csv(csv_path)
    if time_col not in df.columns:
        raise ValueError(f"Thermal history CSV missing column: {time_col}")
    if avg_col not in df.columns:
        raise ValueError(f"Thermal history CSV missing column: {avg_col}")

    time_myr = df[time_col].to_numpy(dtype=float)
    T_k = df[avg_col].to_numpy(dtype=float) + 273.15

    # 保证 time 升序
    if time_myr[0] > time_myr[-1]:
        time_myr = time_myr[::-1]
        T_k = T_k[::-1]

    return time_myr, T_k


def load_test_data(csv_path: Path, d47_col: str, sd_col: str):
    df = pd.read_csv(csv_path)
    if d47_col not in df.columns:
        raise ValueError(f"Test CSV missing column: {d47_col}")
    if sd_col not in df.columns:
        raise ValueError(f"Test CSV missing column: {sd_col}")

    delta47 = df[d47_col].to_numpy(dtype=float)
    delta47_err = df[sd_col].to_numpy(dtype=float)
    return delta47, delta47_err
