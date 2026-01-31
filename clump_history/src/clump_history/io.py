import pandas as pd
import os

def load_thermal_history(filepath, time_col=None, temp_col=None):
    """
    Load thermal history data from a CSV file.
    
    Parameters
    ----------
    filepath : str
        Path to the CSV file containing thermal history data.
    time_col : str, optional
        Column name for time values. If None, defaults to 'Time/Myr'.
    temp_col : str, optional
        Column name for temperature values. If None, defaults to 'Avg_T/Celsius'.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing time and temperature columns.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If required columns are not found in the file.
    """
    # Convert relative path to absolute path
    if not os.path.isabs(filepath):
        filepath = os.path.join(os.getcwd(), filepath)
    
    # Check if file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Read CSV file
    df = pd.read_csv(filepath)
    
    # Set default column names if not provided
    if time_col is None:
        time_col = 'Time/Myr'
    if temp_col is None:
        temp_col = 'Avg_T/Celsius'
    
    # Check if required columns exist
    if time_col not in df.columns:
        raise ValueError(f"Required column '{time_col}' not found in file.")
    if temp_col not in df.columns:
        raise ValueError(f"Required column '{temp_col}' not found in file.")
    
    # Return only the required columns
    return df[[time_col, temp_col]].copy()

def load_test_data(filepath, d47_col=None, sd_col=None):
    """
    Load test data from a CSV file.
    
    Parameters
    ----------
    filepath : str
        Path to the CSV file containing test data.
    d47_col : str, optional
        Column name for Δ47 values. If None, defaults to 'Delta47'.
    sd_col : str, optional
        Column name for standard deviation values. If None, defaults to 'SD'.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing Δ47 and SD columns.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If required columns are not found in the file.
    """
    # Convert relative path to absolute path
    if not os.path.isabs(filepath):
        filepath = os.path.join(os.getcwd(), filepath)
    
    # Check if file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Read CSV file
    df = pd.read_csv(filepath)
    
    # Set default column names if not provided
    if d47_col is None:
        d47_col = 'Delta47'
    if sd_col is None:
        sd_col = 'SD'
    
    # Check if required columns exist
    if d47_col not in df.columns:
        raise ValueError(f"Required column '{d47_col}' not found in file.")
    if sd_col not in df.columns:
        raise ValueError(f"Required column '{sd_col}' not found in file.")
    
    # Return only the required columns
    return df[[d47_col, sd_col]].copy()