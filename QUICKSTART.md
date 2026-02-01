# ClumpIsotope Quick Start Guide

This guide will help you quickly get started with the ClumpIsotope software for carbonate clumped isotope (Δ47) forward/reordering calculations.

## Prerequisites

- Python 3.7 or higher (Python 3.9/3.10 recommended)
- Git (optional, for cloning the repository)

## Installation

### Option 1: Direct Installation from Repository

1. Clone the repository:
   ```bash
   git clone https://github.com/keran-li/ClumpIsotope.git
   cd ClumpIsotope
   ```

2. Navigate to the clump_history directory and install in development mode:
   ```bash
   cd clump_history
   pip install -e .
   ```

### Option 2: Using Conda Environment

1. Clone the repository:
   ```bash
   git clone https://github.com/keran-li/ClumpIsotope.git
   cd ClumpIsotope
   ```

2. Create and activate a conda environment:
   ```bash
   conda env create -f env-run.yaml
   conda activate clump-isotope
   ```

3. Install the package in development mode:
   ```bash
   cd clump_history
   pip install -e .
   ```

## Running the Applications

### Command Line Interface (CLI)

After installation, you can run:

```bash
# Check the version
clump-history --version

# Get help for the run command
clump-history run -h

# Get help for the ufit command
clump-history ufit -h

# Run a forward modeling scenario
clump-history run --thermal datasets/Thermal_History_Hu.csv --test datasets/acutal_test_Hu.csv --outdir results --out fig_smoke
```

### Graphical User Interface (GUI)

Launch the GUI with:
```bash
clump-history-gui
```

## Standalone Executables

If you prefer not to install Python dependencies, you can use the pre-built executables:

1. Download the latest release from the GitHub releases page
2. Extract the archive
3. Run either:
   - `clump-history-gui.exe` for the graphical interface
   - `clump-history.exe` for the command-line interface

### Building Your Own Executables

To build the executables yourself, run:
```bash
# From the project root
build_apps.bat
```

The executables will be created in the `clump_history/dist/` directory:
- `clump-history-gui/clump-history-gui.exe` - GUI application
- `clump-history/clump-history.exe` - CLI application

## Input Data Format

### Thermal History File

Your thermal history CSV file must include:
- `Time/Myr` - Time in millions of years
- `Avg_T/Celsius` - Average temperature in Celsius

Example:
```
Time/Myr,Avg_T/Celsius
0,25
10,30
20,35
...
```

### Test Data File (Optional)

For comparison with measured values:
- `Delta47` - Measured clumped isotope values
- `SD` - Standard deviation of measurements

Example:
```
Delta47,SD
0.450,0.010
0.445,0.008
...
```

## GUI Usage

The GUI has two main tabs:

### Run Tab
- Configure thermal history and test data files
- Set output directory and file name
- Specify column names if different from defaults
- Choose mineral type and reference model
- Define peak window and temperatures
- Adjust plotting options
- Click "Run" to execute the analysis

### Ufit Tab
- Load thermal history file
- Set output directory and CSV name
- Define peak window (start/end times)
- Specify peak temperature
- Click "Ufit" to apply the adjustment

## Theory Overview

The software implements two main models for Δ47 reordering:

1. **Stolper et al. (2015)**: Solid-state isotope-exchange reactions
2. **Hemingway & Henkes (2021)**: Disordered kinetic model for clumped isotope bond reordering

For more details on the theory, see [THEORY.md](THEORY.md).

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError**: Make sure you've installed the package with `pip install -e .`
2. **Permission errors**: Run as administrator or ensure you have write permissions to the directory
3. **GUI doesn't start**: Check that Tkinter is available with your Python installation

### Getting Help

- Check the [README.md](README.md) for detailed instructions
- Review the [THEORY.md](THEORY.md) for scientific background
- Open an issue on the GitHub repository if you encounter problems

## Next Steps

1. Try the example datasets included with the repository
2. Experiment with different thermal histories
3. Adjust parameters to see their effect on Δ47 values
4. Consult the scientific literature for guidance on parameter choices for your samples