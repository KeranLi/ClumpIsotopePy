## The carbonate clumped isotope calulcation by Python

#### Author: Keran Li (Nanjing University)

#### This repository is licensed under the ***MIT License***.

### Introduction

This repository contains code for carbonate clumped isotope (Δ47) forward/reordering calculations in Python, and includes a command-line interface (CLI) for running scenarios and exporting publication-quality figures.

### Dependencies

Core runtime (CLI):
- Python (tested with 3.7+; recommended 3.9/3.10 for long-term support)
- Numpy
- Scipy
- Matplotlib
- Pandas
- isotopylog
- Tkinter (usually bundled with Python on Windows)

Optional (development / notebook):
- Jupyter Notebook / JupyterLab

> If you use conda, you can recreate the environment from `environment.yaml` (or `env-run.yml` if you maintain a minimal runtime environment).

### What has been done?

1. ~~Realize paper "Hemingway, J. D., and Henkes, G. A disordered kinetic model for clumped isotope bond reordering in carbonates, 2021, EPSL."~~
2. ~~郭炀锐, 邓文峰, 韦刚健. 碳酸盐成岩作用中的团簇同位素地球化学研究进展. 2022, 矿物岩石地球化学通报.~~
3. ~~The carbonate clumped isotope reordering calculation by Python (Exchange/diffusion model from Stolper et al., 2015 | "Stolper, D. A., Eiler, J. M., The kinetics of solid-state isotope-exchange reactions for clumped isotopes: A study of inorganic calcites and apatites from natural and experimental samples. 2015. American Journal of Science").~~
4. ~~Added several examples for forward and backward reordering.~~
5. ~~Added a CLI interface (`clump-history`) for running Δ47 forward models and thermal-history peak adjustments.~~
6. ~~Packed the workflow into an .exe file.~~
7. ~~Added a GUI interface (`clump-history-gui`) wrapping the same workflow.~~

### Quickstart (CLI)

#### (1) Install (editable/development mode)
From the repository root (the folder containing `pyproject.toml`):

```bash
pip install -e .
```

You can then run either:

```bash
clump-history --version
clump-history run -h
clump-history ufit -h

clump-history-gui
```

You can also run via module:

```bash
python -m clump_history --version
```

#### (2) Run forward modeling scenarios and save figures

If your workspace layout contains `./datasets/` (recommended), run from the workspace root:

```bash
clump-history run --thermal datasets/Thermal_History_Hu.csv --test datasets/acutal_test_Hu.csv --outdir results --out fig_smoke
```

Outputs:

- `results/fig_smoke.pdf`
- `results/fig_smoke.svg`

#### (3) Apply U-fit peak adjustment and export an adjusted thermal history

```bash
clump-history ufit --thermal datasets/Thermal_History_Hu.csv --peak-window 550 600 --peak-temp 150 --outdir results --out-csv thermal_adjusted.csv
```

Output:

- `results/thermal_adjusted.csv`

> Note: If you run the CLI from a different directory, please pass explicit paths to `--thermal` / `--test`.

### GUI Usage

Start the GUI:

```bash
clump-history-gui
```

The GUI provides two tabs:

- **Run (2×3 scenarios)**: selects Thermal/Test CSV, sets peak window/peak temps, and exports PDF/SVG.
- **Ufit (peak adjust)**: applies constrained U-shaped peak adjustment and exports an adjusted CSV.

---

### Input data format

Thermal history CSV (example: `Thermal_History_Hu.csv`) must include:

- `Time/Myr`
- `Avg_T/Celsius`

Actual test CSV (example: `acutal_test_Hu.csv`) must include:

- `Delta47`
- `SD`

(You can change these column names using CLI flags such as `--time-col`, `--avg-col`, `--d47-col`, and `--sd-col`.)

## PyInstaller (Windows packaging)

### GUI (recommended, no console window)

```bash
pyinstaller run_gui.py ^
  --name clump-history-gui ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --collect-all matplotlib ^
  --collect-all isotopylog ^
  --collect-submodules scipy
```

### CLI

```bash
pyinstaller run_cli.py ^
  --name clump-history ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --console ^
  --collect-all matplotlib ^
  --collect-all isotopylog ^
  --collect-submodules scipy
```

### Results

1. Single initial Δ47 input  
    (1) ours mimic (time bar is converted):
    <div align="center">
    <img width="500" alt="image" src="https://user-images.githubusercontent.com/66153455/280744777-30c49076-eb05-42c7-b12f-0e7e169b6bd9.png">
    </div>

    (2) 刘鑫, 邱楠生, 冯乾乾. 碳酸盐岩团簇同位素约束下的川东地区二叠系热演化. 2023, 地质学报. results:
    <div align="center">
    <img width="500"  alt="image" src="https://user-images.githubusercontent.com/66153455/280745312-d4974462-7839-429e-9a60-23fb31db4722.png">
    </div>

    (3) ours Tibet samples:
    <div align="center">
    <img width="800" alt="image" src="https://private-user-images.githubusercontent.com/66153455/532728226-b5ecb601-fc7f-4410-a9c2-0f2a83f78851.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Njc3NzQ1ODksIm5iZiI6MTc2Nzc3NDI4OSwicGF0aCI6Ii82NjE1MzQ1NS81MzI3MjgyMjYtYjVlY2I2MDEtZmM3Zi00NDEwLWE5YzItMGYyYTgzZjc4ODUxLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNjAxMDclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjYwMTA3VDA4MjQ0OVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWE5YzJjYzRlNTk5MjE0MDBjN2NlZjk2ZGViNjgxZjNmMjRmY2E3ZWMxNzU0MTgyYzFhODYzYTQ1NjM0Yzc0OWImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.LXczP0EWcJoEsGfkSoeBQ-2xPy7BD8qDurH8Fh4Jk70">
    </div>

    (4) ours Guizhou samples:
    <div align="center">
    <img width="800" alt="image" src="https://private-user-images.githubusercontent.com/66153455/532728773-15c6028b-a2be-4f7d-9b64-6d637246af08.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Njc3NzQ2NzcsIm5iZiI6MTc2Nzc3NDM3NywicGF0aCI6Ii82NjE1MzQ1NS81MzI3Mjg3NzMtMTVjNjAyOGItYTJiZS00ZjdkLTliNjQtNmQ2MzcyNDZhZjA4LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNjAxMDclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjYwMTA3VDA4MjYxN1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTM4OWM5YzBhMjQ1YjgwMzE4MTZiNTJkZDJiNjBlM2Q4OGQ3ZWFlZWIxN2E5ODVlYjdkYjhmYWY4NTlmODY5NTgmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.milYsqDn-2XhaD7y-L0YrmxgCs8mkWHYulNC6ZWN2nI">
    </div>

    (5) ours Sichuan samples:
    <div align="center">
    <img width="500" alt="image" src="https://private-user-images.githubusercontent.com/66153455/532729071-835c96f5-f44e-4a94-ad5e-89e6b6ab9de5.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Njc3NzQ3MjYsIm5iZiI6MTc2Nzc3NDQyNiwicGF0aCI6Ii82NjE1MzQ1NS81MzI3MjkwNzEtODM1Yzk2ZjUtZjQ0ZS00YTk0LWFkNWUtODllNmI2YWI5ZGU1LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNjAxMDclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjYwMTA3VDA4MjcwNlomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTg4MThhZjMwMzlmYmQ5NGY5ODQzNmJlYjA4ODI5NDIwNjlhODg3MTM2Y2E3MTAzNGIzODA3ZTg2YjdkMjM3MWMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.DNhTjq9N8UVRCSg1J488kDs5725qu6BnGpWPwYWWYGs">
    </div>

### To do

1. Add Monte-Carlo simulation for different Δ47 inputs (todo)
2. Export scenario results (D, Dstd, Deq) to CSV for each run (todo)
