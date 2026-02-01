## Carbonate Clumped Isotope Geochemistry Modeling Framework

#### Author: Keran Li (Nanjing University)

#### This repository is licensed under the ***MIT License***.

### Introduction

This repository contains a comprehensive Python-based computational framework for carbonate clumped isotope (Δ47) forward/reordering calculations. The software implements solid-state isotope exchange kinetic models to simulate the temperature-dependent reordering of clumped isotopes in geological samples over time, enabling reconstruction of thermal histories and paleotemperature estimates. 

The framework incorporates both theoretical models from Stolper et al. (2015) and the disordered kinetic model of Hemingway & Henkes (2021), providing researchers with multiple approaches for interpreting Δ47 values in geological materials.

### Scientific Background

Carbonate clumped isotope thermometry (Δ47) measures the abundance of multiply substituted isotopologues (e.g., 13C18O16O22-) in carbonates, providing a temperature proxy independent of the isotopic composition of the precipitating fluid. However, post-depositional heating can cause reordering of Δ47 values, potentially erasing primary formation temperatures and recording instead the effects of burial and thermal history. This computational framework quantitatively models these reordering processes.

### Dependencies

Core runtime (CLI):
- Python (tested with 3.7+; recommended 3.9/3.10 for long-term support)
- Numpy
- Scipy
- Matplotlib
- Pandas
- isotopylog
- Tkinter (usually bundled with Python on Windows)

Development & Testing:
- pytest
- jupyter notebook/lab

Optional (development / notebook):
- Jupyter Notebook / JupyterLab

> If you use conda, you can recreate the environment from `env-dev.yaml` (full development) or `env-run.yaml` (minimal runtime).

### Theoretical Models Implemented

1. **Stolper et al. (2015) Model**: Implementation of solid-state isotope-exchange reactions for clumped isotopes based on experimental data from natural and laboratory-reheated samples.
   - Reference: Stolper, D. A., Eiler, J. M. (2015). The kinetics of solid-state isotope-exchange reactions for clumped isotopes: A study of inorganic calcites and apatites from natural and experimental samples. *American Journal of Science*, 315(5), 363-401.

2. **Hemingway & Henkes (2021) Model**: Disordered kinetic model for clumped isotope bond reordering in carbonates with improved treatment of diffusion and disorder effects.
   - Reference: Hemingway, J. D., & Henkes, G. A. (2021). A disordered kinetic model for clumped isotope bond reordering in carbonates. *Earth and Planetary Science Letters*, 575, 117177.

### What has been done?

1. ~~Implement paper "Hemingway, J. D., and Henkes, G. A disordered kinetic model for clumped isotope bond reordering in carbonates, 2021, EPSL."~~
2. ~~Implement Guo Yangrui, Deng Wenfeng, Wei Gangjian. Clumped isotope geochemistry of carbonate diagenesis: Advances in research. 2022, *Bulletin of Mineralogy, Petrology and Geochemistry*.~~
3. ~~Implement the carbonate clumped isotope reordering calculation by Python (Exchange/diffusion model from Stolper et al., 2015 | "Stolper, D. A., Eiler, J. M., The kinetics of solid-state isotope-exchange reactions for clumped isotopes: A study of inorganic calcites and apatites from natural and experimental samples. 2015. American Journal of Science").~~
4. ~~Add examples for forward and backward reordering.~~
5. ~~Add a CLI interface (`clump-history`) for running Δ47 forward models and thermal-history peak adjustments.~~
6. ~~Package the workflow into executable applications.~~
7. ~~Add a GUI interface (`clump-history-gui`) wrapping the same workflow.~~

### Installation & Usage

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

Or from the project root:

```bash
cd clump_history
python run_gui.py
```

The GUI provides an intuitive interface for running forward modeling scenarios:

**Main Features:**

1. **数据文件选择** - 选择热历史CSV和实测Δ47数据CSV文件
2. **CSV列名设置** - 配置时间、温度、Δ47和误差列的名称
3. **模型参数** - 选择矿物类型(calcite/dolomite)和参考文献
4. **潜在热历史设置** - 
   - 设置U-Fit时间窗口 (开始-结束 Myr)
   - 输入多个峰值温度进行批量测试 (空格分隔，如: `150 200 250 300 350`)
5. **绘图设置** - 调整Y轴范围、温度刻度步长
6. **输出设置** - 指定输出目录和文件前缀

**工作流程：**

1. 选择热历史CSV文件 (包含Time/Myr和Avg_T/Celsius列)
2. 选择实测数据CSV文件 (包含Delta47和SD列)
3. 设置潜在热历史的时间窗口和要测试的峰值温度
4. 点击"运行分析"
5. 查看输出目录中的PDF/SVG图表

详细使用说明请参考: [GUI_USAGE.md](clump_history/GUI_USAGE.md)

---

### Input data format

Thermal history CSV (example: `Thermal_History_Hu.csv`) must include:

- `Time/Myr` - Time in millions of years (increasing into past, or decreasing time since present)
- `Avg_T/Celsius` - Average temperature in Celsius during each time segment

Actual test CSV (example: `acutal_test_Hu.csv`) must include:

- `Delta47` - Measured clumped isotope values (absolute, e.g. 0.440‰)
- `SD` - Standard deviation of measurements

(You can change these column names using CLI flags such as `--time-col`, `--avg-col`, `--d47-col`, and `--sd-col`.)

### Validation and Testing

Results have been validated against published datasets and laboratory heating experiments. See examples in the results section of this document.

### Contributing

We welcome contributions to improve the accuracy, efficiency, and applicability of these models. Please submit pull requests or open issues for bugs and feature suggestions.

### Citing this Software

When using this software in academic work, please cite both the original theoretical papers and reference this code repository:

Li, K. (2024). ClumpIsotope: Carbonate Clumped Isotope Geochemistry Modeling Framework [Computer software]. https://github.com/keran-li/ClumpIsotope

### Acknowledgments

This software framework builds upon the theoretical work of many researchers in the clumped isotope community. We acknowledge the pioneering work of the Eiler lab and others in developing this powerful geochemical tool.

## Packaging Applications

To create standalone executables for distribution:

### Automatic Build (Recommended)
Use the provided build script:
```bash
# On Windows
build_apps.bat
```

This will create both CLI and GUI executables in the `clump_history/dist/` directory.

### Manual Build with PyInstaller

#### GUI (recommended, no console window)
```bash
cd clump_history/
pyinstaller clump-history-gui.spec
```

Or using direct command:
```bash
pyinstaller run_gui.py ^
  --name clump-history-gui ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --collect-all matplotlib ^
  --collect-all isotopylog ^
  --collect-submodules scipy ^
  --collect-submodules numpy ^
  --collect-submodules pandas
```

#### CLI
```bash
cd clump_history/
pyinstaller clump-history.spec
```

Or using direct command:
```bash
pyinstaller run_cli.py ^
  --name clump-history ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --console ^
  --collect-all matplotlib ^
  --collect-all isotopylog ^
  --collect-submodules scipy ^
  --collect-submodules numpy ^
  --collect-submodules pandas
```

Executables will be created in the `dist/` directory:
- `clump-history-gui/clump-history-gui.exe` - Graphical user interface
- `clump-history/clump-history.exe` - Command-line interface

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
    <img width="800" alt="image" src="https://private-user-images.githubusercontent.com/66153455/532728773-15c6028b-a2be-4f7d-9b64-6d637246af08.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Njc3NzQ2NzcsIm5iZiI6MTc2Nzc3NDM3NywicGF0aCI6Ii82NjE1MzQ1NS81MzI3Mjg3NzMtMTVjNjAyOGItYTJiZS00ZjdkLTliNjQtNmQ2MzcyNDZhZjA4LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNjAxMDclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjYwMTA3VDA4MjYxN1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTM4OWM5YzBhMjQ1YjgwMzE4MTZiNTJkZDJiNjBlM2Q4OGQ3ZWFlZWIxN2E5ODVlYjdkYjhmYWY4NTlmODY5NTgmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.milYsqDn-2XhaD7y-L0YrmxgCs8mkWHYulNC6ZWN2nI">
    </div>

    (5) ours Sichuan samples:
    <div align="center">
    <img width="500" alt="image" src="https://private-user-images.githubusercontent.com/66153455/532729071-835c96f5-f44e-4a94-ad5e-89e6b6ab9de5.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Njc3NzQ3MjYsIm5iZiI6MTc2Nzc3NDQyNiwicGF0aCI6Ii82NjE1MzQ1NS81MzI3MjkwNzEtODM1Yzk2ZjUtZjQ0ZS00YTk0LWFkNWUtODllNmI2YWI5ZGU1LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNjAxMDclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjYwMTA3VDA4MjcwNlomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTg4MThhZjMwMzlmYmQ5NGY5ODQzNmJlYjA4ODI5NDIwNjlhODg3MTM2Y2E3MTAzNGIzODA3ZTg2YjdkMjM3MWMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.DNhTjq9N8UVRCSg1J488kDs5725qu6BnGpWPwYWWYGs">
    </div>

### To do

1. Add Monte-Carlo simulation for different Δ47 inputs (todo)
2. Export scenario results (D, Dstd, Deq) to CSV for each run (todo)