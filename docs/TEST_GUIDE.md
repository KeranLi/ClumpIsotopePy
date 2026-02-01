# Clump History GUI - 测试指南

## 环境要求

GUI工具需要以下依赖：
- Python 3.7-3.10 (推荐 3.10)
- NumPy 1.23 (注意：NumPy 2.x 与 isotopylog 0.0.8 不兼容)
- pandas
- matplotlib
- scipy
- isotopylog 0.0.8
- tkinter (通常随Python安装)

## 方案1: 使用 Conda 环境（推荐）

### 步骤1: 创建环境

```bash
# 在项目根目录下
conda env create -f create_env.yaml
```

### 步骤2: 激活环境

```bash
conda activate clump-isotope
```

### 步骤3: 测试核心功能

```bash
python test_core_functions.py
```

如果看到 "All tests passed!"，说明核心功能正常。

### 步骤4: 启动GUI

```bash
python test_gui.py
```

或者：

```bash
cd clump_history
python run_gui.py
```

## 方案2: 使用现有 Python 环境

如果你的系统 Python 已经安装了正确版本的依赖：

### 检查依赖版本

```bash
python -c "import numpy; print('numpy:', numpy.__version__)"
```

如果显示 `1.x.x`，可以继续测试。

### 测试核心功能

```bash
python test_core_functions.py
```

### 启动GUI

```bash
python test_gui.py
```

## 方案3: 使用 Jupyter Notebook 环境

如果你的 notebook 环境可以正常运行 `forward_hu.ipynb`，则可以在同一环境中测试GUI：

```python
# 在可以运行 notebook 的环境中
import sys
sys.path.insert(0, 'clump_history/src')
from clump_history.gui import main
main()
```

## 常见问题

### 问题1: "ValueError: setting an array element with a sequence"

**原因**: isotopylog 0.0.8 与 NumPy 2.x 不兼容

**解决方案**: 使用 Conda 环境安装 NumPy 1.23

```bash
conda env create -f create_env.yaml
conda activate clump-isotope
python test_gui.py
```

### 问题2: "ModuleNotFoundError: No module named 'isotopylog'"

**解决方案**: 安装 isotopylog

```bash
pip install isotopylog
```

或在 conda 环境中：

```bash
conda activate clump-isotope
pip install isotopylog
```

### 问题3: "ModuleNotFoundError: No module named 'tkinter'"

**Windows 解决方案**: 重新安装 Python 并勾选 "tcl/tk and IDLE"

**Conda 解决方案**:

```bash
conda install tk
```

### 问题4: GUI窗口不显示

检查是否有报错信息。常见原因：
1. 缺少 tkinter
2. 依赖版本不兼容
3. 需要在图形界面环境下运行（不能在纯SSH终端）

## 测试流程

### 功能测试清单

1. **文件选择**
   - [ ] 能正常选择热历史CSV文件
   - [ ] 能正常选择实测数据CSV文件
   - [ ] 能正常选择输出目录

2. **参数设置**
   - [ ] 可以修改CSV列名
   - [ ] 可以设置矿物类型
   - [ ] 可以设置时间窗口和峰值温度

3. **运行分析**
   - [ ] 点击"运行分析"后状态栏显示进度
   - [ ] 计算完成后弹出成功提示
   - [ ] 输出目录生成PDF和SVG文件

4. **结果查看**
   - [ ] 图表包含多个子图（每个峰值温度一个）
   - [ ] 每个子图显示正演曲线、误差区域、平衡值和实测数据
   - [ ] 右侧Y轴显示温度刻度

### 示例数据测试

使用项目自带的示例数据进行测试：

1. **热历史**: `datasets/Thermal_History_Hu.csv`
2. **实测数据**: `datasets/acutal_test_Hu.csv`
3. **参数设置**:
   - 时间窗口: 550 - 600 Myr
   - 峰值温度: 150 200 250 300 350
   - 矿物类型: calcite

预期结果：大约30-60秒后生成图表，显示6个子图（初始+5个峰值温度）。

## 验证安装

运行以下命令验证所有组件：

```bash
# 1. 检查Python版本
python --version  # 应为 3.7-3.10

# 2. 检查NumPy版本
python -c "import numpy; print('numpy:', numpy.__version__)"  # 应为 1.x

# 3. 检查isotopylog
python -c "import isotopylog; print('isotopylog: OK')"

# 4. 检查tkinter
python -c "import tkinter; print('tkinter: OK')"

# 5. 测试模块导入
python -c "import sys; sys.path.insert(0, 'clump_history/src'); from clump_history.gui import ClumpHistoryGUI; print('GUI module: OK')"
```

## 获取帮助

如果测试过程中遇到问题：

1. 检查错误消息中的具体信息
2. 确认使用的Python环境正确 (`which python` 或 `where python`)
3. 查看项目文档: `GUI_USAGE.md`, `README.md`
4. 在GitHub上提交issue（如果是代码问题）

## 下一步

测试通过后，你可以：

1. 使用自己的数据进行正演模拟
2. 尝试不同的峰值温度组合
3. 对比不同矿物类型（calcite vs dolomite）
4. 导出结果用于论文或报告
