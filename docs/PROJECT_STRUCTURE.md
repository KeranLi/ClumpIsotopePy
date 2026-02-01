# Clump History - 项目结构说明

## 📁 根目录结构

```
ClumpIsotope/
│
├── README.md                      # 项目主README
├── QUICKSTART_GUI.md             # GUI快速开始指南
├── USER_GUIDE.md                 # 用户完全指南
├── RELEASE_NOTES.md              # 发布说明
├── CHECKLIST_BEFORE_RELEASE.md   # 发布前检查清单
├── TEST_GUIDE.md                 # 测试指南
├── THEORY.md                     # 理论基础文档
├── LICENSE                       # MIT许可证
│
├── build_onefile.bat             # 单文件版本构建脚本
├── build_executable.bat          # 完整版本构建脚本
├── create_env.yaml               # Conda环境配置文件
├── test_gui.py                   # GUI测试启动脚本
├── test_core_functions.py        # 核心功能测试
│
├── clump_history/                # 主程序包
│   ├── pyproject.toml            # Python包配置
│   ├── run_gui.py                # GUI入口点
│   ├── run_cli.py                # CLI入口点
│   │
│   ├── clump-history-gui.spec            # PyInstaller GUI配置（文件夹版）
│   ├── clump-history-gui-onefile.spec    # PyInstaller GUI配置（单文件版）
│   └── clump-history.spec                # PyInstaller CLI配置
│   │
│   └── src/
│       └── clump_history/        # Python包源代码
│           ├── __init__.py       # 包初始化
│           ├── __main__.py       # 模块入口点
│           ├── cli.py            # 命令行界面
│           ├── gui.py            # 图形用户界面
│           ├── model.py          # 正演模型计算
│           ├── fit.py            # U-Fit热历史调整
│           ├── io.py             # 数据输入输出
│           └── plot.py           # 绘图功能
│
└── datasets/                     # 示例数据集
    ├── Thermal_History_Hu.csv    # 胡教授研究区热历史
    ├── Thermal_History.csv       # 标准热历史
    ├── acutal_test_Hu.csv        # 实测数据（Hu）
    └── acutal_test.csv           # 实测数据（标准）
```

## 📦 发布包结构

### 单文件版本（推荐用户）
```
release/
├── ClumpHistoryGUI.exe          # 主程序（单文件）
├── 启动程序.bat                  # 中文启动脚本
├── README.txt                    # 发布说明
├── 使用说明.txt                  # 中文使用说明
│
└── datasets/                     # 示例数据
    ├── Thermal_History_Hu.csv
    └── acutal_test_Hu.csv
```

### 文件夹版本（完整功能）
```
ClumpHistoryGUI/
├── clump-history-gui.exe         # 主程序
├── 启动程序.bat                  # 启动脚本
│
├── datasets/                     # 示例数据
│   ├── Thermal_History_Hu.csv
│   └── acutal_test_Hu.csv
│
└── docs/                         # 文档
    ├── README.md
    ├── QUICKSTART_GUI.md
    └── THEORY.md
```

## 🧩 核心模块说明

### gui.py (23KB)
- **功能**: 图形用户界面
- **特点**: 
  - 左右分栏布局
  - 左侧参数设置
  - 右侧实时图表预览
  - 中文界面支持

### model.py (6KB)
- **功能**: 正演模拟计算核心
- **依赖**: isotopylog, numpy
- **接口**: `compute_history(time_myr, T_k, ed, d0_std)`

### fit.py (7KB)
- **功能**: U-Fit热历史调整
- **算法**: 四次多项式拟合U型曲线
- **接口**: `constrained_u_fit(time, values, start, end, max_temp)`

### plot.py (4KB)
- **功能**: 图表生成
- **特点**: 2×3子图布局，支持温度双Y轴
- **输出**: PDF/SVG矢量图

### io.py (3KB)
- **功能**: 数据文件读写
- **支持**: CSV格式，可配置列名

## 🔧 配置文件说明

### pyproject.toml
- 包名称: `clump-history`
- 版本: `0.2.0`
- Python版本: `>=3.7,<3.13`
- NumPy版本: `<2.0` (isotopylog兼容性要求)

### *.spec (PyInstaller配置)
- **clump-history-gui.spec**: 文件夹版本（启动快，体积大）
- **clump-history-gui-onefile.spec**: 单文件版本（启动慢，体积小）
- **clump-history.spec**: CLI版本

## 📝 关键设计决策

### 1. 左右分栏GUI
- 左侧固定宽度550px，容纳所有参数
- 右侧动态宽度800px+，显示图表
- 使用PanedWindow支持调整比例

### 2. 实时图表嵌入
- 使用matplotlib的FigureCanvasTkAgg
- 包含NavigationToolbar2Tk工具栏
- 自动布局子图（根据场景数量）

### 3. 单文件vs文件夹版本
- **单文件**: 用户友好，首次启动慢
- **文件夹**: 启动快，适合频繁使用

### 4. 依赖管理
- 核心依赖: numpy<2.0, pandas, matplotlib, scipy, isotopylog
- 开发依赖: pytest, black, flake8, pyinstaller
- 文档依赖: sphinx

## 🚀 开发工作流

### 开发阶段
```bash
# 1. 创建开发环境
conda env create -f create_env.yaml
conda activate clump-isotope

# 2. 测试修改
cd clump_history
python run_gui.py

# 3. 功能验证
python ../test_core_functions.py
```

### 打包阶段
```bash
# 构建单文件版本（推荐）
build_onefile.bat

# 或构建文件夹版本
build_executable.bat
```

### 发布阶段
```bash
# 1. 版本号检查
# 更新 pyproject.toml 中的 version

# 2. 文档更新
# 更新 RELEASE_NOTES.md

# 3. 构建
build_onefile.bat

# 4. 验证
# 测试生成的exe文件

# 5. 发布
# 上传到GitHub Release
```

## 📊 性能指标

### 启动时间
- 文件夹版本: < 5秒
- 单文件版本: 1-2分钟（首次）

### 计算时间
- 单个场景: ~5秒
- 6个场景（默认）: ~30秒

### 内存占用
- 空闲: ~150MB
- 计算中: ~500MB

### 输出文件大小
- PDF: ~500KB
- SVG: ~200KB

## 🔐 安全考虑

1. **杀毒软件误报**
   - PyInstaller打包的exe可能被误报
   - 建议添加到白名单

2. **数据隐私**
   - 所有计算在本地进行
   - 不连接网络
   - 不收集用户数据

3. **文件权限**
   - 需要读写工作目录
   - 建议以普通用户运行

## 🎯 用户群体

### 主要用户
- 地球化学研究人员
- 地质学研究生
- 石油地质分析师

### 使用场景
- 碳酸盐团簇同位素研究
- 热历史反演分析
- 古温度重建
- 成岩作用研究

---

**维护者**: Keran Li  
**版本**: v0.2.0  
**最后更新**: 2026-02-01
