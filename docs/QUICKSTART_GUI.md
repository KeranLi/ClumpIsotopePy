# 团簇同位素热历史正演模拟 - GUI快速开始指南

## 概述

本工具提供了一个图形用户界面(GUI)，用于模拟碳酸盐团簇同位素(Δ47)在不同热历史条件下的固态重排过程。

## 启动GUI

### 方式1: 直接运行测试脚本
```bash
python test_gui.py
```

### 方式2: 从clump_history目录运行
```bash
cd clump_history
python run_gui.py
```

### 方式3: 安装后运行 (推荐)
```bash
pip install -e ./clump_history
clump-history-gui
```

## 快速入门步骤

### 第1步: 准备数据文件

#### 热历史CSV文件格式示例 (`Thermal_History.csv`):
```csv
Time/Myr,Avg_T/Celsius
0,50
10,55
20,60
...
500,120
```

#### 实测数据CSV文件格式示例 (`test_data.csv`):
```csv
Delta47,SD
0.520,0.015
0.518,0.012
...
```

### 第2步: 配置GUI参数

1. **选择文件**
   - 热历史CSV: 选择包含时间和温度数据的文件
   - 实测数据CSV: 选择包含Δ47和误差的文件

2. **设置潜在热历史**
   - 时间窗口: 输入要调整的起始和结束时间 (Myr)
     - 例如: 550 至 600 Myr
   - 峰值温度: 输入要测试的温度值，空格分隔
     - 例如: `150 200 250 300 350`

3. **选择矿物类型**
   - 方解石: `calcite`
   - 白云石: `dolomite`

### 第3步: 运行分析

点击"运行分析"按钮，等待计算完成。

### 第4步: 查看结果

- 图表将保存到指定的输出目录
- 生成PDF和SVG格式的矢量图
- 图表包含多个子图，每个峰值温度一个

## 输出图表说明

每个子图显示以下内容：
- **蓝色实线**: 正演模拟的Δ47演化曲线
- **浅蓝色阴影**: 模拟结果的不确定性范围
- **橙色实线**: 平衡值曲线 (无动力学障碍的理论值)
- **黑色圆点**: 实测Δ47数据 (在time=0处)
- **右侧Y轴**: 对应的温度刻度

## 示例工作流程

### 示例1: 测试不同峰值温度的影响

1. 加载 `datasets/Thermal_History_Hu.csv`
2. 加载 `datasets/acutal_test_Hu.csv`
3. 设置时间窗口: 550 - 600 Myr
4. 设置峰值温度: `100 150 200 250 300`
5. 运行分析
6. 查看哪个峰值温度场景最符合实测数据

### 示例2: 白云石样品分析

1. 加载热历史和实测数据
2. 修改矿物类型为 `dolomite`
3. 设置合适的峰值温度
4. 运行分析

## 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| 时间列 | 热历史CSV中的时间列名 | Time/Myr |
| 温度列 | 热历史CSV中的温度列名 | Avg_T/Celsius |
| Δ47列 | 实测数据CSV中的Δ47列名 | Delta47 |
| 误差列 | 实测数据CSV中的误差列名 | SD |
| 矿物类型 | calcite 或 dolomite | calcite |
| 初始标准差 | 初始Δ47的不确定性 | 0.02 |
| 时间窗口 | U-Fit调整的时间范围 | 550-600 Myr |
| 峰值温度 | 要测试的最大温度值 | 150 200 250 300 350 |

## 故障排除

| 问题 | 解决方案 |
|------|---------|
| 文件未找到 | 使用绝对路径或确保文件存在于工作目录 |
| 列不存在 | 检查CSV列名设置是否与实际文件匹配 |
| 图表为空 | 检查时间窗口是否在数据范围内 |
| 计算结果异常 | 确认温度单位为摄氏度，矿物类型正确 |

## 理论基础

本工具基于以下理论模型：

1. **Stolper & Eiler (2015)**
   - 固态同位素交换反应动力学
   - American Journal of Science, 315(5), 363-401

2. **Hemingway & Henkes (2021)**
   - 团簇同位素键重排的无序动力学模型
   - Earth and Planetary Science Letters, 575, 117177

3. **U-Fit方法**
   - 在保持地质约束的同时引入温度峰值
   - 使用四次多项式拟合U型温度曲线

## 数据文件示例位置

项目中包含以下示例数据：
- `datasets/Thermal_History_Hu.csv` - 胡教授研究区热历史
- `datasets/Thermal_History.csv` - 标准热历史
- `datasets/acutal_test_Hu.csv` - 实测Δ47数据
- `datasets/acutal_test.csv` - 另一组实测数据

## 获取帮助

- 详细使用说明: [GUI_USAGE.md](clump_history/GUI_USAGE.md)
- 项目README: [README.md](README.md)
- 理论基础: [THEORY.md](THEORY.md)

## 引用

使用本工具请引用：

Li, K. (2024). ClumpIsotope: Carbonate Clumped Isotope Geochemistry Modeling Framework [Computer software]. https://github.com/keran-li/ClumpIsotope

以及原始理论文献：
- Stolper, D. A., & Eiler, J. M. (2015). American Journal of Science, 315(5), 363-401.
- Hemingway, J. D., & Henkes, G. A. (2021). Earth and Planetary Science Letters, 575, 117177.
