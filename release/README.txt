# Clump History GUI v0.2.0 发布说明

## 🎉 版本亮点

### 全新GUI界面
- **左右分栏布局**：左侧参数设置，右侧实时图表预览
- **交互式图表**：支持缩放、平移、保存等操作
- **中文界面**：友好的中文操作提示

### 核心功能
- **正演模拟**：基于Hemingway & Henkes (2021) HH21模型
- **U-Fit热历史调整**：在指定时间窗口内添加温度峰值
- **批量场景计算**：支持多个峰值温度同时计算
- **自动图表生成**：生成2×3子图对比不同场景

## 📦 文件说明

### 单文件版本（推荐普通用户）
```
ClumpHistoryGUI.exe    # 单个可执行文件，双击运行
```

### 文件夹版本（推荐开发者）
```
clump-history-gui/     # 包含所有依赖的完整文件夹
├── clump-history-gui.exe    # 主程序
├── 启动程序.bat              # 快速启动脚本
└── ...
```

## 🚀 快速开始

### 系统要求
- Windows 10/11 (64位)
- 内存：4GB+
- 磁盘空间：500MB+

### 安装步骤

#### 方法1：单文件版本（最简单）
1. 下载 `ClumpHistoryGUI.exe`
2. 双击运行
3. 开始使用！

#### 方法2：文件夹版本
1. 解压 `ClumpHistoryGUI.zip`
2. 运行 `启动程序.bat` 或直接运行 `clump-history-gui.exe`
3. 开始使用！

### 示例工作流程

1. **加载数据**
   - 热历史文件：`datasets/Thermal_History_Hu.csv`
   - 实测数据：`datasets/acutal_test_Hu.csv`

2. **设置潜在热历史**
   - 时间窗口：550 - 600 Myr
   - 峰值温度：`150 200 250 300 350`

3. **运行分析**
   - 点击"▶ 运行分析"
   - 右侧显示结果图表
   - 结果自动保存到 `output_figures/` 目录

## 📁 目录结构

```
工作目录/
├── ClumpHistoryGUI.exe          # 主程序
├── datasets/                     # 数据文件夹
│   ├── Thermal_History_Hu.csv   # 热历史数据
│   └── acutal_test_Hu.csv       # 实测Δ47数据
└── output_figures/              # 输出文件夹
    ├── forward_model.pdf        # PDF图表
    └── forward_model.svg        # SVG图表
```

## 🔬 科学背景

本软件基于以下理论：

1. **Stolper & Eiler (2015)**
   - 固态同位素交换反应动力学
   - American Journal of Science, 315(5), 363-401

2. **Hemingway & Henkes (2021)**
   - 团簇同位素键重排的无序动力学模型
   - Earth and Planetary Science Letters, 575, 117177

## ⚠️ 注意事项

1. **首次启动较慢**：单文件版本首次启动需要解压资源，请耐心等待
2. **数据文件路径**：建议使用相对路径或确保路径正确
3. **内存使用**：大量场景计算可能需要较多内存

## 📝 更新日志

### v0.2.0 (2026-02-01)
- ✨ 全新GUI界面，左右分栏设计
- ✨ 右侧实时图表预览
- ✨ 支持中文界面
- ✨ 单文件exe版本
- 🔧 改进的U-Fit算法
- 🔧 优化图表显示

### v0.1.0 (2026-01-11)
- 🎉 初始版本发布
- ✨ CLI命令行工具
- ✨ 基本GUI功能
- ✨ 正演模拟计算

## 📞 技术支持

如有问题，请通过以下方式联系：
- GitHub Issues: https://github.com/keran-li/ClumpIsotope/issues
- Email: keranli@smail.nju.edu.cn

## 📄 许可证

本项目采用 MIT 许可证
详见 LICENSE 文件

## 🙏 致谢

感谢团簇同位素地球化学社区的理论贡献，特别是Eiler实验室的开创性工作。

---

**发布日期**: 2026年2月1日  
**开发者**: Keran Li (Nanjing University)  
**版本**: v0.2.0
