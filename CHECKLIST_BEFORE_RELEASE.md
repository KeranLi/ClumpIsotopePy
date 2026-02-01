# 发布前检查清单

## ✅ 代码检查

- [ ] `clump_history/src/clump_history/gui.py` - 最新版本，包含右侧图表
- [ ] `clump_history/src/clump_history/model.py` - 正确实现
- [ ] `clump_history/src/clump_history/fit.py` - U-Fit功能正常
- [ ] `clump_history/src/clump_history/plot.py` - 绘图功能正常
- [ ] `clump_history/src/clump_history/io.py` - 数据导入正常

## ✅ 配置文件

- [ ] `clump_history/pyproject.toml` - 版本号正确 (0.2.0)
- [ ] `clump_history/clump-history-gui.spec` - PyInstaller配置正确
- [ ] `clump_history/clump-history-gui-onefile.spec` - 单文件配置正确

## ✅ 数据文件

- [ ] `datasets/Thermal_History_Hu.csv` - 示例热历史数据
- [ ] `datasets/acutal_test_Hu.csv` - 示例实测数据

## ✅ 文档

- [ ] `README.md` - 项目介绍
- [ ] `QUICKSTART_GUI.md` - GUI快速开始
- [ ] `GUI_USAGE.md` - 详细使用说明
- [ ] `RELEASE_NOTES.md` - 发布说明
- [ ] `THEORY.md` - 理论基础

## ✅ 构建脚本

- [ ] `build_executable.bat` - 完整版本构建
- [ ] `build_onefile.bat` - 单文件版本构建
- [ ] `create_env.yaml` - Conda环境配置

## 🧪 测试检查

### 功能测试
- [ ] GUI启动正常
- [ ] 文件选择功能正常
- [ ] 运行分析功能正常
- [ ] 右侧图表显示正常
- [ ] 工具栏（缩放/平移）正常
- [ ] 结果保存正常

### 数据测试
- [ ] 使用示例数据能正常运行
- [ ] 多个峰值温度场景计算正确
- [ ] 图表生成正确

### 打包测试
- [ ] 单文件exe能正常运行
- [ ] 文件夹版本能正常运行
- [ ] 首次启动时间可接受（< 2分钟）

## 📦 发布文件

### 单文件版本
```
release/
├── ClumpHistoryGUI.exe      # 主程序
├── 启动程序.bat              # 启动脚本
├── datasets/                 # 示例数据
│   ├── Thermal_History_Hu.csv
│   └── acutal_test_Hu.csv
├── README.txt               # 发布说明
└── 使用说明.txt              # 使用说明
```

### 文件夹版本
```
release_folder/
└── ClumpHistoryGUI/
    ├── clump-history-gui.exe
    ├── 启动程序.bat
    ├── datasets/            # 示例数据
    └── docs/               # 文档
```

## 🚀 发布步骤

1. **测试通过**
   ```bash
   # 运行核心功能测试
   python test_core_functions.py
   
   # 启动GUI测试
   python test_gui.py
   ```

2. **构建单文件版本**（推荐）
   ```bash
   build_onefile.bat
   ```

3. **验证构建**
   - 检查 `release/ClumpHistoryGUI.exe` 存在
   - 双击运行，测试所有功能

4. **创建发布包**
   - 压缩 `release/` 文件夹
   - 命名为 `ClumpHistoryGUI-v0.2.0.zip`

5. **发布**
   - 上传到GitHub Release
   - 更新README下载链接
   - 通知用户

## 📋 发布后检查

- [ ] GitHub Release创建成功
- [ ] 下载链接可访问
- [ ] 文件大小合理（< 500MB）
- [ ] 用户能正常下载使用

## 🔧 常见问题

### 构建失败
- 确保在正确的conda环境中
- 确保PyInstaller已安装
- 检查spec文件路径正确

### 运行时错误
- 确保所有依赖包含在spec文件中
- 检查matplotlib后端配置
- 检查isotopylog数据文件包含

### 文件过大
- 使用UPX压缩（已在spec中启用）
- 排除不必要的测试模块
- 考虑使用文件夹版本

---

**发布版本**: v0.2.0  
**发布日期**: 2026-02-01  
**负责人**: Keran Li
