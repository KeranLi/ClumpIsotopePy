import threading
from pathlib import Path
import os
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import isotopylog as ipl
from clump_history.model import compute_history
from clump_history.fit import constrained_u_fit
from clump_history.io import load_thermal_history, load_test_data
from clump_history.plot import plot_results

# 导入matplotlib的tkinter嵌入功能
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


class ClumpHistoryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Clump History GUI - 团簇同位素热历史正演模拟")
        self.root.geometry("1400x800")  # 增加宽度以容纳图表
        self.root.minsize(1400, 800)   # 设置最小窗口大小
        
        # 创建主框架，左右分栏
        self.main_paned = tk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 左侧：参数设置面板
        self.left_frame = tk.LabelFrame(self.main_paned, text="参数设置", padx=10, pady=10)
        self.main_paned.add(self.left_frame, width=550)  # 左侧固定宽度
        
        # 右侧：图表显示面板
        self.right_frame = tk.LabelFrame(self.main_paned, text="结果预览", padx=10, pady=10)
        self.main_paned.add(self.right_frame, width=800)  # 右侧较宽用于显示图表
        
        # 初始化UI组件
        self.setup_left_panel()
        self.setup_right_panel()
        
        # 设置默认值
        self.set_default_values()
    
    def setup_left_panel(self):
        """Setup the left panel with all input parameters"""
        # 创建Canvas和Scrollbar以支持滚动
        canvas = tk.Canvas(self.left_frame)
        scrollbar = tk.Scrollbar(self.left_frame, orient="vertical", command=canvas.yview)
        self.run_frame = tk.Frame(canvas)
        
        self.run_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.run_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # 绑定鼠标滚轮
        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        row = 0
        
        # === 文件选择区域 ===
        file_frame = tk.LabelFrame(self.run_frame, text="数据文件", padx=5, pady=5)
        file_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        # Thermal CSV
        tk.Label(file_frame, text="热历史 CSV:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.thermal_csv_var = tk.StringVar(value="datasets/Thermal_History_Hu.csv")
        self.thermal_csv_entry = tk.Entry(file_frame, textvariable=self.thermal_csv_var, width=40)
        self.thermal_csv_entry.grid(row=0, column=1, padx=5, pady=5)
        self.thermal_csv_button = tk.Button(file_frame, text="浏览", command=self.browse_thermal_csv)
        self.thermal_csv_button.grid(row=0, column=2, padx=5, pady=5)
        
        # Test CSV
        tk.Label(file_frame, text="实测数据 CSV:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.test_csv_var = tk.StringVar(value="datasets/acutal_test_Hu.csv")
        self.test_csv_entry = tk.Entry(file_frame, textvariable=self.test_csv_var, width=40)
        self.test_csv_entry.grid(row=1, column=1, padx=5, pady=5)
        self.test_csv_button = tk.Button(file_frame, text="浏览", command=self.browse_test_csv)
        self.test_csv_button.grid(row=1, column=2, padx=5, pady=5)
        
        row += 1
        
        # === 列名设置区域 ===
        col_frame = tk.LabelFrame(self.run_frame, text="CSV列名设置", padx=5, pady=5)
        col_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        tk.Label(col_frame, text="时间列:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.time_col_var = tk.StringVar(value="Time/Myr")
        tk.Entry(col_frame, textvariable=self.time_col_var, width=12).grid(row=0, column=1, padx=5, pady=2)
        
        tk.Label(col_frame, text="温度列:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)
        self.avg_col_var = tk.StringVar(value="Avg_T/Celsius")
        tk.Entry(col_frame, textvariable=self.avg_col_var, width=12).grid(row=0, column=3, padx=5, pady=2)
        
        tk.Label(col_frame, text="Δ47列:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.d47_col_var = tk.StringVar(value="Delta47")
        tk.Entry(col_frame, textvariable=self.d47_col_var, width=12).grid(row=1, column=1, padx=5, pady=2)
        
        tk.Label(col_frame, text="误差列:").grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)
        self.sd_col_var = tk.StringVar(value="SD")
        tk.Entry(col_frame, textvariable=self.sd_col_var, width=12).grid(row=1, column=3, padx=5, pady=2)
        
        row += 1
        
        # === 模型参数区域 ===
        model_frame = tk.LabelFrame(self.run_frame, text="模型参数", padx=5, pady=5)
        model_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        tk.Label(model_frame, text="矿物类型:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.mineral_var = tk.StringVar(value="calcite")
        tk.Entry(model_frame, textvariable=self.mineral_var, width=10).grid(row=0, column=1, padx=5, pady=2)
        
        tk.Label(model_frame, text="参考文献:").grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)
        self.ref_var = tk.StringVar(value="HH21")
        tk.Entry(model_frame, textvariable=self.ref_var, width=10).grid(row=0, column=3, padx=5, pady=2)
        
        tk.Label(model_frame, text="初始标准差:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.d0_std_var = tk.DoubleVar(value=0.02)
        tk.Entry(model_frame, textvariable=self.d0_std_var, width=10).grid(row=1, column=1, padx=5, pady=2)
        
        row += 1
        
        # === 潜在热历史设置 (U-Fit) ===
        peak_frame = tk.LabelFrame(self.run_frame, text="潜在热历史设置 (U-Fit)", padx=5, pady=5)
        peak_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        tk.Label(peak_frame, text="时间窗口 (Myr):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.peak_start_var = tk.DoubleVar(value=550)
        tk.Entry(peak_frame, textvariable=self.peak_start_var, width=10).grid(row=0, column=1, padx=5, pady=2)
        tk.Label(peak_frame, text="至").grid(row=0, column=2, padx=2)
        self.peak_end_var = tk.DoubleVar(value=600)
        tk.Entry(peak_frame, textvariable=self.peak_end_var, width=10).grid(row=0, column=3, padx=5, pady=2)
        
        tk.Label(peak_frame, text="峰值温度 (°C, 空格分隔):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.peak_temps_var = tk.StringVar(value="150 200 250 300 350")
        tk.Entry(peak_frame, textvariable=self.peak_temps_var, width=35).grid(row=1, column=1, columnspan=3, padx=5, pady=5)
        
        self.no_initial_var = tk.BooleanVar(value=False)
        tk.Checkbutton(peak_frame, text="不包含初始场景", variable=self.no_initial_var).grid(row=2, column=0, columnspan=2, padx=5, pady=2)
        
        row += 1
        
        # === 绘图设置区域 ===
        plot_frame = tk.LabelFrame(self.run_frame, text="绘图设置", padx=5, pady=5)
        plot_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        tk.Label(plot_frame, text="Y轴范围 (Δ47):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.ymin_var = tk.DoubleVar(value=0.15)
        tk.Entry(plot_frame, textvariable=self.ymin_var, width=8).grid(row=0, column=1, padx=5, pady=2)
        self.ymax_var = tk.DoubleVar(value=0.68)
        tk.Entry(plot_frame, textvariable=self.ymax_var, width=8).grid(row=0, column=2, padx=5, pady=2)
        
        tk.Label(plot_frame, text="温度刻度:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.tick_step_var = tk.DoubleVar(value=50)
        tk.Entry(plot_frame, textvariable=self.tick_step_var, width=8).grid(row=1, column=1, padx=5, pady=2)
        tk.Label(plot_frame, text="°C").grid(row=1, column=2, sticky=tk.W)
        
        row += 1
        
        # === 输出设置区域 ===
        out_frame = tk.LabelFrame(self.run_frame, text="输出设置", padx=5, pady=5)
        out_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
        
        tk.Label(out_frame, text="输出目录:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.outdir_var = tk.StringVar(value="output_figures")
        tk.Entry(out_frame, textvariable=self.outdir_var, width=30).grid(row=0, column=1, padx=5, pady=2)
        tk.Button(out_frame, text="浏览", command=self.browse_outdir).grid(row=0, column=2, padx=5, pady=2)
        
        tk.Label(out_frame, text="文件前缀:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.out_prefix_var = tk.StringVar(value="forward_model")
        tk.Entry(out_frame, textvariable=self.out_prefix_var, width=15).grid(row=1, column=1, sticky=tk.W, padx=5, pady=2)
        
        row += 1
        
        # === 按钮区域（使用独立Frame确保可见）===
        btn_frame = tk.LabelFrame(self.run_frame, text="操作", padx=10, pady=10)
        btn_frame.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=10)
        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        btn_frame.columnconfigure(2, weight=1)
        
        self.run_button = tk.Button(btn_frame, text="▶ 运行分析", command=self.run_analysis, 
                                    bg="#4CAF50", fg="white", font=("Arial", 11, "bold"),
                                    width=12, height=2, cursor="hand2")
        self.run_button.grid(row=0, column=0, padx=5, pady=5)
        
        self.open_outdir_button = tk.Button(btn_frame, text="📁 打开目录", command=self.open_outdir,
                                            width=12, height=2, font=("Arial", 10))
        self.open_outdir_button.grid(row=0, column=1, padx=5, pady=5)
        
        self.quit_button = tk.Button(btn_frame, text="✕ 退出", command=self.root.destroy,
                                     width=12, height=2, bg="#f0f0f0", font=("Arial", 10))
        self.quit_button.grid(row=0, column=2, padx=5, pady=5)
        
        # 状态栏（固定在底部）
        row += 1
        self.status_var = tk.StringVar(value="就绪 - 请设置参数后点击'运行分析'")
        self.status_bar = tk.Label(self.run_frame, textvariable=self.status_var, 
                                   bd=1, relief=tk.SUNKEN, anchor=tk.W,
                                   bg="#f8f8f8", fg="#333333", font=("Arial", 9))
        self.status_bar.grid(row=row, column=0, columnspan=3, sticky=tk.EW, pady=5)
    
    def setup_right_panel(self):
        """Setup the right panel for figure display"""
        # 创建matplotlib图形
        if MATPLOTLIB_AVAILABLE:
            # 创建图形和轴
            self.fig = Figure(figsize=(10, 8), dpi=100)
            self.ax = self.fig.add_subplot(111)
            
            # 嵌入到tkinter
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
            self.canvas_widget = self.canvas.get_tk_widget()
            self.canvas_widget.pack(fill=tk.BOTH, expand=True)
            
            # 初始显示提示文字
            self.ax.text(0.5, 0.5, '点击"运行分析"\n在此显示结果图表', 
                        ha='center', va='center', fontsize=16, color='gray')
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.canvas.draw()
            
            # 添加工具栏
            from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.right_frame)
            self.toolbar.update()
        else:
            # 如果matplotlib不可用，显示提示
            label = tk.Label(self.right_frame, text="matplotlib不可用\n请安装matplotlib以启用图表显示", 
                           font=("Arial", 14), fg="red")
            label.pack(expand=True)
    
    def browse_thermal_csv(self):
        """Open file dialog to select thermal history CSV"""
        filename = filedialog.askopenfilename(
            title="选择热历史CSV文件",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.thermal_csv_var.set(filename)
    
    def browse_test_csv(self):
        """Open file dialog to select test data CSV"""
        filename = filedialog.askopenfilename(
            title="选择实测数据CSV文件",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.test_csv_var.set(filename)
    
    def browse_outdir(self):
        """Open directory dialog to select output directory"""
        dirname = filedialog.askdirectory(title="选择输出目录")
        if dirname:
            self.outdir_var.set(dirname)
    
    def run_analysis(self):
        """Run the analysis with current parameters"""
        try:
            self.status_var.set("正在运行分析...")
            self.root.update()
            
            # Validate input paths
            thermal_path = self.thermal_csv_var.get().strip()
            test_path = self.test_csv_var.get().strip()
            
            # Convert relative paths to absolute paths
            if not os.path.isabs(thermal_path):
                thermal_path = os.path.join(os.getcwd(), thermal_path)
            
            if not os.path.isabs(test_path):
                test_path = os.path.join(os.getcwd(), test_path)
            
            # Check if files exist
            if not os.path.exists(thermal_path):
                messagebox.showerror("错误", f"热历史文件未找到: {thermal_path}")
                self.status_var.set("错误: 文件未找到")
                return
            
            if not os.path.exists(test_path):
                messagebox.showerror("错误", f"实测数据文件未找到: {test_path}")
                self.status_var.set("错误: 文件未找到")
                return
            
            # Load data
            self.status_var.set("正在加载数据...")
            self.root.update()
            
            thermal_df = pd.read_csv(thermal_path)
            time_col = self.time_col_var.get()
            temp_col = self.avg_col_var.get()
            
            if time_col not in thermal_df.columns:
                raise ValueError(f"时间列 '{time_col}' 不存在于热历史文件中")
            if temp_col not in thermal_df.columns:
                raise ValueError(f"温度列 '{temp_col}' 不存在于热历史文件中")
            
            time_myr = thermal_df[time_col].values
            avg_temp_c = thermal_df[temp_col].values
            avg_temp_k = avg_temp_c + 273.15
            
            # Load test data
            test_df = pd.read_csv(test_path)
            d47_col = self.d47_col_var.get()
            sd_col = self.sd_col_var.get()
            
            if d47_col not in test_df.columns:
                raise ValueError(f"Δ47列 '{d47_col}' 不存在于实测数据文件中")
            if sd_col not in test_df.columns:
                raise ValueError(f"误差列 '{sd_col}' 不存在于实测数据文件中")
            
            delta47 = test_df[d47_col].values
            delta47_err = test_df[sd_col].values
            
            # Parse peak temperatures
            peak_temps_str = self.peak_temps_var.get().strip()
            if not peak_temps_str:
                raise ValueError("请至少输入一个峰值温度")
            
            peak_temps = [float(t) for t in peak_temps_str.split()]
            peak_start = self.peak_start_var.get()
            peak_end = self.peak_end_var.get()
            
            # Get model parameters
            mineral = self.mineral_var.get()
            ref = self.ref_var.get()
            d0_std = self.d0_std_var.get()
            
            # Create output directory
            outdir = self.outdir_var.get().strip()
            if not os.path.isabs(outdir):
                outdir = os.path.join(os.getcwd(), outdir)
            os.makedirs(outdir, exist_ok=True)
            
            # Get EDistribution
            self.status_var.set("正在加载活化能分布...")
            self.root.update()
            
            ed = ipl.EDistribution.from_literature(
                mineral=mineral,
                reference=ref
            )
            
            # Build scenarios
            self.status_var.set("正在构建热历史场景...")
            self.root.update()
            
            scenarios = []
            
            # Initial scenario (no peak adjustment)
            if not self.no_initial_var.get():
                D, Dstd, Deq = compute_history(time_myr, avg_temp_k, ed, d0_std)
                scenarios.append(('initial', D, Dstd, Deq))
            
            # Peak temperature scenarios
            for temp_c in peak_temps:
                temp_k = temp_c + 273.15
                adjusted_temp = constrained_u_fit(
                    time_myr, avg_temp_k, peak_start, peak_end, temp_k, plot=False
                )
                D, Dstd, Deq = compute_history(time_myr, adjusted_temp, ed, d0_std)
                scenarios.append((f'{int(temp_c)}', D, Dstd, Deq))
            
            # Prepare results
            results = {
                'time_myr': time_myr,
                'scenarios': scenarios,
                'delta47': delta47,
                'delta47_err': delta47_err
            }
            
            # Save results and plot
            self.status_var.set("正在生成图表...")
            self.root.update()
            
            out_prefix = self.out_prefix_var.get().strip()
            plot_results(
                results, 
                outdir, 
                out_prefix,
                show_plots=False,  # 不在外部显示
                ymin=self.ymin_var.get(), 
                ymax=self.ymax_var.get(),
                tick_step=self.tick_step_var.get()
            )
            
            # 更新右侧图表显示
            self.status_var.set("正在更新图表显示...")
            self.root.update()
            self.update_figure_display(results)
            
            self.status_var.set(f"分析完成! 结果保存至: {outdir}")
            messagebox.showinfo("成功", f"分析完成!\n结果保存至:\n{outdir}")
            
        except Exception as e:
            import traceback
            error_msg = f"{str(e)}\n\n详细信息:\n{traceback.format_exc()}"
            messagebox.showerror("错误", f"分析过程中出现错误:\n{error_msg}")
            self.status_var.set("错误")
    
    def update_figure_display(self, results):
        """Update the figure display in the right panel"""
        if not MATPLOTLIB_AVAILABLE:
            return
        
        # 清除当前图形
        self.fig.clear()
        
        # 获取参数
        ymin = self.ymin_var.get()
        ymax = self.ymax_var.get()
        tick_step = self.tick_step_var.get()
        
        # 准备数据
        time_myr = results['time_myr']
        scenarios = results['scenarios']
        delta47 = results['delta47']
        delta47_err = results['delta47_err']
        
        # 计算子图布局
        n_scenarios = len(scenarios)
        n_cols = 3
        n_rows = (n_scenarios + n_cols - 1) // n_cols
        
        # 创建子图
        from clump_history.plot import build_secondary_ticks
        d47_ticks, T_labels = build_secondary_ticks(ymin, ymax, tick_step)
        
        for i, (label, D, Dstd, Deq) in enumerate(scenarios):
            ax = self.fig.add_subplot(n_rows, n_cols, i + 1)
            
            # 绘制曲线
            ax.plot(time_myr, D, label=f'Forward-modeled', linewidth=1.5, color='C0')
            ax.fill_between(time_myr, D - Dstd, D + Dstd, alpha=0.35, color='C0')
            ax.plot(time_myr, Deq, label=f'Equilibrium', linewidth=1, color='C1', linestyle='--')
            
            # 绘制实测数据
            ax.errorbar(
                np.zeros_like(delta47),
                delta47, yerr=delta47_err,
                fmt='o',
                label='Actual',
                color='black',
                capsize=4,
                alpha=0.6,
                markerfacecolor='none'
            )
            
            ax.set_ylim(ymin, ymax)
            ax.set_xlabel('Age (Myr)', fontsize=9)
            ax.set_ylabel(r'$\Delta$47 (‰)', fontsize=9)
            ax.set_title(f'{label}°C', fontsize=10, fontweight='bold')
            ax.legend(loc='best', fontsize=7)
            ax.grid(True, alpha=0.3)
            
            # 添加右侧温度轴
            if i % n_cols == n_cols - 1:
                secax = ax.secondary_yaxis('right')
                secax.set_yticks(d47_ticks)
                secax.set_yticklabels(T_labels, fontsize=7)
                if i == n_cols - 1:
                    secax.set_ylabel('Temperature (°C)', fontsize=9)
        
        # 调整布局
        self.fig.tight_layout()
        
        # 刷新画布
        self.canvas.draw()
    
    def open_outdir(self):
        """Open the output directory in file explorer"""
        outdir = self.outdir_var.get().strip()
        if not os.path.isabs(outdir):
            outdir = os.path.join(os.getcwd(), outdir)
        
        if os.path.exists(outdir):
            import subprocess
            subprocess.Popen(['explorer', outdir])
        else:
            messagebox.showerror("错误", f"输出目录不存在: {outdir}")
    
    def set_default_values(self):
        """Set default values for all parameters"""
        # Default values are already set in the initialization
        pass


def main():
    """Main entry point for the GUI application"""
    root = tk.Tk()
    app = ClumpHistoryGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
