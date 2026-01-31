import threading
from pathlib import Path
import os
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
from clump_history.model import compute_history
from clump_history.fit import constrained_u_fit
from clump_history.io import load_thermal_history, load_test_data
from clump_history.plot import plot_results  # 确保此导入正确

class ClumpHistoryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Clump History GUI")
        
        # Create tabs
        self.tabs = tk.Frame(root)
        self.tabs.pack(fill=tk.BOTH, expand=True)
        
        # Run tab
        self.run_frame = tk.Frame(self.tabs)
        self.run_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Ufit tab
        self.ufit_frame = tk.Frame(self.tabs)
        self.ufit_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Initialize UI components
        self.setup_run_tab()
        self.setup_uft_tab()
        
        # Set default values
        self.set_default_values()
    
    def setup_run_tab(self):
        # Thermal CSV
        tk.Label(self.run_frame, text="Thermal CSV:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.thermal_csv_var = tk.StringVar(value="datasets/Thermal_History_Hu.csv")
        self.thermal_csv_entry = tk.Entry(self.run_frame, textvariable=self.thermal_csv_var, width=50)
        self.thermal_csv_entry.grid(row=0, column=1, padx=5, pady=5)
        self.thermal_csv_button = tk.Button(self.run_frame, text="Browse", command=self.browse_thermal_csv)
        self.thermal_csv_button.grid(row=0, column=2, padx=5, pady=5)
        
        # Test CSV
        tk.Label(self.run_frame, text="Test CSV:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.test_csv_var = tk.StringVar(value="datasets/acutal_test_Hu.csv")
        self.test_csv_entry = tk.Entry(self.run_frame, textvariable=self.test_csv_var, width=50)
        self.test_csv_entry.grid(row=1, column=1, padx=5, pady=5)
        self.test_csv_button = tk.Button(self.run_frame, text="Browse", command=self.browse_test_csv)
        self.test_csv_button.grid(row=1, column=2, padx=5, pady=5)
        
        # Outdir
        tk.Label(self.run_frame, text="Outdir:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
        self.outdir_var = tk.StringVar(value="results")
        self.outdir_entry = tk.Entry(self.run_frame, textvariable=self.outdir_var, width=50)
        self.outdir_entry.grid(row=2, column=1, padx=5, pady=5)
        
        # Out prefix
        tk.Label(self.run_frame, text="Out prefix:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)
        self.out_prefix_var = tk.StringVar(value="fig_smoke")
        self.out_prefix_entry = tk.Entry(self.run_frame, textvariable=self.out_prefix_var, width=50)
        self.out_prefix_entry.grid(row=3, column=1, padx=5, pady=5)
        
        # Columns
        tk.Label(self.run_frame, text="Cols (time, avgT):").grid(row=4, column=0, sticky=tk.W, padx=5, pady=5)
        self.time_col_var = tk.StringVar(value="Time/Myr")
        self.time_col_entry = tk.Entry(self.run_frame, textvariable=self.time_col_var, width=20)
        self.time_col_entry.grid(row=4, column=1, padx=5, pady=5)
        self.avg_col_var = tk.StringVar(value="Avg_T/Celsius")
        self.avg_col_entry = tk.Entry(self.run_frame, textvariable=self.avg_col_var, width=20)
        self.avg_col_entry.grid(row=4, column=2, padx=5, pady=5)
        
        tk.Label(self.run_frame, text="Cols (Δ47, SD):").grid(row=5, column=0, sticky=tk.W, padx=5, pady=5)
        self.d47_col_var = tk.StringVar(value="Delta47")
        self.d47_col_entry = tk.Entry(self.run_frame, textvariable=self.d47_col_var, width=20)
        self.d47_col_entry.grid(row=5, column=1, padx=5, pady=5)
        self.sd_col_var = tk.StringVar(value="SD")
        self.sd_col_entry = tk.Entry(self.run_frame, textvariable=self.sd_col_var, width=20)
        self.sd_col_entry.grid(row=5, column=2, padx=5, pady=5)
        
        # Mineral / Ref
        tk.Label(self.run_frame, text="Mineral / Ref:").grid(row=6, column=0, sticky=tk.W, padx=5, pady=5)
        self.mineral_var = tk.StringVar(value="calcite")
        self.mineral_entry = tk.Entry(self.run_frame, textvariable=self.mineral_var, width=20)
        self.mineral_entry.grid(row=6, column=1, padx=5, pady=5)
        self.ref_var = tk.StringVar(value="HH21")
        self.ref_entry = tk.Entry(self.run_frame, textvariable=self.ref_var, width=20)
        self.ref_entry.grid(row=6, column=2, padx=5, pady=5)
        
        # d0_std
        tk.Label(self.run_frame, text="d0_std:").grid(row=7, column=0, sticky=tk.W, padx=5, pady=5)
        self.d0_std_var = tk.DoubleVar(value=0.02)
        self.d0_std_entry = tk.Entry(self.run_frame, textvariable=self.d0_std_var, width=20)
        self.d0_std_entry.grid(row=7, column=1, padx=5, pady=5)
        
        # Peak window
        tk.Label(self.run_frame, text="Peak window (start,end Myr):").grid(row=8, column=0, sticky=tk.W, padx=5, pady=5)
        self.peak_start_var = tk.DoubleVar(value=550)
        self.peak_start_entry = tk.Entry(self.run_frame, textvariable=self.peak_start_var, width=10)
        self.peak_start_entry.grid(row=8, column=1, padx=5, pady=5)
        self.peak_end_var = tk.DoubleVar(value=600)
        self.peak_end_entry = tk.Entry(self.run_frame, textvariable=self.peak_end_var, width=10)
        self.peak_end_entry.grid(row=8, column=2, padx=5, pady=5)
        
        # Peak temps
        tk.Label(self.run_frame, text="Peak temps (°C):").grid(row=9, column=0, sticky=tk.W, padx=5, pady=5)
        self.peak_temps_var = tk.StringVar(value="150 200 250 300 350")
        self.peak_temps_entry = tk.Entry(self.run_frame, textvariable=self.peak_temps_var, width=50)
        self.peak_temps_entry.grid(row=9, column=1, padx=5, pady=5)
        
        # Checkboxes
        self.no_initial_var = tk.BooleanVar(value=False)
        self.no_initial_check = tk.Checkbutton(self.run_frame, text="No initial scenario", variable=self.no_initial_var)
        self.no_initial_check.grid(row=10, column=0, padx=5, pady=5)
        
        self.show_plots_var = tk.BooleanVar(value=False)
        self.show_plots_check = tk.Checkbutton(self.run_frame, text="Show plots (interactive)", variable=self.show_plots_var)
        self.show_plots_check.grid(row=10, column=1, padx=5, pady=5)
        
        # Ylim
        tk.Label(self.run_frame, text="Ylim (Δ47 ymin,ymax):").grid(row=11, column=0, sticky=tk.W, padx=5, pady=5)
        self.ymin_var = tk.DoubleVar(value=0.15)
        self.ymin_entry = tk.Entry(self.run_frame, textvariable=self.ymin_var, width=10)
        self.ymin_entry.grid(row=11, column=1, padx=5, pady=5)
        self.ymax_var = tk.DoubleVar(value=0.68)
        self.ymax_entry = tk.Entry(self.run_frame, textvariable=self.ymax_var, width=10)
        self.ymax_entry.grid(row=11, column=2, padx=5, pady=5)
        
        # Right tick step
        tk.Label(self.run_frame, text="Right tick step (°C):").grid(row=12, column=0, sticky=tk.W, padx=5, pady=5)
        self.tick_step_var = tk.DoubleVar(value=50)
        self.tick_step_entry = tk.Entry(self.run_frame, textvariable=self.tick_step_var, width=10)
        self.tick_step_entry.grid(row=12, column=1, padx=5, pady=5)
        
        # Buttons
        self.run_button = tk.Button(self.run_frame, text="Run", command=self.run_analysis)
        self.run_button.grid(row=13, column=0, padx=5, pady=10)
        
        self.open_outdir_button = tk.Button(self.run_frame, text="Open Outdir", command=self.open_outdir)
        self.open_outdir_button.grid(row=13, column=1, padx=5, pady=10)
        
        self.quit_button = tk.Button(self.run_frame, text="Quit", command=self.root.destroy)
        self.quit_button.grid(row=13, column=2, padx=5, pady=10)
    
    def setup_uft_tab(self):
        # Implementation for Ufit tab
        pass
    
    def browse_thermal_csv(self):
        """Open file dialog to select thermal history CSV"""
        filename = filedialog.askopenfilename(
            title="Select Thermal History CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.thermal_csv_var.set(filename)
    
    def browse_test_csv(self):
        """Open file dialog to select test data CSV"""
        filename = filedialog.askopenfilename(
            title="Select Test Data CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.test_csv_var.set(filename)
    
    def run_analysis(self):
        """Run the analysis with current parameters"""
        try:
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
                messagebox.showerror("Error", f"Thermal history file not found: {thermal_path}")
                return
            
            if not os.path.exists(test_path):
                messagebox.showerror("Error", f"Test data file not found: {test_path}")
                return
            
            # Load data
            thermal_data = load_thermal_history(thermal_path, 
                                            time_col=self.time_col_var.get(),
                                            temp_col=self.avg_col_var.get())
            
            test_data = load_test_data(test_path,
                                     d47_col=self.d47_col_var.get(),
                                     sd_col=self.sd_col_var.get())
            
            # Perform analysis
            results = self.perform_analysis(thermal_data, test_data)
            
            # Save results
            outdir = self.outdir_var.get().strip()
            out_prefix = self.out_prefix_var.get().strip()
            
            if not os.path.isabs(outdir):
                outdir = os.path.join(os.getcwd(), outdir)
            
            os.makedirs(outdir, exist_ok=True)
            
            # Save plots
            plot_results(results, outdir, out_prefix, 
                        show_plots=self.show_plots_var.get(),
                        ymin=self.ymin_var.get(), ymax=self.ymax_var.get(),
                        tick_step=self.tick_step_var.get())
            
            messagebox.showinfo("Success", "Analysis completed successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during analysis: {str(e)}")
    
    def perform_analysis(self, thermal_data, test_data):
        """Perform the actual analysis"""
        # Implementation of analysis logic
        pass
    
    def open_outdir(self):
        """Open the output directory in file explorer"""
        outdir = self.outdir_var.get().strip()
        if not os.path.isabs(outdir):
            outdir = os.path.join(os.getcwd(), outdir)
        
        if os.path.exists(outdir):
            import subprocess
            subprocess.Popen(['explorer', outdir])
        else:
            messagebox.showerror("Error", f"Output directory does not exist: {outdir}")
    
    def set_default_values(self):
        """Set default values for all parameters"""
        # Default values are already set in the initialization
        pass

def main():
    app = App()
    app.mainloop()
