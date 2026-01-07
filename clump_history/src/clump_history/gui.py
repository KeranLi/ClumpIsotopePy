import threading
from pathlib import Path
import os
import traceback

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import isotopylog as ipl

from .io import load_thermal_history, load_test_data
from .fit import constrained_u_fit
from .model import compute_history
from .plot import plot_grid


def _parse_float_list(s: str):
    s = (s or "").strip()
    if not s:
        return []
    parts = s.replace(",", " ").split()
    return [float(x) for x in parts]


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Clump History GUI")
        self.geometry("820x520")

        self._build_ui()

    def _build_ui(self):
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=10, pady=10)

        self.tab_run = ttk.Frame(nb)
        self.tab_ufit = ttk.Frame(nb)
        nb.add(self.tab_run, text="Run (2×3 scenarios)")
        nb.add(self.tab_ufit, text="Ufit (peak adjust)")

        self._build_run_tab(self.tab_run)
        self._build_ufit_tab(self.tab_ufit)

        # status bar
        self.status = tk.StringVar(value="Ready.")
        status_bar = ttk.Label(self, textvariable=self.status, anchor="w")
        status_bar.pack(fill="x", padx=10, pady=(0, 8))

    # ---------- RUN TAB ----------
    def _build_run_tab(self, parent):
        # variables
        self.run_thermal = tk.StringVar(value=str(Path("datasets") / "Thermal_History_Hu.csv"))
        self.run_test = tk.StringVar(value=str(Path("datasets") / "acutal_test_Hu.csv"))
        self.run_outdir = tk.StringVar(value="results")
        self.run_out = tk.StringVar(value="fig_smoke")

        self.run_time_col = tk.StringVar(value="Time/Myr")
        self.run_avg_col = tk.StringVar(value="Avg_T/Celsius")
        self.run_d47_col = tk.StringVar(value="Delta47")
        self.run_sd_col = tk.StringVar(value="SD")

        self.run_mineral = tk.StringVar(value="calcite")
        self.run_reference = tk.StringVar(value="HH21")
        self.run_d0_std = tk.StringVar(value="0.02")

        self.run_peak_start = tk.StringVar(value="550")
        self.run_peak_end = tk.StringVar(value="600")
        self.run_peak_temps = tk.StringVar(value="150 200 250 300 350")
        self.run_no_initial = tk.BooleanVar(value=False)

        self.run_ymin = tk.StringVar(value="0.15")
        self.run_ymax = tk.StringVar(value="0.68")
        self.run_tick_step = tk.StringVar(value="50")
        self.run_show = tk.BooleanVar(value=False)

        frm = ttk.Frame(parent)
        frm.pack(fill="both", expand=True)

        row = 0
        row = self._file_row(frm, row, "Thermal CSV:", self.run_thermal, filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        row = self._file_row(frm, row, "Test CSV:", self.run_test, filetypes=[("CSV", "*.csv"), ("All", "*.*")])

        row = self._text_row(frm, row, "Outdir:", self.run_outdir)
        row = self._text_row(frm, row, "Out prefix:", self.run_out)

        ttk.Separator(frm).grid(row=row, column=0, columnspan=4, sticky="ew", pady=8)
        row += 1

        row = self._text_row(frm, row, "Cols (time, avgT):", self.run_time_col, self.run_avg_col)
        row = self._text_row(frm, row, "Cols (Δ47, SD):", self.run_d47_col, self.run_sd_col)

        ttk.Separator(frm).grid(row=row, column=0, columnspan=4, sticky="ew", pady=8)
        row += 1

        row = self._text_row(frm, row, "Mineral / Ref:", self.run_mineral, self.run_reference)
        row = self._text_row(frm, row, "d0_std:", self.run_d0_std)

        row = self._text_row(frm, row, "Peak window (start,end Myr):", self.run_peak_start, self.run_peak_end)
        row = self._text_row(frm, row, "Peak temps (°C):", self.run_peak_temps)

        ttk.Checkbutton(frm, text="No initial scenario", variable=self.run_no_initial).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Checkbutton(frm, text="Show plots (interactive)", variable=self.run_show).grid(row=row, column=1, sticky="w", pady=4)
        row += 1

        row = self._text_row(frm, row, "Ylim (Δ47 ymin,ymax):", self.run_ymin, self.run_ymax)
        row = self._text_row(frm, row, "Right tick step (°C):", self.run_tick_step)

        btns = ttk.Frame(frm)
        btns.grid(row=row, column=0, columnspan=4, sticky="ew", pady=10)
        btns.columnconfigure(0, weight=1)
        btns.columnconfigure(1, weight=1)
        btns.columnconfigure(2, weight=1)

        self.btn_run = ttk.Button(btns, text="Run", command=self._on_run_click)
        self.btn_run.grid(row=0, column=0, sticky="ew", padx=5)

        ttk.Button(btns, text="Open Outdir", command=lambda: self._open_outdir(self.run_outdir.get())).grid(row=0, column=1, sticky="ew", padx=5)
        ttk.Button(btns, text="Quit", command=self.destroy).grid(row=0, column=2, sticky="ew", padx=5)

    # ---------- UFIT TAB ----------
    def _build_ufit_tab(self, parent):
        self.uf_thermal = tk.StringVar(value=str(Path("datasets") / "Thermal_History_Hu.csv"))
        self.uf_outdir = tk.StringVar(value="results")
        self.uf_outcsv = tk.StringVar(value="thermal_adjusted.csv")

        self.uf_time_col = tk.StringVar(value="Time/Myr")
        self.uf_avg_col = tk.StringVar(value="Avg_T/Celsius")

        self.uf_peak_start = tk.StringVar(value="550")
        self.uf_peak_end = tk.StringVar(value="600")
        self.uf_peak_temp = tk.StringVar(value="150")

        frm = ttk.Frame(parent)
        frm.pack(fill="both", expand=True)

        row = 0
        row = self._file_row(frm, row, "Thermal CSV:", self.uf_thermal, filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        row = self._text_row(frm, row, "Outdir:", self.uf_outdir)
        row = self._text_row(frm, row, "Out CSV name:", self.uf_outcsv)

        ttk.Separator(frm).grid(row=row, column=0, columnspan=4, sticky="ew", pady=8)
        row += 1

        row = self._text_row(frm, row, "Cols (time, avgT):", self.uf_time_col, self.uf_avg_col)
        row = self._text_row(frm, row, "Peak window (start,end Myr):", self.uf_peak_start, self.uf_peak_end)
        row = self._text_row(frm, row, "Peak temp (°C):", self.uf_peak_temp)

        btns = ttk.Frame(frm)
        btns.grid(row=row, column=0, columnspan=4, sticky="ew", pady=10)
        btns.columnconfigure(0, weight=1)
        btns.columnconfigure(1, weight=1)
        btns.columnconfigure(2, weight=1)

        self.btn_ufit = ttk.Button(btns, text="Ufit", command=self._on_ufit_click)
        self.btn_ufit.grid(row=0, column=0, sticky="ew", padx=5)

        ttk.Button(btns, text="Open Outdir", command=lambda: self._open_outdir(self.uf_outdir.get())).grid(row=0, column=1, sticky="ew", padx=5)
        ttk.Button(btns, text="Quit", command=self.destroy).grid(row=0, column=2, sticky="ew", padx=5)

    # ---------- helpers ----------
    def _file_row(self, parent, row, label, var, filetypes):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ent = ttk.Entry(parent, textvariable=var, width=80)
        ent.grid(row=row, column=1, columnspan=2, sticky="ew", padx=5)
        ttk.Button(parent, text="Browse", command=lambda: self._browse_file(var, filetypes)).grid(row=row, column=3, sticky="ew")
        parent.columnconfigure(1, weight=1)
        return row + 1

    def _text_row(self, parent, row, label, var1, var2=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=var1, width=40).grid(row=row, column=1, sticky="ew", padx=5)
        if var2 is not None:
            ttk.Entry(parent, textvariable=var2, width=40).grid(row=row, column=2, sticky="ew", padx=5)
        return row + 1

    def _browse_file(self, var, filetypes):
        p = filedialog.askopenfilename(filetypes=filetypes)
        if p:
            var.set(p)

    def _open_outdir(self, outdir):
        outdir = outdir.strip() or "."
        p = Path(outdir).resolve()
        _ensure_dir(p)
        try:
            os.startfile(str(p))  # Windows
        except Exception:
            messagebox.showinfo("Info", f"Output directory: {p}")

    # ---------- actions ----------
    def _on_run_click(self):
        self._set_busy(True)
        th = threading.Thread(target=self._run_job_wrapper, daemon=True)
        th.start()

    def _on_ufit_click(self):
        self._set_busy(True)
        th = threading.Thread(target=self._ufit_job_wrapper, daemon=True)
        th.start()

    def _set_busy(self, busy: bool):
        state = "disabled" if busy else "normal"
        self.btn_run.configure(state=state)
        self.btn_ufit.configure(state=state)
        self.status.set("Running..." if busy else "Ready.")

    def _run_job_wrapper(self):
        try:
            self._run_job()
            self._ui_ok("Run finished.")
        except Exception as e:
            self._ui_err(e)
        finally:
            self._ui_done()

    def _ufit_job_wrapper(self):
        try:
            self._ufit_job()
            self._ui_ok("Ufit finished.")
        except Exception as e:
            self._ui_err(e)
        finally:
            self._ui_done()

    def _ui_done(self):
        self.after(0, lambda: self._set_busy(False))

    def _ui_ok(self, msg):
        self.after(0, lambda: messagebox.showinfo("OK", msg))

    def _ui_err(self, e: Exception):
        tb = traceback.format_exc()
        self.after(0, lambda: messagebox.showerror("Error", f"{e}\n\n{tb}"))

    # ---------- core jobs ----------
    def _run_job(self):
        thermal = Path(self.run_thermal.get()).expanduser()
        test = Path(self.run_test.get()).expanduser()

        outdir = Path(self.run_outdir.get().strip() or ".")
        outname = Path(self.run_out.get().strip() or "output")
        _ensure_dir(outdir)

        out_prefix = outdir / outname

        start_x = float(self.run_peak_start.get())
        end_x = float(self.run_peak_end.get())
        peak_temps = _parse_float_list(self.run_peak_temps.get())

        time_col = self.run_time_col.get()
        avg_col = self.run_avg_col.get()
        d47_col = self.run_d47_col.get()
        sd_col = self.run_sd_col.get()

        mineral = self.run_mineral.get()
        reference = self.run_reference.get()
        d0_std = float(self.run_d0_std.get())

        ymin = float(self.run_ymin.get())
        ymax = float(self.run_ymax.get())
        tick_step = float(self.run_tick_step.get())

        no_initial = bool(self.run_no_initial.get())
        show = bool(self.run_show.get())

        time_myr, T_avg_k = load_thermal_history(thermal, time_col, avg_col)
        delta47, delta47_err = load_test_data(test, d47_col, sd_col)

        ed = ipl.EDistribution.from_literature(mineral=mineral, reference=reference)

        scenarios = []
        if not no_initial:
            D, Dstd, Deq = compute_history(time_myr, T_avg_k, ed, d0_std)
            scenarios.append(("initial", D, Dstd, Deq))

        for Tpeak_c in peak_temps:
            T_mod_k = constrained_u_fit(time_myr, T_avg_k, start_x, end_x, Tpeak_c + 273.15, plot=False)
            D, Dstd, Deq = compute_history(time_myr, T_mod_k, ed, d0_std)
            scenarios.append((f"{int(Tpeak_c)}", D, Dstd, Deq))

        scenarios = scenarios[:6]

        plot_grid(
            time_myr=time_myr,
            scenarios=scenarios,
            delta47=delta47,
            delta47_err=delta47_err,
            out_prefix=out_prefix,
            ymin=ymin,
            ymax=ymax,
            tick_step_c=tick_step,
            show=show,
        )

    def _ufit_job(self):
        import pandas as pd

        thermal = Path(self.uf_thermal.get()).expanduser()
        outdir = Path(self.uf_outdir.get().strip() or ".")
        outcsv = Path(self.uf_outcsv.get().strip() or "Thermal_History_adjusted.csv")
        _ensure_dir(outdir)
        out_path = outdir / outcsv

        time_col = self.uf_time_col.get()
        avg_col = self.uf_avg_col.get()

        start_x = float(self.uf_peak_start.get())
        end_x = float(self.uf_peak_end.get())
        peak_temp_c = float(self.uf_peak_temp.get())

        time_myr, T_avg_k = load_thermal_history(thermal, time_col, avg_col)
        T_new_k = constrained_u_fit(time_myr, T_avg_k, start_x, end_x, peak_temp_c + 273.15, plot=False)

        out = pd.DataFrame({
            time_col: time_myr,
            avg_col: T_new_k - 273.15,
        })
        out.to_csv(out_path, index=False)


def main():
    app = App()
    app.mainloop()
