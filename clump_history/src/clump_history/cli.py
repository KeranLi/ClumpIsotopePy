import argparse
from pathlib import Path

import isotopylog as ipl

from .io import load_thermal_history, load_test_data
from .fit import constrained_u_fit
from .model import compute_history
from .plot import plot_grid


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="clump-history",
        description="CLI for clumped isotope Δ47 forward modeling along thermal histories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # ---------- run ----------
    run = sub.add_parser("run", help="Run forward models for scenarios and plot a 2x3 grid.")
    run.add_argument("--thermal", type=Path, default=Path("./datasets/Thermal_History_Hu.csv"))
    run.add_argument("--test", type=Path, default=Path("./datasets/acutal_test_Hu.csv"))

    run.add_argument("--time-col", type=str, default="Time/Myr")
    run.add_argument("--avg-col", type=str, default="Avg_T/Celsius")
    run.add_argument("--d47-col", type=str, default="Delta47")
    run.add_argument("--sd-col", type=str, default="SD")

    run.add_argument("--mineral", type=str, default="calcite")
    run.add_argument("--reference", type=str, default="HH21")
    run.add_argument("--d0-std", type=float, default=0.02)

    run.add_argument("--peak-window", type=float, nargs=2, default=[550, 600], metavar=("START", "END"))
    run.add_argument("--peak-temps", type=float, nargs="+", default=[150, 200, 250, 300, 350],
                     help="Peak temperatures in °C for scenarios (initial + these).")
    run.add_argument("--no-initial", action="store_true", help="Do not include initial (unmodified) scenario.")

    run.add_argument("--ylim", type=float, nargs=2, default=[0.15, 0.68], metavar=("YMIN", "YMAX"))
    run.add_argument("--tick-step", type=float, default=50, help="Right-axis temperature tick step (°C).")

    run.add_argument("--out", type=Path, default=Path("output_figure_hu"),
                     help="Output prefix (no extension). Generates .pdf and .svg.")
    run.add_argument("--show", action="store_true", help="Show interactive window (otherwise just save).")

    # ---------- ufit (可选：只做温度曲线U型约束并输出) ----------
    ufit = sub.add_parser("ufit", help="Apply constrained U-fit to a thermal history and export a new CSV.")
    ufit.add_argument("--thermal", type=Path, default=Path("./datasets/Thermal_History_Hu.csv"))
    ufit.add_argument("--time-col", type=str, default="Time/Myr")
    ufit.add_argument("--avg-col", type=str, default="Avg_T/Celsius")

    ufit.add_argument("--peak-window", type=float, nargs=2, default=[550, 600], metavar=("START", "END"))
    ufit.add_argument("--peak-temp", type=float, required=True, help="Peak temperature (°C) for the U-fit.")
    ufit.add_argument("--out-csv", type=Path, default=Path("Thermal_History_adjusted.csv"))

    return p


def cmd_run(args: argparse.Namespace) -> None:
    time_myr, T_avg_k = load_thermal_history(args.thermal, args.time_col, args.avg_col)
    delta47, delta47_err = load_test_data(args.test, args.d47_col, args.sd_col)

    ed = ipl.EDistribution.from_literature(mineral=args.mineral, reference=args.reference)

    start_x, end_x = args.peak_window

    scenarios = []
    if not args.no_initial:
        D, Dstd, Deq = compute_history(time_myr, T_avg_k, ed, args.d0_std)
        scenarios.append(("initial", D, Dstd, Deq))

    for Tpeak_c in args.peak_temps:
        T_mod_k = constrained_u_fit(time_myr, T_avg_k, start_x, end_x, Tpeak_c + 273.15, plot=False)
        D, Dstd, Deq = compute_history(time_myr, T_mod_k, ed, args.d0_std)
        scenarios.append((f"{int(Tpeak_c)}", D, Dstd, Deq))

    # 只画前6个（2x3）
    scenarios = scenarios[:6]

    ymin, ymax = args.ylim
    plot_grid(
        time_myr=time_myr,
        scenarios=scenarios,
        delta47=delta47,
        delta47_err=delta47_err,
        out_prefix=args.out,
        ymin=ymin,
        ymax=ymax,
        tick_step_c=args.tick_step,
        show=args.show,
    )


def cmd_ufit(args: argparse.Namespace) -> None:
    import pandas as pd

    time_myr, T_avg_k = load_thermal_history(args.thermal, args.time_col, args.avg_col)
    start_x, end_x = args.peak_window

    T_new_k = constrained_u_fit(time_myr, T_avg_k, start_x, end_x, args.peak_temp + 273.15, plot=False)

    # 输出一个带新列的CSV（单位还是 Celsius 更方便）
    out = pd.DataFrame({
        args.time_col: time_myr,
        args.avg_col: T_new_k - 273.15,
    })
    out.to_csv(args.out_csv, index=False)
    print(f"[OK] Saved adjusted thermal history: {args.out_csv}")


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "ufit":
        cmd_ufit(args)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")
