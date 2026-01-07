import argparse
from pathlib import Path

import isotopylog as ipl

from .io import load_thermal_history, load_test_data
from .fit import constrained_u_fit
from .model import compute_history
from .plot import plot_grid
from . import __version__

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]   # ...\clump_history
WORKSPACE_ROOT = PROJECT_ROOT.parent                # ...\ClumpIsotope


def resolve_input_path(p: Path) -> Path:
    """Try absolute -> CWD-relative -> project-root-relative -> workspace-root-relative."""
    p = Path(p)
    if p.is_absolute():
        return p
    if p.exists():  # relative to current working dir
        return p
    alt1 = PROJECT_ROOT / p
    if alt1.exists():
        return alt1
    alt2 = WORKSPACE_ROOT / p
    if alt2.exists():
        return alt2
    return p


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="clump-history",
        description="CLI for clumped isotope Δ47 forward modeling along thermal histories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 顶层版本号
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    sub = p.add_subparsers(dest="cmd")
    # Python 3.7 没有 subparsers(required=True) 的兼容问题时也可用；
    # 但为了更友好，在 main() 里手动处理 cmd 为空的情况
    # 保留 required=True

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
    run.add_argument(
        "--peak-temps",
        type=float,
        nargs="+",
        default=[150, 200, 250, 300, 350],
        help="Peak temperatures in °C for scenarios (initial + these).",
    )
    run.add_argument("--no-initial", action="store_true", help="Do not include initial (unmodified) scenario.")

    run.add_argument("--ylim", type=float, nargs=2, default=[0.15, 0.68], metavar=("YMIN", "YMAX"))
    run.add_argument("--tick-step", type=float, default=50, help="Right-axis temperature tick step (°C).")

    # ---------- add out-fir ----------
    run.add_argument(
    "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory (will be created if missing).",
    )

    run.add_argument(
        "--out",
        type=Path,
        default=Path("output_figure_hu"),
        help="Output prefix (no extension). Generates .pdf and .svg.",
    )
    run.add_argument("--show", action="store_true", help="Show interactive window (otherwise just save).")

    # ---------- ufit ----------
    ufit = sub.add_parser("ufit", help="Apply constrained U-fit to a thermal history and export a new CSV.")
    ufit.add_argument("--thermal", type=Path, default=Path("./datasets/Thermal_History_Hu.csv"))
    ufit.add_argument("--time-col", type=str, default="Time/Myr")
    ufit.add_argument("--avg-col", type=str, default="Avg_T/Celsius")

    ufit.add_argument("--peak-window", type=float, nargs=2, default=[550, 600], metavar=("START", "END"))
    ufit.add_argument("--peak-temp", type=float, required=True, help="Peak temperature (°C) for the U-fit.")
    ufit.add_argument("--out-csv", type=Path, default=Path("Thermal_History_adjusted.csv"))
    ufit.add_argument("--outdir", type=Path, default=Path("."), help="Output directory (will be created if missing).")

    return p


def cmd_run(args: argparse.Namespace) -> None:

    if not args.thermal.exists():
        raise FileNotFoundError(
            f"Thermal file not found: {args.thermal}\n"
            f"Current working dir: {Path.cwd()}\n"
            f"Tip: run from workspace root (where ./datasets exists) "
            f"or pass --thermal with an absolute/relative path."
        )


    thermal_path = resolve_input_path(args.thermal)
    test_path = resolve_input_path(args.test)

    time_myr, T_avg_k = load_thermal_history(thermal_path, args.time_col, args.avg_col)
    delta47, delta47_err = load_test_data(test_path, args.d47_col, args.sd_col)

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

    # --------- output path handling ---------
    out_prefix = args.out

    # 如果用户给的 --out 没带目录（比如 fig_smoke），就放进 --outdir
    # 如果用户给的 --out 自带目录（比如 results/fig），则优先尊重它
    if out_prefix.parent == Path("."):
        out_prefix = args.outdir / out_prefix

    # 自动创建输出目录
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ymin, ymax = args.ylim
    plot_grid(
        time_myr=time_myr,
        scenarios=scenarios,
        delta47=delta47,
        delta47_err=delta47_err,
        out_prefix=out_prefix,   # contact path
        ymin=ymin,
        ymax=ymax,
        tick_step_c=args.tick_step,
        show=args.show,
    )


def cmd_ufit(args: argparse.Namespace) -> None:
    import pandas as pd

    if not args.thermal.exists():
        raise FileNotFoundError(
            f"Thermal file not found: {args.thermal}\n"
            f"Current working dir: {Path.cwd()}\n"
            f"Tip: run from workspace root (where ./datasets exists) "
            f"or pass --thermal with an absolute/relative path."
        )

    thermal_path = resolve_input_path(args.thermal)
    time_myr, T_avg_k = load_thermal_history(thermal_path, args.time_col, args.avg_col)
    start_x, end_x = args.peak_window

    T_new_k = constrained_u_fit(time_myr, T_avg_k, start_x, end_x, args.peak_temp + 273.15, plot=False)

    # 输出一个带新列的CSV（单位是 Celsius）
    out = pd.DataFrame({
        args.time_col: time_myr,
        args.avg_col: T_new_k - 273.15,
    })

    out_csv = args.out_csv
    if out_csv.parent == Path("."):
        out_csv = args.outdir / out_csv
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    out.to_csv(out_csv, index=False)
    print(f"[OK] Saved adjusted thermal history: {out_csv}")



def main():
    parser = build_parser()
    args = parser.parse_args()

    # 未提供子命令时，显示帮助
    if args.cmd is None:
        parser.print_help()
        raise SystemExit(2)

    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "ufit":
        cmd_ufit(args)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")

