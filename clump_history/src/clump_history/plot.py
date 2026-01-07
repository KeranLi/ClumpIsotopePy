import numpy as np
import matplotlib.pyplot as plt
import isotopylog as ipl


def build_secondary_ticks(ymin, ymax, tick_step_c):
    d47_map = np.linspace(ymin, ymax, 800)
    T_map = np.array([ipl.T_from_Deq(d) - 273.15 for d in d47_map], dtype=float)

    ok = np.isfinite(T_map)
    d47_map = d47_map[ok]
    T_map = T_map[ok]

    Tmin, Tmax = float(np.min(T_map)), float(np.max(T_map))
    tick_start = np.floor(Tmin / tick_step_c) * tick_step_c
    tick_end = np.ceil(Tmax / tick_step_c) * tick_step_c
    T_ticks = np.arange(tick_start, tick_end + 0.1, tick_step_c)

    pairs = []
    for T in T_ticks:
        idx = int(np.abs(T_map - T).argmin())
        pairs.append((float(d47_map[idx]), int(round(T))))

    tmp = {}
    for d, T in pairs:
        tmp[d] = T
    pairs_sorted = sorted(tmp.items(), key=lambda x: x[0])

    d47_ticks = [d for d, _ in pairs_sorted]
    T_labels = [f"{T}" for _, T in pairs_sorted]
    return d47_ticks, T_labels


def plot_grid(time_myr, scenarios, delta47, delta47_err, out_prefix, ymin, ymax, tick_step_c, show=False):
    plt.rcParams["font.family"] = "Arial"

    d47_ticks, T_labels = build_secondary_ticks(ymin, ymax, tick_step_c)

    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True, sharey=True, constrained_layout=True)

    for i, (ax, (label, D, Dstd, Deq)) in enumerate(zip(axes.flatten(), scenarios)):
        ax.plot(time_myr, D, label=f'Forward-modeled at {label}°C')
        ax.fill_between(time_myr, D - Dstd, D + Dstd, alpha=0.35)
        ax.plot(time_myr, Deq, label=f'Equilibrium values at {label}°C')

        ax.errorbar(
            np.zeros_like(delta47),
            delta47, yerr=delta47_err,
            fmt='o',
            label='Actual $\Delta$47',
            color='black',
            capsize=4,
            alpha=0.6,
            markerfacecolor='none'
        )

        ax.set_ylim(ymin, ymax)

        if i >= 3:
            ax.set_xlabel('Age (Myr)')
        if i % 3 == 0:
            ax.set_ylabel('$\Delta$47 (‰)')

        secax = ax.secondary_yaxis('right')
        secax.set_yticks(d47_ticks)
        if i % 3 == 2:
            secax.set_yticklabels(T_labels)
            secax.set_ylabel('Temperature (°C)')
        else:
            secax.tick_params(labelright=False)

        xmax = float(np.max(time_myr))
        ax.set_xticks(np.arange(0, xmax + 1, 5), minor=True)
        ax.set_yticks(np.arange(ymin, ymax + 1e-9, 0.05), minor=True)
        ax.tick_params(axis='x', which='minor', length=4)
        ax.tick_params(axis='y', which='minor', length=4)

        ax.text(0.95, 0.05, f'({chr(97 + i)})', transform=ax.transAxes,
                fontsize=12, va='bottom', ha='right')

        ax.legend(loc='best', fontsize=9)

    svg_path = out_prefix.with_suffix(".svg")
    pdf_path = out_prefix.with_suffix(".pdf")
    fig.savefig(svg_path, format="svg")
    fig.savefig(pdf_path, format="pdf")

    if show:
        plt.show()
    else:
        plt.close(fig)

    print(f"[OK] Saved: {svg_path}")
    print(f"[OK] Saved: {pdf_path}")
