"""Run the configured nonlinear optimization."""

from pathlib import Path

import numpy as np

import lattice_config as lattice_cfg
import nonlinear_config as nonlinear_cfg
import optimization_config as cfg
import optimization as opt


# =============================================================================
# 1. REPORT
# =============================================================================


def _f(value):
    return f"{float(value):.16e}"


def build_report(result):
    start_p = result["start_parameters"]
    final_p = result["final_parameters"]
    start_m = result["start_snapshot"]
    final_m = result["final_snapshot"]

    lines = []
    lines.append("=" * 90)
    lines.append("NONLINEAR OPTIMIZATION REPORT")
    lines.append("=" * 90)
    lines.append("")

    lines.append("OBJECTIVE")
    lines.append("-" * 90)
    lines.append(f"gradient weight = {cfg.GRADIENT_WEIGHT}")
    lines.append(f"initial J       = {_f(result['start_details']['objective'])}")
    lines.append(f"final J         = {_f(result['final_details']['objective'])}")
    lines.append("")

    lines.append("OPTIMIZED PARAMETERS")
    lines.append("-" * 90)
    lines.append(f"{'parameter':<14}{'initial':>24}{'final':>24}{'change':>24}")
    for name in cfg.VARY:
        a = float(start_p[name])
        b = float(final_p[name])
        lines.append(f"{name:<14}{_f(a):>24}{_f(b):>24}{_f(b-a):>24}")
    lines.append("")

    lines.append("MAGNET CHANGES")
    lines.append("-" * 90)
    lines.append(f"{'parameter':<12}{'magnet':<12}{'field':<10}{'initial':>22}{'final':>22}{'change':>22}")
    for parameter in cfg.VARY:
        for magnet, field in lattice_cfg.PARAMETER_MAP.get(parameter, []):
            key = {
                "LENGTH": "length",
                "ANGLE": "angle",
                "K": "K",
                "S": "S",
                "O": "O",
            }[field.upper()]
            a = start_m[magnet][key]
            b = final_m[magnet][key]
            lines.append(
                f"{parameter:<12}{magnet:<12}{field:<10}{_f(a):>22}{_f(b):>22}{_f(b-a):>22}"
            )
    lines.append("")

    lines.append("DIPOLE ANGLE CHANGES")
    lines.append("-" * 90)
    found = False
    for parameter in cfg.VARY:
        for magnet, field in lattice_cfg.PARAMETER_MAP.get(parameter, []):
            if field.upper() != "ANGLE":
                continue
            found = True
            a = start_m[magnet]["angle"]
            b = final_m[magnet]["angle"]
            lines.append(f"{parameter:<12}{magnet:<12}{_f(a):>22}{_f(b):>22}{_f(b-a):>22}")
    if not found:
        lines.append("No optimized variable changes a dipole angle.")
    lines.append("")

    for label, correction in (
        ("INITIAL CHROMATIC CORRECTION", result["start_correction"]),
        ("FINAL CHROMATIC CORRECTION", result["final_correction"]),
    ):
        lines.append(label)
        lines.append("-" * 90)
        if correction is None:
            lines.append("No chromatic correction.")
        else:
            lines.append(f"{correction[0]} = {_f(correction[1])}")
            lines.append(f"{correction[2]} = {_f(correction[3])}")
            lines.append(f"chrom_x = {_f(correction[4])}")
            lines.append(f"chrom_y = {_f(correction[5])}")
        lines.append("")

    lines.append("SLICE PLOTS")
    lines.append("-" * 90)
    if result["start_plots"] is not None:
        lines.append(f"start: {Path(result['start_plots']['folder']).resolve()}")
    if result["final_plots"] is not None:
        lines.append(f"end:   {Path(result['final_plots']['folder']).resolve()}")
    lines.append("=" * 90)

    return "\n".join(lines)


# =============================================================================
# 2. RUN
# =============================================================================


def main():
    cfg.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    context = opt.create_context(
        cfg.BASE_PARAMETERS,
        ring_names=lattice_cfg.RING_NAMES,
        magnet_builder=lattice_cfg.define_magnets,
        energy_parameter=lattice_cfg.ENERGY_PARAMETER,
        correction_parameter_map=lattice_cfg.CORRECTION_PARAMETER_MAP,
        correct_chromatic=cfg.CORRECT_CHROMATICITY,
        family1=lattice_cfg.CHROMATIC_FAMILY1,
        family2=lattice_cfg.CHROMATIC_FAMILY2,
        target_chrom_x=lattice_cfg.TARGET_CHROM_X,
        target_chrom_y=lattice_cfg.TARGET_CHROM_Y,
        repetitions=lattice_cfg.REPETITIONS,
        step=lattice_cfg.STEP,
        linear_variables=lattice_cfg.LINEAR_VARIABLES,
        chromatic_variables=lattice_cfg.CHROMATIC_VARIABLES,
        parameter_map=lattice_cfg.PARAMETER_MAP,
        m=nonlinear_cfg.ORDER,
        d=nonlinear_cfg.DELTA_ORDER,
        hamiltonian=nonlinear_cfg.HAMILTONIAN,
        a_box=nonlinear_cfg.A_BOX,
        variables=nonlinear_cfg.VARIABLES,
        field_symbols=nonlinear_cfg.FIELD_SYMBOLS,
        n_planes=nonlinear_cfg.N_PLANES,
    )

    v0 = opt.initial_vector(cfg.VARY, context["parameters"])

    slice_settings = {
        "y_values": cfg.SLICE_Y_VALUES,
        "x_values": cfg.SLICE_X_VALUES,
        "delta_values": cfg.SLICE_DELTA_VALUES,
        "frozen_momentum": cfg.SLICE_FROZEN_MOMENTUM,
        "levels": nonlinear_cfg.PLOT_LEVELS,
        "grid_points": nonlinear_cfg.PLOT_GRID_POINTS,
        "rmin": nonlinear_cfg.PLOT_RMIN,
        "rmax": nonlinear_cfg.PLOT_RMAX,
        "x_max": nonlinear_cfg.PLOT_X_MAX,
        "px_max": nonlinear_cfg.PLOT_PX_MAX,
        "y_max": nonlinear_cfg.PLOT_Y_MAX,
        "py_max": nonlinear_cfg.PLOT_PY_MAX,
    }

    print("=" * 80)
    print("STARTING NONLINEAR OPTIMIZATION")
    print("=" * 80)
    print("Varying:", cfg.VARY)
    print("Initial vector:", np.array2string(v0, precision=8))
    print("=" * 80)

    result = opt.hybrid_optimize(
        context,
        v0,
        cfg.VARY,
        gradient_weight=cfg.GRADIENT_WEIGHT,
        tol=nonlinear_cfg.LEAST_SQUARES_TOL,
        invalid_penalty=cfg.INVALID_PENALTY,
        sigma=cfg.CMA_SIGMA,
        cma_iters=cfg.CMA_ITERS,
        popsize=cfg.CMA_POPSIZE,
        print_every=cfg.PRINT_EVERY,
        powell_iters=cfg.POWELL_ITERS,
        powell_maxfev=cfg.POWELL_MAXFEV,
        plot_start_end_slices=cfg.PLOT_START_END_SLICES,
        plot_root=cfg.PLOT_ROOT,
        slice_settings=slice_settings,
    )

    final_path = opt.save_final_lattice(cfg.FINAL_LATTICE_FILE, result["context"])
    report = build_report(result)
    cfg.REPORT_FILE.write_text(report, encoding="utf-8")

    print("\n" + report)
    print(f"\nReport: {cfg.REPORT_FILE.resolve()}")
    print(f"Final lattice: {final_path.resolve()}")


if __name__ == "__main__":
    main()
