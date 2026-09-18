"""Executable nonlinear-model test, report, and plotting runner.

Dependencies expected in the same project:
    linear_lattice.py
    lattice_config.py
    nonlinear.py
    nonlinear_config.py

GLOBAL variables below are run choices only.  Accelerator globals remain in
lattice_config.py and nonlinear-model globals remain in nonlinear_config.py.
Everything created inside main() is local runtime state.
"""

from pathlib import Path
import io

import numpy as np
import sympy as sp

import linear_lattice as lin
import lattice_config as lcfg
import nonlinear as nl
import nonlinear_config as ncfg


# =============================================================================
# 1. GLOBAL RUN SETTINGS -- USER MAY EDIT THESE
# =============================================================================

RUN_CHECKS = True
RUN_PLOTS = True
PRINT_REPORT = True
PRINT_INVARIANT_POLYNOMIALS = True
PRINT_HAMILTONIAN = True

# Always generated.  The file is written in the current working directory.
REPORT_FILE = "nonlinear_report.txt"

# Horizontal shape objective:
#   J = sqrt(||Ix-Sx||_G^2 + lambda * ||grad_x,px(Ix-Sx)||_G^2)
# lambda=0 reproduces the old |Ix-Sx| objective.
GRADIENT_WEIGHT = 0.1

# Batch invariant slices.
# Ix: x-px plane with py=0 and y fixed to each value below.
# Iy: y-py plane with px=0 and x fixed to each value below.

PLOT_Y_SLICES = (0.0, 0.5e-3, 1.0e-3)
PLOT_X_SLICES = (0.0, 0.5e-3, 1.0e-3, 2e-3, 3e-3, 4e-3)
PLOT_DELTA_SLICES = (0.0, 0.005e-2)
PLOT_FROZEN_MOMENTUM = 0.0


# =============================================================================
# 2. REPORT HELPERS -- NOT USER CONFIGURATION
# =============================================================================


def _write_matrix(out, name, matrix):
    out.write(f"{name}\n")
    out.write(np.array2string(np.asarray(matrix), precision=10, suppress_small=False))
    out.write("\n\n")


def build_report(
    magnets,
    lattice,
    linear_data,
    correction,
    parameters,
    state,
    Sx,
    Sy,
    transfer,
    result,
    shape_objective,
    value_objective,
    gradient_objective,
    Ix_pol,
    Iy_pol,
    state_checks,
    transfer_check,
    plot_paths,
):
    out = io.StringIO()
    line = "=" * 88

    out.write(line + "\n")
    out.write("NONLINEAR MODEL TEST REPORT\n")
    out.write(line + "\n\n")

    out.write("ARCHITECTURE\n")
    out.write("linear_lattice.py    : general linear functions\n")
    out.write("lattice_config.py    : accelerator configuration\n")
    out.write("nonlinear.py         : general nonlinear functions\n")
    out.write("nonlinear_config.py  : nonlinear symbolic/numerical configuration\n")
    out.write("running_nonlinear.py : executable test/report/plot runner\n\n")

    out.write("LINEAR MODEL INPUT\n")
    out.write(f"Defined magnets       : {len(magnets)}\n")
    out.write(f"Elements in full ring : {len(lattice)}\n")
    out.write(f"Energy [GeV]          : {parameters[lcfg.ENERGY_PARAMETER]}\n")
    out.write(f"Tune x                : {lin.linear_data(linear_data, 'TUNE_X')}\n")
    out.write(f"Tune y                : {lin.linear_data(linear_data, 'TUNE_Y')}\n")
    if correction is not None:
        out.write(f"Corrected family 1    : {correction[0]} = {correction[1]}\n")
        out.write(f"Corrected family 2    : {correction[2]} = {correction[3]}\n")
    out.write("\n")

    out.write("NONLINEAR CONFIGURATION\n")
    out.write(f"Order m               : {ncfg.ORDER}\n")
    out.write(f"Delta order d         : {ncfg.DELTA_ORDER}\n")
    out.write(f"Canonical planes      : {ncfg.N_PLANES}\n")
    out.write(f"Variables             : {ncfg.VARIABLES}\n")
    out.write(f"Field symbols         : {ncfg.FIELD_SYMBOLS}\n")
    out.write(f"a_box                 : {np.asarray(ncfg.A_BOX)}\n")
    out.write(f"Least-squares tol     : {ncfg.LEAST_SQUARES_TOL}\n")
    out.write("Thin multipoles       : integrated kick, no artificial length\n")
    out.write(f"Cache repeated maps   : {ncfg.CACHE_REPEATED_MAGNET_MAPS}\n\n")

    out.write("SYMBOLIC HAMILTONIAN\n")
    out.write(str(sp.expand(ncfg.HAMILTONIAN)) + "\n\n")

    out.write("POLYNOMIAL BASIS\n")
    out.write(f"Total monomials       : {len(state['idx_to_vec'])}\n")
    out.write(f"Quadratic block size  : {state['quad_size']}\n")
    out.write(f"Nonquadratic size     : {state['nonquad_size']}\n")
    out.write(f"Hamiltonian terms kept: {len(state['H_dict'])}\n")
    out.write(f"Lie basis matrices    : {len(state['M_basis'])}\n")
    out.write(f"Poisson B matrices    : {len(state['B'])}\n\n")

    out.write("HAMILTONIAN BASIS TERMS\n")
    for idx in state["order"]:
        powers = state["idx_to_vec"][idx]
        coeff = state["H_dict"][idx]
        out.write(f"index={idx:5d}  powers={powers}  coefficient={coeff}\n")
    out.write("\n")

    out.write("COURANT-SNYDER QUADRATIC VECTORS\n")
    out.write(f"Sx = {np.array2string(Sx, precision=10)}\n")
    out.write(f"Sy = {np.array2string(Sy, precision=10)}\n\n")

    out.write("NONLINEAR TRANSFER\n")
    out.write(f"Transfer shape        : {transfer.shape}\n")
    q = state["quad_size"]
    out.write(f"Tqq shape             : {transfer[:q, :q].shape}\n")
    out.write(f"Tqn shape             : {transfer[:q, q:].shape}\n")
    out.write(f"Tnq shape             : {transfer[q:, :q].shape}\n")
    out.write(f"Tnn shape             : {transfer[q:, q:].shape}\n\n")

    out.write("INVARIANT / LEAST-SQUARES RESULTS\n")
    labels = [
        "||hx||_G", "||hy||_G", "err_x", "err_y",
        "||Ux||_G", "||Uy||_G", "relative err x",
        "relative err y", "||{Ix,Iy}||_G",
    ]
    for label, value in zip(labels, result[:9]):
        out.write(f"{label:22s}: {value:.16e}\n")
    out.write(f"Ix vector length      : {len(result[-2])}\n")
    out.write(f"Iy vector length      : {len(result[-1])}\n\n")

    out.write("HORIZONTAL SHAPE OBJECTIVE\n")
    out.write(f"Gradient weight lambda: {GRADIENT_WEIGHT}\n")
    out.write(f"||Ix-Sx||_G           : {value_objective:.16e}\n")
    out.write(f"gradient term         : {gradient_objective:.16e}\n")
    out.write(f"combined objective J  : {shape_objective:.16e}\n\n")

    out.write("INTERNAL STATE CHECKS\n")
    for key, value in state_checks.items():
        out.write(f"{key:30s}: {value}\n")
    out.write("\nTRANSFER CHECKS\n")
    for key, value in transfer_check.items():
        out.write(f"{key:30s}: {value}\n")
    out.write("\n")

    if PRINT_INVARIANT_POLYNOMIALS:
        out.write("HORIZONTAL INVARIANT POLYNOMIAL\n")
        out.write(str(Ix_pol) + "\n\n")
        out.write("VERTICAL INVARIANT POLYNOMIAL\n")
        out.write(str(Iy_pol) + "\n\n")

    if plot_paths:
        out.write("PLOTS\n")
        for label, path in plot_paths.items():
            out.write(f"{label}: {path}\n")
        out.write("\n")

    # Full transfer is intentionally last because it may be very large.
    _write_matrix(out, "FULL NONLINEAR TRANSFER MATRIX", transfer)
    out.write(line + "\n")
    return out.getvalue()


def save_report(report_text, filename):
    """Write the complete report to the current working directory."""
    report_path = Path.cwd() / filename

    with open(report_path, "w", encoding="utf-8") as file:
        file.write(report_text)

    if not report_path.is_file():
        raise RuntimeError(f"Report was not created: {report_path}")

    return report_path


# =============================================================================
# 3. MAIN -- ALL VARIABLES CREATED HERE ARE LOCAL RUNTIME VARIABLES
# =============================================================================


def main():
    # -------------------------------------------------------------------------
    # A. PREPARE THE EXISTING LINEAR MODEL
    # -------------------------------------------------------------------------
    magnets, lattice, linear_data, correction, parameters = lin.prepare_lattice(
        parameters=lcfg.PARAMETERS,
        ring_names=lcfg.RING_NAMES,
        magnet_builder=lcfg.define_magnets,
        energy_parameter=lcfg.ENERGY_PARAMETER,
        correction_parameter_map=lcfg.CORRECTION_PARAMETER_MAP,
        correct_chromatic=True,
        family1=lcfg.CHROMATIC_FAMILY1,
        family2=lcfg.CHROMATIC_FAMILY2,
        target_chrom_x=lcfg.TARGET_CHROM_X,
        target_chrom_y=lcfg.TARGET_CHROM_Y,
        repetitions=lcfg.REPETITIONS,
        step=lcfg.STEP,
    )

    # -------------------------------------------------------------------------
    # B. BUILD THE NONLINEAR SYMBOLIC/NUMERICAL STATE FROM nonlinear_config.py
    # -------------------------------------------------------------------------
    state = nl.initialize_nonlinear(
        data=linear_data,
        m=ncfg.ORDER,
        d=ncfg.DELTA_ORDER,
        hamiltonian=ncfg.HAMILTONIAN,
        a_box=ncfg.A_BOX,
        variables=ncfg.VARIABLES,
        field_symbols=ncfg.FIELD_SYMBOLS,
        n=ncfg.N_PLANES,
    )

    # -------------------------------------------------------------------------
    # C. TEST QUADRATIC PART + FULL NONLINEAR TRANSFER + INVARIANT SOLVE
    # -------------------------------------------------------------------------
    Sx, Sy = nl.quadratic_invariants(linear_data, state)

    Ix, Iy, result, transfer = nl.invariant_vectors(
        lattice=lattice,
        data=linear_data,
        state=state,
        tol=ncfg.LEAST_SQUARES_TOL,
        cache=ncfg.CACHE_REPEATED_MAGNET_MAPS,
    )

    # Fast horizontal shape objective.  This uses only coefficient-space
    # operations already prepared inside nonlinear.py; no grids or SymPy are
    # involved in the objective evaluation.
    shape_objective, value_objective, gradient_objective = (
        nl.horizontal_shape_objective(
            Ix=Ix,
            Sx=Sx,
            state=state,
            gradient_weight=GRADIENT_WEIGHT,
        )
    )

    Ix_pol = nl.vector_to_poly(np.asarray(Ix) * state["C"], state["monomial_basis"])
    Iy_pol = nl.vector_to_poly(np.asarray(Iy) * state["C"], state["monomial_basis"])

    # -------------------------------------------------------------------------
    # D. INTERNAL TESTS
    # -------------------------------------------------------------------------
    state_checks = nl.check_nonlinear_state(state) if RUN_CHECKS else {}
    transfer_check = nl.transfer_checks(transfer, state) if RUN_CHECKS else {}

    # -------------------------------------------------------------------------
    # E. OPTIONAL PLOTS
    # -------------------------------------------------------------------------
    # One timestamped folder is created per run. It contains 18 plots:
    #   Ix(x,px) for 3 fixed y values x 3 delta values
    #   Iy(y,py) for 3 fixed x values x 3 delta values
    # The inactive momentum is fixed to zero by default.
    plot_paths = {}
    plot_folder = None

    if RUN_PLOTS:
        plots = nl.plot_invariant_slices(
            Ix=Ix,
            Iy=Iy,
            state=state,
            y_values=PLOT_Y_SLICES,
            x_values=PLOT_X_SLICES,
            delta_values=PLOT_DELTA_SLICES,
            frozen_momentum=PLOT_FROZEN_MOMENTUM,
            levels=ncfg.PLOT_LEVELS,
            grid_points=ncfg.PLOT_GRID_POINTS,
            rmin=ncfg.PLOT_RMIN,
            rmax=ncfg.PLOT_RMAX,
            folder=ncfg.PLOT_FOLDER,
            x_max=ncfg.PLOT_X_MAX,
            px_max=ncfg.PLOT_PX_MAX,
            y_max=ncfg.PLOT_Y_MAX,
            py_max=ncfg.PLOT_PY_MAX,
        )

        plot_folder = plots["folder"]

        for delta0, y_slices in plots["Ix"].items():
            for y0, info in y_slices.items():
                label = f"Ix | y={y0:g}, py={PLOT_FROZEN_MOMENTUM:g}, delta={delta0:g}"
                plot_paths[label] = info[-1]

        for delta0, x_slices in plots["Iy"].items():
            for x0, info in x_slices.items():
                label = f"Iy | x={x0:g}, px={PLOT_FROZEN_MOMENTUM:g}, delta={delta0:g}"
                plot_paths[label] = info[-1]

    # -------------------------------------------------------------------------
    # F. COMPLETE TEXT REPORT -- ALWAYS GENERATED
    # -------------------------------------------------------------------------
    report_text = build_report(
        magnets,
        lattice,
        linear_data,
        correction,
        parameters,
        state,
        Sx,
        Sy,
        transfer,
        result,
        shape_objective,
        value_objective,
        gradient_objective,
        Ix_pol,
        Iy_pol,
        state_checks,
        transfer_check,
        plot_paths,
    )
    report_path = save_report(report_text, REPORT_FILE)

    if PRINT_REPORT:
        print(report_text)
    elif PRINT_HAMILTONIAN:
        print("Hamiltonian:")
        print(sp.expand(ncfg.HAMILTONIAN))

    print(f"Report saved to: {report_path}")
    print(
        "Horizontal shape objective: "
        f"J={shape_objective:.6e}  "
        f"value={value_objective:.6e}  "
        f"gradient={gradient_objective:.6e}"
    )
    if plot_folder is not None:
        print(f"Plot folder: {plot_folder}")
    for label, path in plot_paths.items():
        print(f"{label}: {path}")

    return {
        "magnets": magnets,
        "lattice": lattice,
        "linear_data": linear_data,
        "correction": correction,
        "parameters": parameters,
        "state": state,
        "Sx": Sx,
        "Sy": Sy,
        "Ix": Ix,
        "Iy": Iy,
        "Ix_pol": Ix_pol,
        "Iy_pol": Iy_pol,
        "transfer": transfer,
        "result": result,
        "shape_objective": shape_objective,
        "value_objective": value_objective,
        "gradient_objective": gradient_objective,
        "gradient_weight": GRADIENT_WEIGHT,
        "state_checks": state_checks,
        "transfer_checks": transfer_check,
        "plot_folder": plot_folder,
        "plot_paths": plot_paths,
        "report_path": report_path,
    }


if __name__ == "__main__":
    main()
