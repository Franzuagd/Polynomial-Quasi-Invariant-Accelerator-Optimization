"""Functions used by the nonlinear optimizer.

The lattice is prepared once.  Every candidate is applied with
linear_lattice.update_linear().  The nonlinear polynomial state is also built
once and reused during the optimization.
"""

from pathlib import Path
import json
import math

import numpy as np
from scipy.optimize import minimize

import linear_lattice as lin
import nonlinear as nl


# =============================================================================
# 1. SETUP / UPDATE
# =============================================================================


def initial_vector(vary, parameters):
    return np.asarray([parameters[name] for name in vary], dtype=float)


def create_context(
    parameters,
    *,
    ring_names,
    magnet_builder,
    energy_parameter,
    correction_parameter_map,
    correct_chromatic,
    family1,
    family2,
    target_chrom_x,
    target_chrom_y,
    repetitions,
    step,
    linear_variables,
    chromatic_variables,
    parameter_map,
    m,
    d,
    hamiltonian,
    a_box,
    variables,
    field_symbols,
    n_planes=2,
):
    magnets, lattice, data, correction, p = lin.prepare_lattice(
        parameters=parameters,
        ring_names=ring_names,
        magnet_builder=magnet_builder,
        energy_parameter=energy_parameter,
        correction_parameter_map=correction_parameter_map,
        correct_chromatic=correct_chromatic,
        family1=family1,
        family2=family2,
        target_chrom_x=target_chrom_x,
        target_chrom_y=target_chrom_y,
        repetitions=repetitions,
        step=step,
    )

    state = nl.initialize_nonlinear(
        data,
        m=m,
        d=d,
        hamiltonian=hamiltonian,
        a_box=a_box,
        variables=variables,
        field_symbols=field_symbols,
        n=n_planes,
    )

    return {
        "magnets": magnets,
        "lattice": lattice,
        "data": data,
        "correction": correction,
        "parameters": dict(p),
        "state": state,
        "map_cache": {},
        "settings": {
            "magnet_builder": magnet_builder,
            "energy_parameter": energy_parameter,
            "correction_parameter_map": correction_parameter_map,
            "correct_chromatic": correct_chromatic,
            "family1": family1,
            "family2": family2,
            "target_chrom_x": target_chrom_x,
            "target_chrom_y": target_chrom_y,
            "repetitions": repetitions,
            "step": step,
            "linear_variables": set(linear_variables),
            "chromatic_variables": set(chromatic_variables),
            "parameter_map": parameter_map,
        },
    }


def _refresh_nonlinear_normalization(state, data):
    """Update only the part of the nonlinear state that depends on CS0."""
    cs0 = np.asarray(lin.linear_data(data, "CS0"), dtype=float)
    bx0, ax0, gx0, _, _, _ = cs0

    vec_to_idx = state["vec_to_idx"]
    epsilon = np.asarray(state["epsilon"], dtype=float)

    ix2 = vec_to_idx[(0, 2, 0, 0, 0)]
    ixpx = vec_to_idx[(0, 1, 0, 1, 0)]
    ipx2 = vec_to_idx[(0, 0, 0, 2, 0)]

    arg = (
        bx0 * gx0 / (epsilon[ipx2] * epsilon[ix2])
        - ax0**2 / epsilon[ixpx]**2
    )
    if arg <= 0.0:
        raise ValueError("Nonlinear normalization is not positive.")

    state["C"] = epsilon * math.sqrt(arg)
    state["linear_cs0"] = cs0.copy()
    state["D_x"] = nl.build_derivative_matrix(state, 1)
    state["D_px"] = nl.build_derivative_matrix(state, 3)


def apply_candidate(context, v, vary):
    """Apply one optimizer vector using linear_lattice.update_linear()."""
    v = np.asarray(v, dtype=float)
    if len(v) != len(vary):
        raise ValueError(f"Expected {len(vary)} variables, received {len(v)}.")

    p = dict(context["parameters"])
    edited = []
    for name, value in zip(vary, v):
        value = float(value)
        if value != float(p[name]):
            edited.append(name)
        p[name] = value

    if not edited:
        return

    settings = context["settings"]
    linear_changed = bool(set(edited) & settings["linear_variables"])

    lattice, data, correction, p = lin.update_linear(
        lattice=context["lattice"],
        data=context["data"],
        parameters=p,
        edited_variables=edited,
        correct_chromatic=settings["correct_chromatic"],
        family1=settings["family1"],
        family2=settings["family2"],
        target_chrom_x=settings["target_chrom_x"],
        target_chrom_y=settings["target_chrom_y"],
        repetitions=settings["repetitions"],
        step=settings["step"],
        linear_variables=settings["linear_variables"],
        chromatic_variables=settings["chromatic_variables"],
        parameter_map=settings["parameter_map"],
        magnet_builder=settings["magnet_builder"],
        correction_parameter_map=settings["correction_parameter_map"],
        energy_parameter=settings["energy_parameter"],
    )

    context["lattice"] = lattice
    context["data"] = data
    context["parameters"] = dict(p)
    if correction is not None:
        context["correction"] = correction

    if linear_changed:
        _refresh_nonlinear_normalization(context["state"], data)


# =============================================================================
# 2. NONLINEAR TRANSFER / OBJECTIVE
# =============================================================================


def _magnet_signature(elem):
    return (
        lin.magnet_field(elem, "TYPE"),
        float(lin.magnet_field(elem, "LENGTH")),
        float(lin.magnet_field(elem, "ANGLE")),
        float(lin.magnet_field(elem, "K")),
        float(lin.magnet_field(elem, "S")),
        float(lin.magnet_field(elem, "O")),
    )


def nonlinear_transfer(context, tol):
    """Build the ring transfer, reusing maps of unchanged magnets."""
    state = context["state"]
    transfer = np.eye(len(state["idx_to_vec"]), dtype=float)
    cache = context["map_cache"]

    for elem in context["lattice"]:
        key = _magnet_signature(elem)
        tmatrix = cache.get(key)
        if tmatrix is None:
            tmatrix = nl.element_transfer(elem, state, tol=tol)[0]
            cache[key] = tmatrix
        transfer = tmatrix @ transfer

    q = state["quad_size"]
    return transfer, transfer[q:, q:], transfer[q:, :q]


def _horizontal_invariant(context, tol):
    state = context["state"]
    Sx, _ = nl.quadratic_invariants(context["data"], state)
    _, tnn, tnq = nonlinear_transfer(context, tol)

    Gnn = state["Gnn"]
    D = np.eye(state["nonquad_size"], dtype=float) - tnn
    Ux = tnq @ Sx

    Lg = np.linalg.cholesky(Gnn)
    hx, *_ = np.linalg.lstsq(Lg.T @ D, Lg.T @ Ux, rcond=tol)
    Ix = np.concatenate((Sx, hx))
    return Ix, Sx


def evaluate_candidate(context, v, vary, gradient_weight, tol, invalid_penalty):
    """Update one candidate and return the new horizontal objective."""
    try:
        apply_candidate(context, v, vary)
        Ix, Sx = _horizontal_invariant(context, tol)
        objective, value_norm, gradient_norm = nl.horizontal_shape_objective(
            Ix,
            Sx,
            context["state"],
            gradient_weight=gradient_weight,
        )
        if not np.isfinite(objective):
            return float(invalid_penalty)
        return float(objective)
    except (ValueError, KeyError, FloatingPointError, np.linalg.LinAlgError):
        return float(invalid_penalty)


def full_diagnostics(context, gradient_weight, tol):
    """Compute Ix and Iy for the start/end report and plots."""
    state = context["state"]
    Sx, Sy = nl.quadratic_invariants(context["data"], state)
    transfer, tnn, tnq = nonlinear_transfer(context, tol)
    result = nl.invariant(tnn, tnq, Sx, Sy, state, tol=tol)
    Ix, Iy = result[-2], result[-1]

    objective, value_norm, gradient_norm = nl.horizontal_shape_objective(
        Ix,
        Sx,
        state,
        gradient_weight=gradient_weight,
    )

    return {
        "objective": float(objective),
        "value_norm": float(value_norm),
        "gradient_norm": float(gradient_norm),
        "Ix": Ix,
        "Iy": Iy,
    }


# =============================================================================
# 3. PLOTS / SNAPSHOTS / SAVE
# =============================================================================


def plot_slices(details, state, folder, settings):
    return nl.plot_invariant_slices(
        details["Ix"],
        details["Iy"],
        state,
        y_values=settings["y_values"],
        x_values=settings["x_values"],
        delta_values=settings["delta_values"],
        frozen_momentum=settings["frozen_momentum"],
        levels=settings["levels"],
        grid_points=settings["grid_points"],
        rmin=settings["rmin"],
        rmax=settings["rmax"],
        folder=str(folder),
        x_max=settings["x_max"],
        px_max=settings["px_max"],
        y_max=settings["y_max"],
        py_max=settings["py_max"],
    )


def magnet_snapshot(lattice):
    result = {}
    for elem in lin.unique_magnets(lattice):
        name = str(lin.magnet_field(elem, "NAME"))
        result[name] = {
            "type": str(lin.magnet_field(elem, "TYPE")),
            "length": float(lin.magnet_field(elem, "LENGTH")),
            "angle": float(lin.magnet_field(elem, "ANGLE")),
            "K": float(lin.magnet_field(elem, "K")),
            "S": float(lin.magnet_field(elem, "S")),
            "O": float(lin.magnet_field(elem, "O")),
        }
    return result


def save_final_lattice(file_name, context):
    path = Path(file_name)
    path.parent.mkdir(parents=True, exist_ok=True)

    magnets = []
    for elem in lin.unique_magnets(context["lattice"]):
        magnets.append({
            "name": lin.magnet_field(elem, "NAME"),
            "type": lin.magnet_field(elem, "TYPE"),
            "length": float(lin.magnet_field(elem, "LENGTH")),
            "angle": float(lin.magnet_field(elem, "ANGLE")),
            "K": float(lin.magnet_field(elem, "K")),
            "S": float(lin.magnet_field(elem, "S")),
            "O": float(lin.magnet_field(elem, "O")),
        })

    payload = {
        "parameters": context["parameters"],
        "chromatic_correction": context["correction"],
        "magnets": magnets,
        "ring": [lin.magnet_field(elem, "NAME") for elem in context["lattice"]],
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


# =============================================================================
# 4. CMA-ES + POWELL
# =============================================================================


def hybrid_optimize(
    context,
    v0,
    vary,
    *,
    gradient_weight,
    tol,
    invalid_penalty,
    sigma,
    cma_iters,
    popsize,
    print_every,
    powell_iters,
    powell_maxfev,
    plot_start_end_slices,
    plot_root,
    slice_settings,
):
    try:
        import cma
    except ImportError as exc:
        raise ImportError("Install CMA-ES with: pip install cma") from exc

    v0 = np.asarray(v0, dtype=float)

    start_parameters = dict(context["parameters"])
    start_snapshot = magnet_snapshot(context["lattice"])
    start_correction = context["correction"]
    start_details = full_diagnostics(context, gradient_weight, tol)

    start_plots = None
    if plot_start_end_slices:
        start_plots = plot_slices(
            start_details,
            context["state"],
            Path(plot_root) / "start",
            slice_settings,
        )

    best_x = v0.copy()
    best_f = float(start_details["objective"])

    def objective(x):
        return evaluate_candidate(
            context,
            x,
            vary,
            gradient_weight,
            tol,
            invalid_penalty,
        )

    options = {"verb_disp": 0}
    if popsize is not None:
        options["popsize"] = int(popsize)

    es = cma.CMAEvolutionStrategy(v0, sigma, options)

    for iteration in range(1, int(cma_iters) + 1):
        solutions = es.ask()
        values = [objective(x) for x in solutions]
        es.tell(solutions, values)

        i = int(np.argmin(values))
        if values[i] < best_f:
            best_f = float(values[i])
            best_x = np.asarray(solutions[i], dtype=float).copy()

        if print_every and iteration % int(print_every) == 0:
            print(f"[CMA-ES] {iteration}/{cma_iters}   best J = {best_f:.6e}")

    powell_options = {"maxiter": int(powell_iters), "disp": False}
    if powell_maxfev is not None:
        powell_options["maxfev"] = int(powell_maxfev)

    res = minimize(objective, best_x, method="Powell", options=powell_options)

    if float(res.fun) < best_f:
        x_final = np.asarray(res.x, dtype=float)
    else:
        x_final = best_x

    apply_candidate(context, x_final, vary)
    final_details = full_diagnostics(context, gradient_weight, tol)
    final_parameters = dict(context["parameters"])
    final_snapshot = magnet_snapshot(context["lattice"])
    final_correction = context["correction"]

    final_plots = None
    if plot_start_end_slices:
        final_plots = plot_slices(
            final_details,
            context["state"],
            Path(plot_root) / "end",
            slice_settings,
        )

    return {
        "x_final": np.asarray(x_final, dtype=float),
        "f_final": float(final_details["objective"]),
        "start_details": start_details,
        "final_details": final_details,
        "start_parameters": start_parameters,
        "final_parameters": final_parameters,
        "start_snapshot": start_snapshot,
        "final_snapshot": final_snapshot,
        "start_correction": start_correction,
        "final_correction": final_correction,
        "start_plots": start_plots,
        "final_plots": final_plots,
        "context": context,
    }
