"""Track particles through the physical full ring and test the horizontal invariant Ix.

Ix may have been computed from one configured cell, several configured cells, or
the full ring.  The tracking test is deliberately different: one tracking turn
is always one complete physical 360-degree ring, supplied by FMA.py.

The invariant vector is stored in the scaled polynomial basis used by
nonlinear.py.  make_ix_callable() converts it to a numerical function of the
physical variables ordered as [delta, x, y, px, py].

The map metric is

    max_turn |Ix(turn) - Ix(0)| / (Ix_scale * turn)

where ``turn`` counts completed physical rings and ``Ix_scale`` is normally
|Ix(0)|, with a small configurable floor near the origin.  The plotted value is
log10 of this maximum relative drift per completed ring.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# 1. IX VECTOR -> CALLABLE FUNCTION
# =============================================================================


def make_ix_callable(Ix, state):
    """Return a vectorized callable Ix(delta, x, y, px, py)."""
    Ix = np.asarray(Ix, dtype=float)
    C = np.asarray(state["C"], dtype=float)
    idx_to_vec = state["idx_to_vec"]

    if len(Ix) != len(C):
        raise ValueError("Ix and state['C'] must have the same length.")

    terms = []
    for k, coefficient in enumerate(Ix * C):
        if coefficient == 0.0:
            continue
        powers = tuple(int(value) for value in idx_to_vec[k])
        terms.append((float(coefficient), powers))

    def evaluate(delta, x, y, px, py):
        arrays = np.broadcast_arrays(
            np.asarray(delta, dtype=float),
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(px, dtype=float),
            np.asarray(py, dtype=float),
        )
        result = np.zeros(arrays[0].shape, dtype=float)

        for coefficient, powers in terms:
            term = coefficient
            for values, power in zip(arrays, powers):
                if power:
                    term = term * values**power
            result = result + term

        return result

    return evaluate


# =============================================================================
# 2. SAME LAUNCH GRID USED BY PYAT FMA
# =============================================================================


def _linear_fma_grid(coords_mm, steps):
    """Reproduce the linear x-y launch grid used by fmap_parallel_track()."""
    if len(coords_mm) != 4:
        raise ValueError("coords_mm must contain [xmin, xmax, ymin, ymax].")
    if len(steps) != 2 or int(steps[0]) <= 0 or int(steps[1]) <= 0:
        raise ValueError("steps must contain two positive integers.")

    xmin = min(float(coords_mm[0]), float(coords_mm[1]))
    xmax = max(float(coords_mm[0]), float(coords_mm[1]))
    ymin = min(float(coords_mm[2]), float(coords_mm[3]))
    ymax = max(float(coords_mm[2]), float(coords_mm[3]))

    xsteps = int(steps[0])
    ysteps = int(steps[1])
    xstep = (xmax - xmin) / xsteps
    ystep = (ymax - ymin) / ysteps

    x_values = np.arange(xmin, xmax + 1.0e-6, xstep)
    y_values = np.arange(ymin, ymax + 1.0e-6, ystep)
    return x_values, y_values


def _initial_row(x_values_mm, y_value_mm, orbit, delta):
    """Build one y-row of initial coordinates exactly in PyAT coordinate order."""
    orbit = np.asarray(orbit, dtype=float).reshape(6)
    z0 = np.zeros((len(x_values_mm), 6), dtype=float)
    z0 += orbit
    z0[:, 4] += float(delta)

    # fmap_parallel_track() uses millimetres and adds 1 nm in both x and y
    # to avoid exact zero coordinates for ideal lattices.
    z0[:, 0] += 1.0e-3 * np.asarray(x_values_mm, dtype=float) + 1.0e-9
    z0[:, 2] += 1.0e-3 * float(y_value_mm) + 1.0e-9
    return z0


def _evaluate_at_coordinates(ix_callable, coordinates):
    """Evaluate Ix on PyAT coordinates [x, px, y, py, delta, ct]."""
    coordinates = np.asarray(coordinates, dtype=float)
    return ix_callable(
        coordinates[4],
        coordinates[0],
        coordinates[2],
        coordinates[1],
        coordinates[3],
    )


# =============================================================================
# 3. FULL-RING TRACKING AND IX DRIFT
# =============================================================================


def track_ix_invariance(
    ring,
    Ix,
    state,
    *,
    coords_mm,
    steps,
    nturns,
    orbit,
    delta=0.0,
    pool_size=None,
    normalization_floor_fraction=1.0e-12,
    verbose=True,
):
    """Track the FMA grid and measure Ix drift once per completed full ring."""
    try:
        from at.tracking import patpass
    except ImportError as exc:
        raise ImportError(
            'Install Accelerator Toolbox with: pip install "accelerator-toolbox[plot]"'
        ) from exc

    nturns = int(nturns)
    if nturns < 1:
        raise ValueError("nturns must be a positive integer.")

    normalization_floor_fraction = float(normalization_floor_fraction)
    if not np.isfinite(normalization_floor_fraction) or normalization_floor_fraction <= 0.0:
        raise ValueError("normalization_floor_fraction must be positive and finite.")

    ix_callable = make_ix_callable(Ix, state)
    x_values, y_values = _linear_fma_grid(coords_mm, steps)

    initial_rows = [
        _initial_row(x_values, y_value, orbit, delta)
        for y_value in y_values
    ]
    initial_ix_rows = [
        np.asarray(
            _evaluate_at_coordinates(ix_callable, row.T),
            dtype=float,
        ).reshape(-1)
        for row in initial_rows
    ]

    all_initial = np.concatenate(initial_ix_rows)
    reference_ix = float(np.max(np.abs(all_initial))) if all_initial.size else 0.0
    if reference_ix == 0.0 or not np.isfinite(reference_ix):
        reference_ix = 1.0
    normalization_floor = normalization_floor_fraction * reference_ix

    rows = []
    total_rows = len(y_values)

    for iy_index, (y_value, z0, ix0_row) in enumerate(
        zip(y_values, initial_rows, initial_ix_rows)
    ):
        if verbose:
            print(
                f"Ix tracking {100.0 * iy_index / max(total_rows, 1):.1f} %"
            )

        kwargs = {"losses": True}
        if pool_size is not None:
            kwargs["pool_size"] = int(pool_size)

        tracked, loss_map = patpass(
            ring,
            np.asfortranarray(z0.T),
            nturns,
            **kwargs,
        )
        tracked = np.asarray(tracked, dtype=float)

        if tracked.ndim != 4 or tracked.shape[0] != 6:
            raise RuntimeError(
                "Unexpected PyAT tracking output shape: "
                f"{tracked.shape}; expected (6, particles, refpts, turns)."
            )

        tracks = tracked[:, :, 0, :]
        lost_flags = np.asarray(
            loss_map.get("islost", np.zeros(len(x_values), dtype=bool)),
            dtype=bool,
        )

        for ix_index, x_value in enumerate(x_values):
            particle = tracks[:, ix_index, :]
            finite_turns = np.all(np.isfinite(particle), axis=0)

            if np.all(finite_turns):
                completed_turns = nturns
            else:
                bad = np.flatnonzero(~finite_turns)
                completed_turns = int(bad[0]) if len(bad) else nturns

            completed_turns = max(0, min(completed_turns, nturns))
            survived = bool(
                completed_turns == nturns
                and not bool(lost_flags[ix_index])
            )

            ix0 = float(ix0_row[ix_index])
            denominator = max(abs(ix0), normalization_floor)

            if completed_turns > 0:
                particle_valid = particle[:, :completed_turns]
                ix_turns = np.asarray(
                    _evaluate_at_coordinates(ix_callable, particle_valid),
                    dtype=float,
                ).reshape(-1)

                delta_ix = ix_turns - ix0
                relative = np.abs(delta_ix) / denominator
                turn_numbers = np.arange(1, completed_turns + 1, dtype=float)
                relative_rate = relative / turn_numbers

                ix_last = float(ix_turns[-1])
                max_abs_delta = float(np.max(np.abs(delta_ix)))
                max_relative = float(np.max(relative))
                rms_relative = float(np.sqrt(np.mean(relative**2)))
                max_relative_rate = float(np.max(relative_rate))
                log_rate = float(np.log10(max(max_relative_rate, 1.0e-300)))
            else:
                ix_last = np.nan
                max_abs_delta = np.nan
                max_relative = np.nan
                rms_relative = np.nan
                max_relative_rate = np.nan
                log_rate = np.nan

            rows.append([
                float(x_value),
                float(y_value),
                ix0,
                ix_last,
                max_abs_delta,
                max_relative,
                rms_relative,
                max_relative_rate,
                log_rate,
                float(completed_turns),
                float(survived),
            ])

    if verbose:
        print("Ix tracking 100.0 %")

    return np.asarray(rows, dtype=float)


# =============================================================================
# 4. SAVE AND PLOT
# =============================================================================


def save_ix_data(data, output_directory, case_label):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / f"{case_label}_Ix_invariance.csv"

    header = (
        "x_mm,y_mm,Ix_initial,Ix_last,max_abs_delta_Ix,"
        "max_relative_excursion,rms_relative_excursion,"
        "max_relative_drift_per_turn,log10_relative_drift_per_turn,"
        "completed_full_ring_turns,survived"
    )
    np.savetxt(path, data, delimiter=",", header=header, comments="")
    return path


def plot_ix_invariance(
    data,
    *,
    output_directory,
    case_label,
    native,
    coords_mm,
    delta,
    log_min,
    log_max,
    show_plot,
):
    """Plot log10 maximum relative Ix drift per completed physical ring."""
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / f"{case_label}_Ix_invariance_map.png"

    data = np.asarray(data, dtype=float)
    valid = (
        (data[:, 10] > 0.5)
        & np.isfinite(data[:, 8])
    )

    if not np.any(valid):
        raise RuntimeError(
            "No particle survived the full Ix-invariance tracking interval."
        )

    shown = data[valid]

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    scatter = ax.scatter(
        shown[:, 0],
        shown[:, 1],
        c=shown[:, 8],
        s=34,
        marker="s",
        vmin=float(log_min),
        vmax=float(log_max),
    )

    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label(
        r"$\log_{10}$ max relative $I_x$ drift / full-ring turn"
    )

    ax.set_xlabel(r"$x_0$ [mm]")
    ax.set_ylabel(r"$y_0$ [mm]")
    ax.set_title(
        "Horizontal Invariant Tracking Test\n"
        f"{case_label} | physical ring={native['n_cells']} cells | "
        f"delta={float(delta):g} | valid={len(shown)}/{len(data)}"
    )
    ax.set_xlim(float(coords_mm[0]), float(coords_mm[1]))
    ax.set_ylim(float(coords_mm[2]), float(coords_mm[3]))
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.20)

    fig.tight_layout()
    fig.savefig(path, dpi=220)

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return path


# =============================================================================
# 5. BUILD IX FROM THE CONFIGURED ANALYSIS LATTICE
# =============================================================================


def configured_ix(parameters=None):
    """Compute Ix with lattice_config.RING_NAMES, which may be only part of the ring."""
    import lattice_config as lattice_cfg
    import nonlinear_config as nonlinear_cfg
    import optimization_config as cfg
    import optimization as opt

    base_parameters = cfg.BASE_PARAMETERS if parameters is None else parameters

    context = opt.create_context(
        base_parameters,
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
    details = opt.full_diagnostics(
        context,
        cfg.GRADIENT_WEIGHT,
        nonlinear_cfg.LEAST_SQUARES_TOL,
    )
    return details["Ix"], context["state"], context, details


# =============================================================================
# 6. STANDALONE TEST
# =============================================================================


def main():
    import optimization_config as cfg
    import FMA
    from at.physics import find_orbit

    native = FMA.prepare_native_ring()
    ring = FMA.build_at_ring(native)
    orbit, _ = find_orbit(ring)

    Ix, state, _, _ = configured_ix(native["parameters"])

    data = track_ix_invariance(
        ring,
        Ix,
        state,
        coords_mm=cfg.FMA_COORDS_MM,
        steps=cfg.FMA_STEPS,
        nturns=2 * int(cfg.FMA_TURNS),
        orbit=orbit,
        delta=cfg.FMA_DELTA,
        pool_size=cfg.FMA_POOL_SIZE,
        normalization_floor_fraction=cfg.IX_INVARIANCE_NORM_FLOOR_FRACTION,
        verbose=True,
    )

    output_directory = cfg.FMA_OUTPUT_DIRECTORY / "Ix_test"
    data_path = save_ix_data(data, output_directory, cfg.FMA_CASE_LABEL)
    plot_path = plot_ix_invariance(
        data,
        output_directory=output_directory,
        case_label=cfg.FMA_CASE_LABEL,
        native=native,
        coords_mm=cfg.FMA_COORDS_MM,
        delta=cfg.FMA_DELTA,
        log_min=cfg.IX_INVARIANCE_LOG_MIN,
        log_max=cfg.IX_INVARIANCE_LOG_MAX,
        show_plot=cfg.FMA_SHOW_PLOT,
    )

    print(f"Ix data: {data_path}")
    print(f"Ix plot: {plot_path}")
    return data, data_path, plot_path


if __name__ == "__main__":
    main()
