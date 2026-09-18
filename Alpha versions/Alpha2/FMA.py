"""Frequency Map Analysis (FMA) for the accelerator lattice in lattice_config.py.

This script:
    1. Reads the native lattice through lattice_config.py + linear_lattice.py.
    2. Applies the same chromatic correction used by the native model.
    3. Converts the native magnets to a PyAT ring.
    4. Runs Accelerator Toolbox fmap_parallel_track().
    5. Saves the numerical FMA data and a tune-diffusion map.
    6. Tracks the same launch grid to test turn-by-turn Ix invariance.

Required package:
    pip install "accelerator-toolbox[plot]"

Native magnet convention:
    [name, type, length, angle_deg, K, S, O, M, M5]

The native nonlinear Hamiltonian uses S and O directly as normal multipole
coefficients, so the PyAT mapping is:
    K -> PolynomB[1]
    S -> PolynomB[2]
    O -> PolynomB[3]

For zero-length multipoles, O is already an integrated strength and is passed
straight to at.ThinMultipole().
"""

from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np

import at
from at.physics.frequency_maps import fmap_parallel_track
from at.physics import find_orbit

import lattice_config as cfg
import optimization_config as opt_cfg
import linear_lattice as lin
import test_Ix as ix_test


# =============================================================================
# 1. USER SETTINGS FOR THE FMA
# =============================================================================
# All FMA user settings are defined in optimization_config.py.


# =============================================================================
# 2. BUILD THE PHYSICAL NATIVE RING
# =============================================================================


def _field(elem, name):
    return lin.magnet_field(elem, name)


def _infer_physical_ring_cells(parameters):
    """Infer the identical-cell repetition count from the net cell bend."""
    magnets = cfg.define_magnets(parameters)
    by_name = {_field(elem, "NAME"): elem for elem in magnets}

    cell_bend_deg = 0.0
    for name in cfg.CELL_NAMES:
        elem = by_name[name]
        if _field(elem, "TYPE") == "bending":
            cell_bend_deg += float(_field(elem, "ANGLE"))

    if abs(cell_bend_deg) < 1.0e-12:
        raise ValueError(
            "Cannot infer a full ring: the configured cell has zero net bending. "
            "Set PHYSICAL_RING_CELLS manually."
        )

    cells_float = 360.0 / abs(cell_bend_deg)
    cells = int(round(cells_float))

    if cells < 1 or not np.isclose(cells_float, cells, rtol=0.0, atol=1.0e-7):
        raise ValueError(
            f"The configured cell bends by {cell_bend_deg:.12g} deg, which does not "
            "give an integer number of cells for a 360-deg ring. "
            "Set PHYSICAL_RING_CELLS manually."
        )

    return cells, cell_bend_deg


def prepare_native_ring(parameters=None, apply_parameter_overrides=True):
    """Build the full native ring and apply the configured chromatic correction."""
    parameters = dict(opt_cfg.BASE_PARAMETERS if parameters is None else parameters)
    if apply_parameter_overrides:
        parameters.update(opt_cfg.FMA_PARAMETER_OVERRIDES)

    if opt_cfg.FMA_PHYSICAL_RING_CELLS is None:
        n_cells, cell_bend_deg = _infer_physical_ring_cells(parameters)
    else:
        n_cells = int(opt_cfg.FMA_PHYSICAL_RING_CELLS)
        if n_cells < 1:
            raise ValueError("PHYSICAL_RING_CELLS must be a positive integer.")
        _, cell_bend_deg = _infer_physical_ring_cells(parameters)

    ring_names = list(cfg.CELL_NAMES) * n_cells

    magnets, lattice, data, correction, parameters = lin.prepare_lattice(
        parameters=parameters,
        ring_names=ring_names,
        magnet_builder=cfg.define_magnets,
        energy_parameter=cfg.ENERGY_PARAMETER,
        correction_parameter_map=cfg.CORRECTION_PARAMETER_MAP,
        correct_chromatic=opt_cfg.FMA_CORRECT_CHROMATICITY,
        family1=cfg.CHROMATIC_FAMILY1,
        family2=cfg.CHROMATIC_FAMILY2,
        target_chrom_x=cfg.TARGET_CHROM_X,
        target_chrom_y=cfg.TARGET_CHROM_Y,
        repetitions=1,
        step=cfg.STEP,
    )

    total_bend_deg = sum(
        float(_field(elem, "ANGLE"))
        for elem in lattice
        if _field(elem, "TYPE") == "bending"
    )

    return {
        "magnets": magnets,
        "lattice": lattice,
        "data": data,
        "correction": correction,
        "parameters": parameters,
        "n_cells": n_cells,
        "cell_bend_deg": cell_bend_deg,
        "total_bend_deg": total_bend_deg,
    }


# =============================================================================
# 3. CONVERT THE NATIVE MAGNETS TO ACCELERATOR TOOLBOX
# =============================================================================


def native_element_to_at(elem):
    """Convert one native [name,type,L,angle,K,S,O,...] magnet to PyAT."""
    name = str(_field(elem, "NAME"))
    element_type = str(_field(elem, "TYPE")).lower()
    length = float(_field(elem, "LENGTH"))
    angle_deg = float(_field(elem, "ANGLE"))
    K = float(_field(elem, "K"))
    S = float(_field(elem, "S"))
    O = float(_field(elem, "O"))

    if element_type == "drift":
        return at.Drift(name, length)

    if element_type == "quadrupole":
        if abs(S) > 0.0 or abs(O) > 0.0:
            poly_b = np.array([0.0, K, S, O], dtype=float)
            return at.Multipole(
                name,
                length,
                np.zeros_like(poly_b),
                poly_b,
                NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS,
            )
        return at.Quadrupole(name, length, k=K, NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS)

    if element_type == "sextupole":
        if abs(K) > 0.0 or abs(O) > 0.0:
            poly_b = np.array([0.0, K, S, O], dtype=float)
            return at.Multipole(
                name,
                length,
                np.zeros_like(poly_b),
                poly_b,
                NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS,
            )
        # Native S is already the AT normal sextupole coefficient b_2.
        return at.Sextupole(name, length, h=S, NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS)

    if element_type == "multipole":
        poly_b = np.array([0.0, K, S, O], dtype=float)
        poly_a = np.zeros_like(poly_b)

        if length == 0.0:
            # K/S/O are integrated strengths for the native zero-length element.
            return at.ThinMultipole(name, poly_a, poly_b)

        return at.Multipole(
            name,
            length,
            poly_a,
            poly_b,
            NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS,
        )

    if element_type == "bending":
        if abs(S) > 0.0 or abs(O) > 0.0:
            raise NotImplementedError(
                f"Bending magnet {name} contains S={S} or O={O}. "
                "The present converter assumes bends contain only curvature and K."
            )

        return at.Dipole(
            name,
            length,
            bending_angle=math.radians(angle_deg),
            k=K,
            NumIntSteps=opt_cfg.FMA_NUM_INT_STEPS,
        )

    raise ValueError(f"Unknown native magnet type: {element_type!r} for {name}")


def build_at_ring(native):
    """Convert the prepared native full ring to a PyAT Lattice."""
    elements = [native_element_to_at(elem) for elem in native["lattice"]]

    # Native energy is stored in GeV; PyAT lattice energy is in eV.
    energy_eV = float(native["parameters"][cfg.ENERGY_PARAMETER]) * 1.0e9

    return at.Lattice(
        elements,
        name=f"FMA_{opt_cfg.FMA_CASE_LABEL}",
        energy=energy_eV,
    )


# =============================================================================
# 4. FREQUENCY MAP ANALYSIS
# =============================================================================


def run_fma(ring, orbit=None):
    """Run PyAT frequency-map analysis and return (fmap, lossmap)."""
    offset6d = np.zeros(6, dtype=float)
    offset6d[4] = float(opt_cfg.FMA_DELTA)

    kwargs = {}
    if opt_cfg.FMA_POOL_SIZE is not None:
        kwargs["pool_size"] = int(opt_cfg.FMA_POOL_SIZE)

    return fmap_parallel_track(
        ring,
        coords=list(opt_cfg.FMA_COORDS_MM),
        steps=list(opt_cfg.FMA_STEPS),
        turns=int(opt_cfg.FMA_TURNS),
        orbit=orbit,
        add_offset6D=offset6d,
        verbose=True,
        lossmap=True,
        **kwargs,
    )


# =============================================================================
# 5. SAVE AND PLOT
# =============================================================================


def save_fma_data(fmap, output_directory=None, case_label=None):
    output_directory = (
        opt_cfg.FMA_OUTPUT_DIRECTORY
        if output_directory is None
        else Path(output_directory)
    )
    case_label = opt_cfg.FMA_CASE_LABEL if case_label is None else str(case_label)
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / f"{case_label}_fma.csv"

    header = "x_mm,y_mm,nux,nuy,dnux,dnuy,log10_tune_diffusion"
    np.savetxt(path, fmap, delimiter=",", header=header, comments="")
    return path


def plot_frequency_map(fmap, native, output_directory=None, case_label=None):
    """Plot tune diffusion in the launched x-y plane."""
    output_directory = (
        opt_cfg.FMA_OUTPUT_DIRECTORY
        if output_directory is None
        else Path(output_directory)
    )
    case_label = opt_cfg.FMA_CASE_LABEL if case_label is None else str(case_label)
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / f"{case_label}_frequency_map.png"

    if fmap.size == 0:
        raise RuntimeError(
            "No particle survived with a valid frequency estimate. "
            "Reduce opt_cfg.FMA_COORDS_MM or inspect the lattice conversion."
        )

    x_values, y_values = ix_test._linear_fma_grid(
        opt_cfg.FMA_COORDS_MM,
        opt_cfg.FMA_STEPS,
    )
    total_launched = len(x_values) * len(y_values)
    valid = len(fmap)
    valid_fraction = valid / total_launched

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    scatter = ax.scatter(
        fmap[:, 0],
        fmap[:, 1],
        c=fmap[:, 6],
        s=34,
        marker="s",
        vmin=-10.0,
        vmax=-2.0,
    )

    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label(r"$\log_{10}$ tune diffusion")

    ax.set_xlabel(r"$x_0$ [mm]")
    ax.set_ylabel(r"$y_0$ [mm]")
    ax.set_title(
        "Frequency Map Analysis\n"
        f"{case_label} | {native['n_cells']} cells | "
        f"delta={opt_cfg.FMA_DELTA:g} | valid={valid}/{total_launched} ({100*valid_fraction:.1f}%)"
    )
    ax.set_xlim(opt_cfg.FMA_COORDS_MM[0], opt_cfg.FMA_COORDS_MM[1])
    ax.set_ylim(opt_cfg.FMA_COORDS_MM[2], opt_cfg.FMA_COORDS_MM[3])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.20)

    fig.tight_layout()
    fig.savefig(path, dpi=220)

    if opt_cfg.FMA_SHOW_PLOT:
        plt.show()
    else:
        plt.close(fig)

    return path


# =============================================================================
# 6. MAIN
# =============================================================================


def main(parameters=None, ix_vector=None, ix_state=None, stage=None):
    print("=" * 78)
    print("FREQUENCY MAP + IX INVARIANCE ANALYSIS")
    print("=" * 78)

    supplied_ix = ix_vector is not None or ix_state is not None
    if (ix_vector is None) != (ix_state is None):
        raise ValueError("ix_vector and ix_state must be supplied together.")

    native = prepare_native_ring(
        parameters=parameters,
        apply_parameter_overrides=not supplied_ix,
    )

    case_label = opt_cfg.FMA_CASE_LABEL
    output_directory = opt_cfg.FMA_OUTPUT_DIRECTORY
    if stage is not None:
        stage = str(stage)
        case_label = f"{case_label}_{stage}"
        output_directory = output_directory / stage

    print(f"Case                     : {case_label}")
    print(f"Physical ring cells      : {native['n_cells']}")
    print(f"Net bend / cell [deg]    : {native['cell_bend_deg']:.12g}")
    print(f"Total ring bend [deg]    : {native['total_bend_deg']:.12g}")
    print(f"Energy [GeV]             : {native['parameters'][cfg.ENERGY_PARAMETER]}")
    print(f"Native tune x            : {lin.linear_data(native['data'], 'TUNE_X')}")
    print(f"Native tune y            : {lin.linear_data(native['data'], 'TUNE_Y')}")

    if native["correction"] is not None:
        correction = native["correction"]
        print(
            f"Chromatic correction     : {correction[0]}={correction[1]:.10g}, "
            f"{correction[2]}={correction[3]:.10g}"
        )
        print(
            f"Corrected chromaticity   : ({correction[4]:.6g}, {correction[5]:.6g})"
        )

    if not np.isclose(abs(native["total_bend_deg"]), 360.0, atol=1.0e-6):
        raise ValueError(
            "FMA and Ix tracking must use one physical turn, but the converted "
            f"native lattice has total bend {native['total_bend_deg']:.12g} deg "
            "instead of 360 deg."
        )

    ring = build_at_ring(native)
    orbit, _ = find_orbit(ring)

    if ix_vector is None:
        ix_vector, ix_state, _, _ = ix_test.configured_ix(native["parameters"])

    print(f"PyAT elements            : {len(ring)}")
    print(f"FMA window [mm]          : {opt_cfg.FMA_COORDS_MM}")
    print(f"Grid                     : {opt_cfg.FMA_STEPS[0]} x {opt_cfg.FMA_STEPS[1]}")
    print(f"Tune-analysis turns      : {opt_cfg.FMA_TURNS} + {opt_cfg.FMA_TURNS}")
    print(f"Ix full-ring turns       : {2 * int(opt_cfg.FMA_TURNS)}")
    print(f"Momentum offset delta    : {opt_cfg.FMA_DELTA}")
    print("Running PyAT FMA...")

    fmap, losses = run_fma(ring, orbit=orbit)

    fma_data_path = save_fma_data(
        fmap,
        output_directory=output_directory,
        case_label=case_label,
    )
    frequency_plot_path = plot_frequency_map(
        fmap,
        native,
        output_directory=output_directory,
        case_label=case_label,
    )

    print("Running PyAT Ix-invariance tracking on the physical full ring...")
    ix_data = ix_test.track_ix_invariance(
        ring,
        ix_vector,
        ix_state,
        coords_mm=opt_cfg.FMA_COORDS_MM,
        steps=opt_cfg.FMA_STEPS,
        nturns=2 * int(opt_cfg.FMA_TURNS),
        orbit=orbit,
        delta=opt_cfg.FMA_DELTA,
        pool_size=opt_cfg.FMA_POOL_SIZE,
        normalization_floor_fraction=opt_cfg.IX_INVARIANCE_NORM_FLOOR_FRACTION,
        verbose=True,
    )
    ix_data_path = ix_test.save_ix_data(
        ix_data,
        output_directory,
        case_label,
    )
    ix_plot_path = ix_test.plot_ix_invariance(
        ix_data,
        output_directory=output_directory,
        case_label=case_label,
        native=native,
        coords_mm=opt_cfg.FMA_COORDS_MM,
        delta=opt_cfg.FMA_DELTA,
        log_min=opt_cfg.IX_INVARIANCE_LOG_MIN,
        log_max=opt_cfg.IX_INVARIANCE_LOG_MAX,
        show_plot=opt_cfg.FMA_SHOW_PLOT,
    )

    total_launched = len(ix_data)
    print("-" * 78)
    print(f"Valid FMA points         : {len(fmap)} / {total_launched}")
    if len(fmap):
        print(f"Median log diffusion     : {np.median(fmap[:, 6]):.6g}")
        print(f"Best log diffusion       : {np.min(fmap[:, 6]):.6g}")
        print(f"Worst log diffusion      : {np.max(fmap[:, 6]):.6g}")

    ix_survived = ix_data[:, 10] > 0.5
    ix_valid = ix_survived & np.isfinite(ix_data[:, 8])
    print(f"Valid Ix points          : {np.count_nonzero(ix_valid)} / {total_launched}")
    if np.any(ix_valid):
        print(f"Median log Ix drift      : {np.median(ix_data[ix_valid, 8]):.6g}")
        print(f"Best log Ix drift        : {np.min(ix_data[ix_valid, 8]):.6g}")
        print(f"Worst log Ix drift       : {np.max(ix_data[ix_valid, 8]):.6g}")

    print(f"FMA data                 : {fma_data_path}")
    print(f"FMA plot                 : {frequency_plot_path}")
    print(f"Ix data                  : {ix_data_path}")
    print(f"Ix plot                  : {ix_plot_path}")
    print("=" * 78)

    return {
        "fmap": fmap,
        "losses": losses,
        "ix_data": ix_data,
        "ring": ring,
        "native": native,
        "fma_data_path": fma_data_path,
        "frequency_plot_path": frequency_plot_path,
        "ix_data_path": ix_data_path,
        "ix_plot_path": ix_plot_path,
        "output_directory": output_directory,
        "case_label": case_label,
    }


if __name__ == "__main__":
    main()
