"""Executable runner for the configured linear lattice.

VARIABLE SCOPE RULE
-------------------
1) Everything in the section "GLOBAL RUN SETTINGS" is GLOBAL to this file.
   These are user choices controlling what this run will do.

2) The accelerator definition and physical model parameters are GLOBAL to
   lattice_config.py, not this file. They are accessed as cfg.<name>.

3) Everything created inside main() is LOCAL to main(). Runtime objects such
   as magnets, lattice, data, correction, parameters, and checks should remain
   local. Other modules should NOT import runtime objects from this file.
"""

from contextlib import redirect_stdout
from pathlib import Path
import io

import numpy as np

import linear_lattice as lin
import lattice_config as cfg


# =============================================================================
# 1. GLOBAL RUN SETTINGS -- USER MAY EDIT THESE
# =============================================================================
# These names are module-level GLOBAL variables because they control THIS run.
# They are not accelerator-physics parameters; those belong in lattice_config.py.

RUN_CHECKS = True
RUN_PLOTS = True
RUN_UPDATE_TEST = False

PRINT_PARAMETERS = True
PRINT_MAGNETS = True

# running.py ALWAYS writes a complete text report in the current working directory.
# Printing that same report in the terminal is optional.
REPORT_FILE = "linear_lattice_report.txt"
PRINT_REPORT = True

# Chromatic correction is a run choice. Families/targets are defined in cfg.
CORRECT_CHROMATICITY = True

# Only used when RUN_UPDATE_TEST = True.
# This dictionary is GLOBAL because it is a user-selected test configuration.
UPDATE_TEST_VALUES = {
    # Example: "ko1": 5.0,
    # Example: "X1": 3.60,
}


# =============================================================================
# 2. REPORT FUNCTION
# =============================================================================
# This is a helper function, not configuration. Every variable created inside
# this function is LOCAL to the function.


def print_full_report(
    magnets,
    lattice,
    data,
    correction,
    parameters,
    checks,
    update_result=None,
):
    """Print a detailed linear-lattice test report."""

    cs0 = lin.linear_data(data, "CS0")
    disp0 = lin.linear_data(data, "DISP0")
    M4 = lin.linear_data(data, "LATTICE_M4")
    M5 = lin.linear_data(data, "LATTICE_M5")
    s_values = lin.linear_data(data, "S_VALUES")
    cs_values = lin.linear_data(data, "CS_VALUES")
    disp_values = lin.linear_data(data, "DISP_VALUES")

    print("=" * 80)
    print("LINEAR LATTICE TEST REPORT")
    print("=" * 80)

    print("\nGENERAL")
    print("-" * 80)
    print(f"Defined magnets       : {len(magnets)}")
    print("Scope                  : full ring")
    print(f"Analyzed elements      : {len(lattice)}")
    print(f"Elements in one cell   : {len(cfg.CELL_NAMES)}")
    print(f"Number of cells        : {cfg.N_CELLS}")
    print(f"Elements in full ring  : {len(cfg.RING_NAMES)}")
    print(f"Energy [GeV]           : {parameters[cfg.ENERGY_PARAMETER]}")
    print(f"Integration step [m]   : {cfg.STEP}")
    print(f"Repetitions            : {cfg.REPETITIONS}")

    print("\nCURRENT PARAMETERS")
    print("-" * 80)
    for name, value in parameters.items():
        print(f"{name:12s} = {value}")

    print("\nVARIABLE GROUPS")
    print("-" * 80)
    print(f"VARY                  = {cfg.VARY}")
    print(f"LINEAR_VARIABLES      = {sorted(cfg.LINEAR_VARIABLES)}")
    print(f"CHROMATIC_VARIABLES   = {sorted(cfg.CHROMATIC_VARIABLES)}")
    print(f"NONLINEAR_VARIABLES   = {sorted(cfg.NONLINEAR_VARIABLES)}")

    print("\nPERIODIC TWISS AT s = 0")
    print("-" * 80)
    print(f"beta_x  = {cs0[0]:.12g}")
    print(f"alpha_x = {cs0[1]:.12g}")
    print(f"gamma_x = {cs0[2]:.12g}")
    print(f"beta_y  = {cs0[3]:.12g}")
    print(f"alpha_y = {cs0[4]:.12g}")
    print(f"gamma_y = {cs0[5]:.12g}")

    print("\nPERIODIC DISPERSION AT s = 0")
    print("-" * 80)
    print(f"Dx      = {disp0[0]:.12g}")
    print(f"Dpx     = {disp0[1]:.12g}")
    print(f"Dy      = {disp0[2]:.12g}")
    print(f"Dpy     = {disp0[3]:.12g}")
    print(f"delta   = {disp0[4]:.12g}")

    print("\nRING OPTICS")
    print("-" * 80)
    print(f"Tune x                 = {lin.linear_data(data, 'TUNE_X'):.12g}")
    print(f"Tune y                 = {lin.linear_data(data, 'TUNE_Y'):.12g}")
    print(f"Natural chromaticity x = {lin.linear_data(data, 'CHROM_X'):.12g}")
    print(f"Natural chromaticity y = {lin.linear_data(data, 'CHROM_Y'):.12g}")
    print(f"Natural emittance      = {lin.linear_data(data, 'EMITTANCE'):.12g}")
    print(f"Circumference [m]      = {lin.linear_data(data, 'CIRCUMFERENCE'):.12g}")

    print("\nCHROMATIC CORRECTION")
    print("-" * 80)
    if correction is None:
        print("Chromatic correction disabled or not performed.")
    else:
        print(f"{correction[0]:8s} = {correction[1]:.12g}")
        print(f"{correction[2]:8s} = {correction[3]:.12g}")
        print(f"Corrected chromaticity x = {correction[4]:.12g}")
        print(f"Corrected chromaticity y = {correction[5]:.12g}")

    print("\nFULL-RING 4x4 MATRIX")
    print("-" * 80)
    print(np.array2string(M4, precision=12, suppress_small=False))

    print("\nFULL-RING 5x5 MATRIX")
    print("-" * 80)
    print(np.array2string(M5, precision=12, suppress_small=False))

    print("\nCONSISTENCY CHECKS")
    print("-" * 80)
    if checks is None:
        print("Consistency checks disabled.")
    else:
        print(f"Symplectic error       = {checks[0]:.12e}")
        print(f"Twiss closure          = {checks[1]:.12e}")
        print(f"Dispersion closure     = {checks[2]:.12e}")
        print(f"CS identity x error    = {checks[3]:.12e}")
        print(f"CS identity y error    = {checks[4]:.12e}")
        print(f"Total bending [deg]    = {checks[5]:.12f}")

    print("\nSAMPLING")
    print("-" * 80)
    print(f"Number of s samples    = {len(s_values)}")
    print(f"s range [m]            = {s_values[0]:.12g} -> {s_values[-1]:.12g}")
    print(f"s_values shape         = {np.shape(s_values)}")
    print(f"cs_values shape        = {np.shape(cs_values)}")
    print(f"disp_values shape      = {np.shape(disp_values)}")

    print("\nMAGNET TABLE")
    print("-" * 80)
    print(
        f"{'NAME':<8} {'TYPE':<12} {'LENGTH':>12} {'ANGLE':>12} "
        f"{'K':>16} {'S':>16} {'O':>16}"
    )
    print("-" * 96)
    for elem in magnets:
        print(
            f"{lin.magnet_field(elem, 'NAME'):<8} "
            f"{lin.magnet_field(elem, 'TYPE'):<12} "
            f"{lin.magnet_field(elem, 'LENGTH'):>12.6g} "
            f"{lin.magnet_field(elem, 'ANGLE'):>12.6g} "
            f"{lin.magnet_field(elem, 'K'):>16.9g} "
            f"{lin.magnet_field(elem, 'S'):>16.9g} "
            f"{lin.magnet_field(elem, 'O'):>16.9g}"
        )

    print("\nUPDATE_LINEAR TEST")
    print("-" * 80)
    if update_result is None:
        print("Update test not requested.")
    else:
        edited_variables, test_parameters, update_data, update_correction, update_checks = update_result
        print(f"Edited variables       = {edited_variables}")
        for name in edited_variables:
            print(f"{name:12s} = {test_parameters[name]}")
        print(f"Tune x                 = {lin.linear_data(update_data, 'TUNE_X'):.12g}")
        print(f"Tune y                 = {lin.linear_data(update_data, 'TUNE_Y'):.12g}")
        print(f"Natural chromaticity x = {lin.linear_data(update_data, 'CHROM_X'):.12g}")
        print(f"Natural chromaticity y = {lin.linear_data(update_data, 'CHROM_Y'):.12g}")
        if update_correction is not None:
            print(f"{update_correction[0]:8s} = {update_correction[1]:.12g}")
            print(f"{update_correction[2]:8s} = {update_correction[3]:.12g}")
            print(f"Corrected chromaticity x = {update_correction[4]:.12g}")
            print(f"Corrected chromaticity y = {update_correction[5]:.12g}")
        print(f"Symplectic error       = {update_checks[0]:.12e}")
        print(f"Twiss closure          = {update_checks[1]:.12e}")
        print(f"Dispersion closure     = {update_checks[2]:.12e}")

    print("\n" + "=" * 80)
    print("LINEAR LATTICE TEST COMPLETED SUCCESSFULLY")
    print("=" * 80)


def save_report(
    magnets,
    lattice,
    data,
    correction,
    parameters,
    checks,
    update_result=None,
    filename=REPORT_FILE,
):
    """Generate and save the complete linear-lattice report as a text file."""

    report_buffer = io.StringIO()

    with redirect_stdout(report_buffer):
        print_full_report(
            magnets=magnets,
            lattice=lattice,
            data=data,
            correction=correction,
            parameters=parameters,
            checks=checks,
            update_result=update_result,
        )

    report_text = report_buffer.getvalue()

    # =========================================================================
    # ACTUAL TEXT-FILE CREATION AND SAVING
    # =========================================================================
    # Path.cwd() is the folder from which you run:  python running.py
    report_path = Path.cwd() / filename

    with open(report_path, "w", encoding="utf-8") as file:
        file.write(report_text)

    # Verify that the file really exists on disk.
    if not report_path.is_file():
        raise RuntimeError(f"Report was not created: {report_path}")

    return report_path.resolve(), report_text


# =============================================================================
# 3. MAIN EXECUTION
# =============================================================================


def main():
    # -------------------------------------------------------------------------
    # LOCAL RUNTIME VARIABLES
    # -------------------------------------------------------------------------
    # Every variable assigned below is LOCAL to main().
    # Do not import these runtime variables from running.py in nonlinear/optimization.
    # If another module needs the model, it should call lin.prepare_lattice(...)
    # using the GLOBAL machine definitions in lattice_config.py.

    magnets, lattice, data, correction, parameters = lin.prepare_lattice(
        parameters=cfg.PARAMETERS,
        ring_names=cfg.RING_NAMES,
        magnet_builder=cfg.define_magnets,
        energy_parameter=cfg.ENERGY_PARAMETER,
        correction_parameter_map=cfg.CORRECTION_PARAMETER_MAP,
        correct_chromatic=CORRECT_CHROMATICITY,
        family1=cfg.CHROMATIC_FAMILY1,
        family2=cfg.CHROMATIC_FAMILY2,
        target_chrom_x=cfg.TARGET_CHROM_X,
        target_chrom_y=cfg.TARGET_CHROM_Y,
        repetitions=cfg.REPETITIONS,
        step=cfg.STEP,
    )

    checks = lin.check_linear_lattice(lattice, data) if RUN_CHECKS else None
    update_result = None

    # -------------------------------------------------------------------------
    # OPTIONAL UPDATE_LINEAR TEST
    # -------------------------------------------------------------------------

    if RUN_UPDATE_TEST:
        if not UPDATE_TEST_VALUES:
            print("UPDATE TEST SKIPPED: UPDATE_TEST_VALUES is empty.")
        else:
            test_parameters = parameters.copy()
            test_parameters.update(UPDATE_TEST_VALUES)
            edited_variables = list(UPDATE_TEST_VALUES)

            lattice, update_data, update_correction, test_parameters = lin.update_linear(
                lattice=lattice,
                data=data,
                parameters=test_parameters,
                edited_variables=edited_variables,
                correct_chromatic=CORRECT_CHROMATICITY,
                family1=cfg.CHROMATIC_FAMILY1,
                family2=cfg.CHROMATIC_FAMILY2,
                target_chrom_x=cfg.TARGET_CHROM_X,
                target_chrom_y=cfg.TARGET_CHROM_Y,
                repetitions=cfg.REPETITIONS,
                step=cfg.STEP,
                linear_variables=cfg.LINEAR_VARIABLES,
                chromatic_variables=cfg.CHROMATIC_VARIABLES,
                parameter_map=cfg.PARAMETER_MAP,
                magnet_builder=cfg.define_magnets,
                correction_parameter_map=cfg.CORRECTION_PARAMETER_MAP,
                energy_parameter=cfg.ENERGY_PARAMETER,
            )

            update_checks = lin.check_linear_lattice(lattice, update_data)
            update_result = [
                edited_variables,
                test_parameters,
                update_data,
                update_correction,
                update_checks,
            ]

    # -------------------------------------------------------------------------
    # OPTIONAL SIMPLE TERMINAL DUMPS
    # -------------------------------------------------------------------------

    if PRINT_PARAMETERS:
        print("=" * 80)
        print("CURRENT PARAMETERS")
        print("=" * 80)
        for name, value in parameters.items():
            print(f"{name:12s} = {value}")
        print()

    if PRINT_MAGNETS:
        print("=" * 80)
        print("DEFINED MAGNETS")
        print("=" * 80)
        for elem in magnets:
            print(elem)
        print()

    # -------------------------------------------------------------------------
    # COMPLETE TEXT REPORT -- ALWAYS GENERATED
    # -------------------------------------------------------------------------
    # The report file is a standard output of running.py. It is generated on
    # every run regardless of the optional terminal/plot/test choices above.

    report_path, report_text = save_report(
        magnets=magnets,
        lattice=lattice,
        data=data,
        correction=correction,
        parameters=parameters,
        checks=checks,
        update_result=update_result,
        filename=REPORT_FILE,
    )

    if PRINT_REPORT:
        print(report_text, end="")

    print(f"\nReport saved to: {report_path}")

    # -------------------------------------------------------------------------
    # PLOTS
    # -------------------------------------------------------------------------

    if RUN_PLOTS:
        lin.plot_linear_functions(data)


if __name__ == "__main__":
    main()
