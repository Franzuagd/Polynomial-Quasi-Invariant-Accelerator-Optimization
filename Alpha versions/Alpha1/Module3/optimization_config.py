"""User-editable settings for nonlinear optimization."""

from pathlib import Path

import lattice_config as lattice_cfg
import nonlinear_config as nonlinear_cfg


# =============================================================================
# 1. VARIABLES TO OPTIMIZE
# =============================================================================

VARY = list(lattice_cfg.VARY)
BASE_PARAMETERS = dict(lattice_cfg.PARAMETERS)


# =============================================================================
# 2. OBJECTIVE
# =============================================================================

GRADIENT_WEIGHT = 0.10
INVALID_PENALTY = 1.0e30


# =============================================================================
# 3. CHROMATIC CORRECTION
# =============================================================================

CORRECT_CHROMATICITY = True


# =============================================================================
# 4. OPTIMIZER SETTINGS
# =============================================================================

CMA_SIGMA = 0.50

CMA_ITERS = 200
CMA_POPSIZE = 12
PRINT_EVERY = 30

POWELL_ITERS = 200
POWELL_MAXFEV = 250

# =============================================================================
# 5. START / END SLICES
# =============================================================================

PLOT_START_END_SLICES = True

SLICE_Y_VALUES = (
    0.0,
    0.5 * nonlinear_cfg.PLOT_Y_MAX,
    nonlinear_cfg.PLOT_Y_MAX,
)
SLICE_X_VALUES = (
    0.0,
    0.5 * nonlinear_cfg.PLOT_X_MAX,
    nonlinear_cfg.PLOT_X_MAX,
)
SLICE_DELTA_VALUES = (
    0.0,
    0.5 * float(nonlinear_cfg.A_BOX[0]),
    float(nonlinear_cfg.A_BOX[0]),
)
SLICE_FROZEN_MOMENTUM = 0.0


# =============================================================================
# 6. OUTPUTS
# =============================================================================

OUTPUT_ROOT = Path("optimization_output")
REPORT_FILE = OUTPUT_ROOT / "optimization_report.txt"
FINAL_LATTICE_FILE = OUTPUT_ROOT / "final_lattice.json"
PLOT_ROOT = OUTPUT_ROOT / "slices"
