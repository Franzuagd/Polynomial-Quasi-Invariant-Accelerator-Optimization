"""Settings for nonlinear optimization."""
from pathlib import Path
import lattice_config as lattice_cfg
import nonlinear_config as nonlinear_cfg

# 1. VARIABLES TO OPTIMIZE
VARY = list(lattice_cfg.VARY)     #magnets to edit names.
BASE_PARAMETERS = dict(lattice_cfg.PARAMETERS)   #parameters of the lattice (varables in the OPA file)   #sus it might be redundant

# 2. OBJECTIVE
GRADIENT_WEIGHT = 0.10    #how much influence is the derivate part gonna have
INVALID_PENALTY = 1.0e30    #???????????????????

# 3. CHROMATIC CORRECTION
CORRECT_CHROMATICITY = True   #Do you want it to be corrected? ofc you do, but just in case you got the choice.

# 4. OPTIMIZER SETTINGS
CMA_SIGMA = 0.50  #IDK check google
CMA_ITERS = 2#00    #
CMA_POPSIZE = 3#2   #
PRINT_EVERY = 1#0   
POWELL_ITERS = 2#00
POWELL_MAXFEV = 2#50
SCALES = {name: 100.0 for name in VARY}
CMA_TIME = 60.0
POWELL_TIME_FRACTION = 0.25

# 5. START / END SLICES  
PLOT_START_END_SLICES = True   #Do you want it to be plot at the beggining and at the end? 
SLICE_Y_VALUES = (0.0, 0.1 * nonlinear_cfg.PLOT_Y_MAX, 0.2*nonlinear_cfg.PLOT_Y_MAX)  #Amount of Ix slices to plot: (list of values of y)
SLICE_X_VALUES = (0.0) #Amount of Iy slices to plot: (list of values of x)
SLICE_DELTA_VALUES = (0.0, 0.5 * float(nonlinear_cfg.A_BOX[0]), float(nonlinear_cfg.A_BOX[0])) #  Different values of delta for the slices: (list of values of delta)
SLICE_FROZEN_MOMENTUM = 0.0 #IDK what this is for?

# 6. OUTPUTS
OUTPUT_ROOT = Path("optimization_output")
REPORT_FILE = OUTPUT_ROOT / "optimization_report.txt"
FINAL_LATTICE_FILE = OUTPUT_ROOT / "final_lattice.json"
PLOT_ROOT = OUTPUT_ROOT / "slices"

# 7. FREQUENCY MAP ANALYSIS
RUN_FMA_START_END = True

FMA_CASE_LABEL = "current_lattice"

# Frequency-map window. Accelerator Toolbox expects these coordinates in mm:
# [xmin, xmax, ymin, ymax].
FMA_COORDS_MM = [-15.0, 15.0, -15.0, 15.0]

# Number of initial conditions in x and y.
# Start small for testing; increase to [60, 60] or [100, 100] for final maps.
FMA_STEPS = [121, 121]

# PyAT tracks 2*TURNS internally: one block for the first tune estimate and
# one block for the second tune estimate.
FMA_TURNS = 256

# Number of integration slices used by PyAT for thick multipoles.
FMA_NUM_INT_STEPS = 10

# Use the same chromatic correction as the native lattice before conversion.
FMA_CORRECT_CHROMATICITY = CORRECT_CHROMATICITY

# The working lattice_config currently may use fewer cells for fast internal
# calculations. FMA should represent one physical turn of the machine.
# If None, infer the number of identical cells needed for 360 degrees of net bend.
FMA_PHYSICAL_RING_CELLS = None

# Optional parameter replacement. Leave empty to use BASE_PARAMETERS exactly
# when FMA.py is run directly. Optimization start/end diagnostics deliberately
# ignore these overrides so they test the exact start and final lattices.
# Example:
# FMA_PARAMETER_OVERRIDES = {"kse1": -30.1, "ko1": 5.0}
FMA_PARAMETER_OVERRIDES = {}

# Optional common momentum offset. 0.0 means on-momentum FMA.
FMA_DELTA = 0.0

# Parallel tracking. None lets PyAT choose its default process count.
FMA_POOL_SIZE = None

FMA_OUTPUT_DIRECTORY = OUTPUT_ROOT / "FMA"
FMA_SHOW_PLOT = True

# Ix invariance tracking uses the same launch grid and total tracking length as
# the FMA diagnostic (2 * FMA_TURNS complete physical-ring turns).
IX_INVARIANCE_NORM_FLOOR_FRACTION = 1.0e-12
IX_INVARIANCE_LOG_MIN = -14.0
IX_INVARIANCE_LOG_MAX = 0.0

