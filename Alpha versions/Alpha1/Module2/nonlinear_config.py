"""User-editable configuration for the nonlinear polynomial model.

This file contains ONLY nonlinear-model choices.  The accelerator itself remains
in lattice_config.py.  The symbolic Hamiltonian is intentionally defined here so
that changing the physical nonlinear model never requires editing nonlinear.py.
"""

import numpy as np
import sympy as sp


# =============================================================================
# 1. SYMBOLIC PHASE-SPACE VARIABLES
# =============================================================================
# Keep this order unless nonlinear.py is generalized beyond the current model:
#     [delta, x, y, px, py]

delta, x, y, px, py = sp.symbols("delta x y px py")
VARIABLES = [delta, x, y, px, py]


# =============================================================================
# 2. SYMBOLIC ELEMENT / HAMILTONIAN COEFFICIENTS
# =============================================================================
# nonlinear.py feeds these symbols, in this order, with:
#     b1 = curvature
#     b2 = quadrupole strength K
#     b3 = sextupole strength S
#     b4 = nonlinear multipole strength O
#     b5 = reserved fifth coefficient (currently 0 for every lattice element)

b1, b2, b3, b4, b5 = sp.symbols("b1 b2 b3 b4 b5")
FIELD_SYMBOLS = [b1, b2, b3, b4, b5]


# =============================================================================
# 3. USER-EDITABLE HAMILTONIAN
# =============================================================================
# IMPORTANT:
# This polynomial defines the symbolic nonlinear model used to build the Lie
# matrices.  Edit it here -- never inside nonlinear.py.
#

HAMILTONIAN = (
    sp.Rational(1, 2) * (px**2 + py**2) * (1 - delta + delta**2)
    - b1 * x * delta
    + sp.Rational(1, 2) * b1**2 * x**2
    + sp.Rational(1, 2) * b2 * (x**2 - y**2)
    + sp.Rational(1, 3) * b3 * (x**3 - 3 * x * y**2)
    + sp.Rational(1, 4) * b4 * (x**4 - 6 * x**2 * y**2 + y**4)
)


# =============================================================================
# 4. POLYNOMIAL SPACE
# =============================================================================
# m = maximum transverse degree
# d = maximum delta degree

ORDER = 8
DELTA_ORDER = 1
N_PLANES = 2


# =============================================================================
# 5. NORMALIZATION / PHYSICAL BOX
# =============================================================================
# Entries follow VARIABLES exactly:
#     [delta, x, y, px, py]

A_BOX = np.array(
    [0.05e-2, 4.5e-3, 3.0e-3, 0.8e-3, 0.6e-3],
    dtype=float,
)


# =============================================================================
# 6. NUMERICAL SETTINGS
# =============================================================================

LEAST_SQUARES_TOL = 1.0e-14
THIN_MULTIPOLE_LENGTH = 1.0e-8
CACHE_REPEATED_MAGNET_MAPS = True
CHECK_ELEMENT_UPPER_RIGHT = False


# =============================================================================
# 7. PLOT SETTINGS
# =============================================================================

PLOT_PLANE = "both"
PLOT_LEVELS = 40
PLOT_GRID_POINTS = 350
PLOT_RMIN = 0.06
PLOT_RMAX = 0.95
PLOT_DELTA = 0.0
PLOT_FOLDER = "nonlinear_plots"

# Physical plotting aperture.
PLOT_X_MAX = 5.0e-3
PLOT_PX_MAX = 1.0e-3
PLOT_Y_MAX = 3.0e-3
PLOT_PY_MAX = 1.0e-3
