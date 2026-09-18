"""Editable definition of the accelerator lattice used.
"""
import linear_lattice as lin

# 1. PHYSICAL PARAMETERS
PARAMETERS = {
    "energy": 3.0,
    "LSD": 0.1,
    "F": 0.8,

    # Linear / geometric variables
    "X1": 3.633167514008421,
    "X2": -4.258277621861492,
    "X3": -2.690860661253351,
    "X4": 2.754505457254375,
    "X5": -3.336431720192718,
    "X6": -1.492176552197721,
    "X7": 2.943728718874649e-3,
    "X8": 6.010287762639232e-1,
    "X9": 5.370545478298007e-1,

    # Sextupole strengths
    "kse1": -31.67808026513335,
    "kfd2": -51.10107383195414,
    "kfd3": -93.324704323109,
    "ks1": -31.98668968024124,
    "ks2": 41.64059482880812,
    "ksd3": -11.2882623813059,
    "ks1s": 696.35209098583,
    "ks2s": 488.0292055866639,
    "ksf1": 29.726909817,
    "ksd1": -127.01,

    # Thin higher-multipole strengths
    "ko1": 4.72312991538625,
    "ko2": 46.6898786404116,
    "ko3": 23.51531647572548,
}

# 2. OPTIMIZATION VARIABLES
VARY = [
    "kse1", "kfd2", "kfd3", "ks1", "ks2", "ksd3",
    "ks1s", "ks2s", "ko1", "ko2", "ko3",
]

# 3. MAGNET DEFINITIONS
def define_magnets(parameters):
    """Return the unique magnet definitions for the current parameter values."""
    p = parameters
    LSD = p["LSD"]
    F = p["F"]

    return [   #Name/Type/strength
        # Drifts
        lin.magnet("D1", "drift", 2.654400 - LSD),
        lin.magnet("D4", "drift", 0.081240),
        lin.magnet("D11", "drift", 0.063628),
        lin.magnet("D12", "drift", 0.0099526),
        lin.magnet("D5D6", "drift", p["X8"]),
        lin.magnet("D9D10", "drift", p["X9"]),
        # Quadrupoles
        lin.magnet("QF1", "quadrupole", 0.349140, K_value=p["X1"]),
        lin.magnet("QD2", "quadrupole", 0.222950, K_value=p["X2"]),
        lin.magnet("QD3", "quadrupole", 0.194780, K_value=p["X3"]),
        lin.magnet("QF4", "quadrupole", 0.224580, K_value=p["X4"]),
        lin.magnet("QD5", "quadrupole", 0.210950, K_value=p["X5"]),
        lin.magnet("QF7", "quadrupole", 0.020986, K_value=p["X6"]),
        # Sextupoles
        lin.magnet("SE1", "sextupole", LSD, S_value=p["kse1"]),
        lin.magnet("FD2", "sextupole", 0.094502, S_value=p["kfd2"]),
        lin.magnet("FD3", "sextupole", p["X7"], S_value=p["kfd3"]),
        lin.magnet("S1", "sextupole", LSD, S_value=p["ks1"]),
        lin.magnet("S2", "sextupole", LSD, S_value=p["ks2"]),
        lin.magnet("SD3", "sextupole", 0.010176, S_value=p["ksd3"]),
        lin.magnet("S1S", "sextupole", 0.002964, S_value=p["ks1s"]),
        lin.magnet("S2S", "sextupole", 0.172130, S_value=p["ks2s"]),
        lin.magnet("SF1", "sextupole", 0.220440, S_value=p["ksf1"]),
        lin.magnet("SD1", "sextupole", LSD, S_value=p["ksd1"]),
        # Thin higher multipoles. O is the integrated nonlinear strength.
        lin.magnet("O1", "multipole", 0.0, O_value=p["ko1"]),
        lin.magnet("O2", "multipole", 0.0, O_value=p["ko2"]),
        lin.magnet("O3", "multipole", 0.0, O_value=p["ko3"]),
        # Bending / combined-function magnets
        lin.magnet("DQ6", "bending", 0.275390, angle=-0.73179259 * F, K_value=2.692600),
        lin.magnet("A1", "bending", 0.075497, angle=0.0021719 * F),
        lin.magnet("A2", "bending", 0.384040, angle=0.53380 * F),
        lin.magnet("A3", "bending", 0.001995, angle=0.00032534 * F),
        lin.magnet("A4", "bending", 0.913400, angle=2.0382 * F),
        lin.magnet("A5", "bending", 0.152490, angle=0.93133 * F),
        lin.magnet("B1", "bending", 0.400570, angle=0.63294 * F),
        lin.magnet("B2", "bending", 0.563170, angle=1.1254 * F),
        lin.magnet("B3", "bending", 0.362720, angle=1.1741 * F),
        lin.magnet("B4", "bending", 0.285610, angle=1.4465 * F),
        lin.magnet("B5", "bending", 0.240960, angle=0.58358 * F),
        lin.magnet("B1S", "bending", 0.015767, angle=0.080780 * F),
        lin.magnet("B2S", "bending", 0.001644, angle=-0.00041155 * F),
        lin.magnet("B3S", "bending", 0.212550, angle=1.7586 * F),
        lin.magnet("DQ1S", "bending", 0.257080, angle=0.81690 * F, K_value=-5.135300),
        lin.magnet("ABQ1", "bending", 0.215990, angle=-0.60542 * F, K_value=6.191000),
    ]


# 4. FULL-RING
DA1 = ["A1", "A2", "A3", "A4", "A5"]
IDA1 = DA1[::-1]

DBA = [
    "D1", "SE1", "QF1", "FD2", "QD2", "FD3",
    *IDA1,
    "D4", "QD3", "SD1", "O2", "D5D6", "S1", "QF4", "SF1",
    "O1", "QF4", "S2", "D9D10", "O3", "SD1", "QD5", "D11",
    "B1", "B2", "B3", "B4", "B5", "D12", "QF7", "SD3", "DQ6",
]

CELA = [
    "S1S", "ABQ1", "S2S", "DQ1S", "B1S", "B2S", "B3S",
    "B2S", "B1S", "DQ1S", "S2S", "ABQ1", "S1S",
]

CELL_NAMES = DBA + CELA + CELA + CELA + DBA[::-1]
N_CELLS = 1#20
RING_NAMES = CELL_NAMES * N_CELLS

# 5. PARAMETER DEPENDENCIES  (Errors here are important to avoid)

LINEAR_VARIABLES = { #Variables that if changed, then we need to recompute the linear part.
    "energy", "LSD", "F", "X1", "X2", "X3", "X4", "X5", "X6", 
    "X7", "X8", "X9",
} 

CHROMATIC_VARIABLES = { #Variables that if changed, we need to recompute the chromatic sextupoles.
    "kse1", "kfd2", "kfd3", "ks1", "ks2", "ksd3", "ks1s", "ks2s",
    "ksf1", "ksd1",
}

NONLINEAR_VARIABLES = {"ko1", "ko2", "ko3"}  #Variables that if changed, they only changed the non-linear transport matrix.

PARAMETER_MAP = {  #Which magnets do each variable affects. (name,parameter affected)
    "energy": [],
    "X1": [("QF1", "K")],
    "X2": [("QD2", "K")],
    "X3": [("QD3", "K")],
    "X4": [("QF4", "K")],
    "X5": [("QD5", "K")],
    "X6": [("QF7", "K")],
    "X7": [("FD3", "LENGTH")],
    "X8": [("D5D6", "LENGTH")],
    "X9": [("D9D10", "LENGTH")],
    "kse1": [("SE1", "S")],
    "kfd2": [("FD2", "S")],
    "kfd3": [("FD3", "S")],
    "ks1": [("S1", "S")],
    "ks2": [("S2", "S")],
    "ksd3": [("SD3", "S")],
    "ks1s": [("S1S", "S")],
    "ks2s": [("S2S", "S")],
    "ksf1": [("SF1", "S")],
    "ksd1": [("SD1", "S")],
    "ko1": [("O1", "O")],
    "ko2": [("O2", "O")],
    "ko3": [("O3", "O")],
    "LSD": [
        ("D1", "LENGTH"), ("SE1", "LENGTH"), ("S1", "LENGTH"),
        ("S2", "LENGTH"), ("SD1", "LENGTH"),
    ],
    "F": [
        (name, "ANGLE")
        for name in (
            "DQ6", "A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3",
            "B4", "B5", "B1S", "B2S", "B3S", "DQ1S", "ABQ1",
        )
    ],
}

CORRECTION_PARAMETER_MAP = {  #SExtupoles used for chromatic correction
    "SF1": "ksf1",
    "SD1": "ksd1",
}

# 6. LINEAR-MODEL SETTINGS
ENERGY_PARAMETER = "energy"
REPETITIONS = 1  #what is this used for?
STEP = 0.01 #subdivition used for linear computations
CHROMATIC_FAMILY1 = "SF1"  #Chromatic term 1   #Maybe I can make this more general.
CHROMATIC_FAMILY2 = "SD1"  #Chromatic term 2
TARGET_CHROM_X = 0.0
TARGET_CHROM_Y = 0.0
