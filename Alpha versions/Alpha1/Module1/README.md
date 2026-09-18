Linear Lattice Model

This project separates the linear accelerator model into three files:

linear_lattice.py
lattice_config.py
running.py

The separation is intentional:

lattice_config.py
    = machine definition and editable global configuration

linear_lattice.py
    = reusable linear-optics functions

running.py
    = executable test / plotting / reporting program

The main rule is:

Edit the accelerator in lattice_config.py, run and test it with running.py, and do not place machine-specific configuration inside linear_lattice.py.

1. Project Structure

linear_lattice.py

This is the general linear-optics library.

It contains:

magnet construction helpers,

linear transfer matrices,

lattice construction,

periodic Twiss and dispersion calculations,

tunes,

natural chromaticities,

natural horizontal emittance,

first-order chromatic correction,

efficient linear updates,

consistency checks,

plotting and summary helpers.

It intentionally does not contain:

a specific accelerator definition,

editable physical parameters,

lattice ordering,

optimization-variable lists,

global model settings,

an executable main().

The same linear_lattice.py can therefore be used with a different machine by supplying a different configuration file.

lattice_config.py

This is the machine-specific configuration.

Anything that describes the current accelerator or may need to be changed by the user belongs here.

The file contains the global configuration for:

PARAMETERS
VARY

magnet definitions
cell definition
full-ring definition

LINEAR_VARIABLES
CHROMATIC_VARIABLES
NONLINEAR_VARIABLES

PARAMETER_MAP
CORRECTION_PARAMETER_MAP

ENERGY_PARAMETER
REPETITIONS
STEP

CHROMATIC_FAMILY1
CHROMATIC_FAMILY2
TARGET_CHROM_X
TARGET_CHROM_Y

These are module-level global variables.

From another file, use:

import lattice_config as cfg

print(cfg.PARAMETERS)
print(cfg.VARY)
print(cfg.STEP)

This is the preferred way to access the global lattice configuration.

running.py

This is the executable linear-lattice program.

Run it with:

python running.py

It:

reads the accelerator definition from lattice_config.py,

builds the magnets and full ring,

calculates the linear optics,

optionally performs chromatic correction,

optionally runs consistency checks,

optionally tests update_linear(),

always generates a complete text report,

optionally prints the report in the terminal,

optionally generates the linear-optics plot.

running.py is intended to let the user test essentially everything in the linear model without editing linear_lattice.py.

2. Requirements

Install:

pip install numpy matplotlib

3. First Run

Keep the three files in the same directory:

project/
│
├── linear_lattice.py
├── lattice_config.py
└── running.py

Then run:

python running.py

The linear model will be prepared using the current configuration.

A text report is generated automatically:

linear_lattice_report.txt

The report is saved in the current working directory, meaning the directory from which:

python running.py

is executed.

At the end of the run, the exact location is printed:

Report saved to: ...

4. Global Variables vs Runtime Variables

The project deliberately separates global configuration from runtime data.

Global accelerator configuration

Machine-dependent globals belong in:

lattice_config.py

For example:

PARAMETERS = {...}

VARY = [...]

N_CELLS = 20
STEP = 0.01

CHROMATIC_FAMILY1 = "SF1"
CHROMATIC_FAMILY2 = "SD1"

Other modules should access them through:

import lattice_config as cfg

cfg.PARAMETERS
cfg.VARY
cfg.N_CELLS
cfg.STEP

Global run choices

running.py also has a small set of global variables, but these control only what the current run should do.

For example:

RUN_CHECKS = True
RUN_PLOTS = True
RUN_UPDATE_TEST = False

PRINT_PARAMETERS = False
PRINT_MAGNETS = False
PRINT_REPORT = True

CORRECT_CHROMATICITY = True

REPORT_FILE = "linear_lattice_report.txt"

These are run settings, not accelerator parameters.

Local runtime variables

Objects constructed inside:

main()

are local runtime variables.

Typical examples are:

magnets
lattice
data
correction
parameters
checks

They should not be imported from running.py.

If nonlinear.py, optimization.py, or another module needs a prepared lattice, it should prepare the model using linear_lattice.py and lattice_config.py.

5. Editing the Accelerator

Most normal changes should be made only in:

lattice_config.py

The parameter dictionary is:

PARAMETERS = {
    "energy": 3.0,
    "LSD": 0.1,
    "F": 0.8,

    "X1": ...,
    "X2": ...,
    "X3": ...,

    "kse1": ...,
    "ks1": ...,
    "ks2": ...,

    "ko1": ...,
    "ko2": ...,
    "ko3": ...,
}

To change the initial machine, edit these values directly.

For example:

PARAMETERS["X1"] = 3.8
PARAMETERS["ko1"] = 5.0

or simply change their values in the dictionary definition.

6. Magnet Definition

Magnets are defined inside:

lattice_config.py

through:

def define_magnets(parameters):
    ...

using the generic constructor from linear_lattice.py:

lin.magnet(...)

The internal magnet representation is:

[
    name,
    type,
    length,
    angle,
    K,
    S,
    O,
    M,
    M5,
]

where:

Field

Meaning

NAME

magnet/family name

TYPE

magnet type

LENGTH

physical length

ANGLE

total bending angle

K

quadrupole strength

S

sextupole strength

O

thin nonlinear multipole strength

M

cached 4×4 transverse map

M5

cached 5×5 dispersion map

Supported magnet types include:

drift
quadrupole
sextupole
bending
multipole

Repeated occurrences of the same family in the ring share the same magnet object.

7. Reading Magnet Data

The list positions are internal to linear_lattice.py.

Instead of using numeric indices in other modules, use:

lin.magnet_field(element, "NAME")
lin.magnet_field(element, "TYPE")
lin.magnet_field(element, "LENGTH")
lin.magnet_field(element, "K")
lin.magnet_field(element, "S")
lin.magnet_field(element, "O")

For example:

qf1 = lin.get_magnet(lattice, "QF1")

print(lin.magnet_field(qf1, "K"))
print(lin.magnet_field(qf1, "LENGTH"))

This keeps the internal list representation controlled by linear_lattice.py.

8. Lattice Definition

The ring ordering is defined in:

lattice_config.py

through objects such as:

DBA
CELA
CELL_NAMES
N_CELLS
RING_NAMES

Conceptually:

magnet definitions
        ↓
CELL_NAMES
        ↓
repeat N_CELLS times
        ↓
RING_NAMES
        ↓
full lattice

Changing the machine layout should therefore be done in lattice_config.py, not in linear_lattice.py.

9. Main Calculation Workflow

The general calculation is:

cfg.PARAMETERS
      ↓
cfg.define_magnets()
      ↓
build full ring from cfg.RING_NAMES
      ↓
compute magnet linear matrices
      ↓
periodic Twiss + dispersion
      ↓
propagate linear functions
      ↓
tunes + natural chromaticities
      ↓
radiation integrals + emittance
      ↓
optional chromatic correction
      ↓
prepared linear model

running.py performs this complete workflow automatically.

10. Preparing the Lattice From Another Module

For nonlinear.py, optimization.py, or another program, import:

import linear_lattice as lin
import lattice_config as cfg

Then prepare the lattice with:

magnets, lattice, data, correction, parameters = lin.prepare_lattice(
    parameters=cfg.PARAMETERS,
    ring_names=cfg.RING_NAMES,
    magnet_builder=cfg.define_magnets,
    energy_parameter=cfg.ENERGY_PARAMETER,
    correction_parameter_map=cfg.CORRECTION_PARAMETER_MAP,
    correct_chromatic=True,
    family1=cfg.CHROMATIC_FAMILY1,
    family2=cfg.CHROMATIC_FAMILY2,
    target_chrom_x=cfg.TARGET_CHROM_X,
    target_chrom_y=cfg.TARGET_CHROM_Y,
    repetitions=cfg.REPETITIONS,
    step=cfg.STEP,
)

This is the standard programmatic interface.

linear_lattice.py does not import the machine configuration itself. The caller explicitly supplies the configuration it wants to use.

11. Main Returned Objects

prepare_lattice() returns:

magnets, lattice, data, correction, parameters

with:

magnets
    unique magnet definitions

lattice
    ordered full-ring lattice

data
    calculated linear-optics data

correction
    chromatic-correction result, or None

parameters
    parameter dictionary actually used

The returned parameters dictionary is a copy of the input configuration and includes corrected chromatic-family strengths when correction is enabled.

12. Reading Linear Data

The internal layout of data is controlled by linear_lattice.py.

Use:

lin.linear_data(data, "CS0")
lin.linear_data(data, "DISP0")

lin.linear_data(data, "LATTICE_M4")
lin.linear_data(data, "LATTICE_M5")

lin.linear_data(data, "TUNE_X")
lin.linear_data(data, "TUNE_Y")

lin.linear_data(data, "CHROM_X")
lin.linear_data(data, "CHROM_Y")

lin.linear_data(data, "EMITTANCE")
lin.linear_data(data, "CIRCUMFERENCE")

lin.linear_data(data, "S_VALUES")
lin.linear_data(data, "CS_VALUES")
lin.linear_data(data, "DISP_VALUES")

For example:

print(lin.linear_data(data, "TUNE_X"))
print(lin.linear_data(data, "TUNE_Y"))
print(lin.linear_data(data, "EMITTANCE"))

This avoids depending on numeric positions in the result list.

13. Parameter Groups

lattice_config.py classifies editable parameters into three groups.

Linear variables

LINEAR_VARIABLES

These change the linear lattice.

Examples include:

energy
quadrupole strengths
bend scaling
magnet lengths
geometry parameters

A change in one of these variables requires recomputing the linear optics.

Chromatic variables

CHROMATIC_VARIABLES

These change sextupole strengths.

In the current on-momentum linear model, sextupoles have the same linear map as drifts.

Therefore these changes normally reuse the existing Twiss and dispersion but require a new chromatic correction.

Nonlinear variables

NONLINEAR_VARIABLES

These currently contain the thin nonlinear multipole strengths:

ko1
ko2
ko3

They have no linear effect.

Therefore a nonlinear-only change can reuse both the existing linear optics and the existing first-order chromatic correction.

14. Recalculation Rules

The intended update behavior is:

Changed quantity

Linear optics

Chromatic correction

linear/geometric variable

recompute

recompute

sextupole/chromatic variable

reuse

recompute

thin nonlinear multipole

reuse

reuse

In short:

linear change
    → recompute linear optics
    → recompute chromatic correction

sextupole change
    → reuse linear optics
    → recompute chromatic correction

thin nonlinear change
    → reuse linear optics
    → reuse chromatic correction

This is why the variables are classified explicitly in lattice_config.py.

15. Efficient Updates

During optimization it is usually unnecessary to rebuild everything after every parameter change.

The generic update function is:

lin.update_linear(...)

It receives:

current lattice
current linear data
new parameter values
list of edited variables
configuration describing parameter dependencies

and decides what must be recalculated from the variable groups in lattice_config.py.

A typical pattern is:

new_parameters = parameters.copy()
new_parameters["ko1"] = 5.0

lattice, data, correction, parameters = lin.update_linear(
    lattice=lattice,
    data=data,
    parameters=new_parameters,
    edited_variables=["ko1"],
    correct_chromatic=True,
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

This is the path intended for later use by the optimization code.

16. Chromaticity Correction

The correction families are configured globally in:

lattice_config.py

using:

CHROMATIC_FAMILY1 = "SF1"
CHROMATIC_FAMILY2 = "SD1"

TARGET_CHROM_X = 0.0
TARGET_CHROM_Y = 0.0

linear_lattice.py contains only the generic correction algorithm.

The machine-specific choice of families and targets remains in lattice_config.py.

When correction is enabled, the solved family strengths are written into the shared magnet objects and into the returned parameter dictionary.

17. running.py Options

The top of running.py contains the user-editable run controls:

RUN_CHECKS = True
RUN_PLOTS = True
RUN_UPDATE_TEST = False

PRINT_PARAMETERS = False
PRINT_MAGNETS = False
PRINT_REPORT = True

CORRECT_CHROMATICITY = True

REPORT_FILE = "linear_lattice_report.txt"

For example, to run without plots:

RUN_PLOTS = False

To run the consistency tests but avoid printing the entire report in the terminal:

RUN_CHECKS = True
PRINT_REPORT = False

The text report is still generated.

18. Testing update_linear()

running.py can also test the efficient update path.

Enable:

RUN_UPDATE_TEST = True

and specify temporary test values:

UPDATE_TEST_VALUES = {
    "X1": 3.60,
    "ko1": 5.0,
}

These values are only part of the test performed by running.py.

They do not replace the original machine configuration stored in lattice_config.py.

The update-test results are included in the generated text report.

19. Text Report

Every execution of:

python running.py

generates:

linear_lattice_report.txt

unless REPORT_FILE is changed.

The report contains the main information useful for validating the model, including:

general ring information
current parameters
variable groups
periodic Twiss values
periodic dispersion
tunes
natural chromaticities
natural emittance
circumference
chromatic correction
full-ring 4×4 matrix
full-ring 5×5 matrix
consistency checks
sampling information
magnet table
optional update_linear test

The report is intended to be the standard diagnostic output when testing changes to the linear model.

20. Validation

The consistency checks include:

symplecticity
periodic Twiss closure
periodic dispersion closure
Courant-Snyder identities
total ring bending

They can be enabled in running.py with:

RUN_CHECKS = True

For model validation, the principal linear quantities can also be compared with an independent accelerator-optics code.

21. Interface With nonlinear.py

The intended dependency structure is:

lattice_config.py
        │
        │ machine definition
        ▼
linear_lattice.py
        │
        │ prepared lattice + linear data
        ▼
nonlinear.py
        │
        ▼
optimization.py

nonlinear.py should not reconstruct the linear accelerator.

Instead, it should receive or prepare:

lattice
data
parameters

using the same linear_lattice.py + lattice_config.py configuration.

When only nonlinear multipole strengths change, the existing linear data should be reused.

22. Recommended Working Style

For normal lattice editing:

edit lattice_config.py
        ↓
run running.py
        ↓
inspect linear_lattice_report.txt
        ↓
inspect plots/checks

For nonlinear calculations:

import linear_lattice
import lattice_config
        ↓
prepare lattice once
        ↓
pass lattice + data to nonlinear.py

For optimization:

prepare lattice once
        ↓
change selected cfg.VARY parameters
        ↓
update_linear()
        ↓
recompute only what the edited variables require
        ↓
evaluate nonlinear objective

23. Summary

The three files have separate responsibilities:

File

Responsibility

User normally edits it?

linear_lattice.py

general linear-optics functions

No

lattice_config.py

accelerator definition and global model configuration

Yes

running.py

run settings, testing, plots, and text report

Yes, for run options

The central design rule is:

linear_lattice.py
    knows how to calculate

lattice_config.py
    defines what machine to calculate

running.py
    decides what to run and what to display/save

This keeps the computational model reusable while making both the accelerator definition and the testing workflow easy to modify.
