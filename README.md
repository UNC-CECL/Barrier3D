# Barrier3D
A spatially explicit exploratory model of barrier island evolution in three dimensions.

To see a full listing of the requirements, have a look at the project's
*requirements.txt* file.

If you are a developer of *Barrier3D* you will also want to install
additional dependencies for running *Barrier3D*'s tests to make sure
that things are working as they should. These dependencies are listed
in *requirements-testing.txt* and *requirements-development.txt*.

Installation
------------

To install *Barrier3D*, first create a new environment in
which *Barrier3D* will be installed. This, although not necessary, will
isolate the installation so that there won't be conflicts with your
base *Python* installation. This can be done with *conda* as::

  $ conda create -n barrier3d-env python=3
  $ conda activate barrier3d-env

From Source
+++++++++++

After downloading the *Barrier3D* source code, run the following from
*Barrier3D*'s top-level folder (the one that contains *setup.py*) to
install *Barrier3D* into the current environment::

  $ pip install -e .

Input Files
-----------

Barrier3D Parameter File
++++++++++++++++++++

*Barrier3D* has several input files. The main input file is a yaml-formatted text file that lists
parameter values for the various components. Running the following will
print a sample *Barrier3D* parameter file::

  $ b3d show parameters

.. code:: yaml

TMAX: 150                       # [y] Duration of simulation
StormStart: 2                   # [y] Year when storm can start occurring
BarrierLength: 500.0            # [m] Static length (alongshore) of island segment
DuneWidth: 20.0                 # [m] Width (cross-shore) of island dune field; for illustration purposes only
LShoreface: 500.0               # [m] Length of shoreface
DShoreface: 10.0                # [m] Height of shoreface
BayDepth: 3.0                   # [m] Depth of bay benind island segment
MHW: 0.46                       # [m] Elevation of Mean High Water
Dstart: 0.5                     # [m] Initial height of dune domain above berm elevation
BermEl: 1.9                     # [m] Static elevation of berm; berm elevation + dune height = dune elevation
rmin: 0.35                      # Minimum growth rate for logistic dune growth
rmax: 0.85                      # Maximum growth rate for logistic dune growth
HdDiffu: 0.75                   # [m] Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
Dmaxel: 3.4                     # [m] Maximum elevation of dunes
C1: 8.8                         # [m] Empirical dune erosion parameter
C2: 4.6                         # [m] Empirical dune erosion parameter
DuneRestart: 0.075              # [m] Restart height for dunes lowered to essentially zero
Rat: 0.0                        # [m / y] Rate of shoreline retreat attributed to alongshore transport; (-) = erosion, (+) = accretion
RSLR_Constant: true             # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
RSLR_const: 0.004               # [m / y] Relative sea-level rise rate
beta: 0.04                      # Beach slope for runup calculations
StormSeries: []                 # Time series of storms
nn: 0.5                         # Flow routing constant
mm: 2.0                         # Exponent constant for sediment transport
Rin_r: 2.0                      # Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Rin_i: 0.25                     # Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Qs_min: 1.0                     # [m^3 / hr] Minimum discharge needed for sediment transport
MaxUpSlope: 0.25                # [m / m] Maximum slope water can flow upward
threshold_in: 0.25              # [m^3 / hr] Threshold to determine if in inundation regime
Kr: 7.5e-05                     # Sediment flux constant, run-up regime
Ki: 7.5e-06                     # Sediment flux constant, inundation regime
Cbb_r: 0.5                      # Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime
Cbb_i: 0.8                      # Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime
Qs_bb_min: 1.0                  # [m^3 / hr] Minimum sediment flux in back-barrier bay (below which sediment won't flux)
Cx: 10.0                        # Multiplier with the average slope of the interior for constant C in inundation transport rule
OWss_i: 2                       # Overwash substep
OWss_r: 1                       # Overwash substep
k_sf: 5000.0                    # [m^3 / m / y] Shoreface flux rate constant
s_sf_eq: 0.02                   # Equilibrium shoreface slope
Shrub_ON: false                 # 1 = shrubs on in simulation, 0 = shrubs off
Seedmin: 100.0                  # [1 / yr] Seeds produced per shrub per year (fecundity)
Seedmax: 1000.0                 # [1 / yr] Seeds produced per shrub per year (fecundity)
disp_mu: -0.721891              # For lognormal probability distribution of seed dispersal distance
disp_sigma: 1.5                 # For lognormal probability distribution of seed dispersal distance
Dshrub: 2.0                     # [m] Minimum elevation of fronting dune for shrub growth
GermRate: 0.6                   # Germination rate
TimeFruit: 5.0                  # [yr] Age shrubs need to be before they start fruiting
Female: 0.5                     # Percentage of shrubs that are female
ShrubEl_min: 0.6                # [m] Elevation range for shrub growth, minimum bound
ShrubEl_max: 2.3                # [m] Elevation range for shrub growth, maximum bound
TideAmp: 1.2                    # [m] Tidal amplitude
SprayDist: 170.0                # [m] Distance from ocean shoreline that shrubs can establish
BurialLimit: 0.5                # [m] Shrubs buried beyond this limit killed
UprootLimit: -0.3               # [m] Shrubs eroded beyond this limit killed
SalineLimit: 0.05               # [m^3 / hr] Dishcharge limit to determine shrub mortality via saline flooding
Qshrub_max: 0.15                # Maximum percentage of overwash reduction through a shrub cell with full percent cover
DuneParamStart: true            # Dune height will come from external file
GrowthParamStart: true          # Dune growth parameters will come from external file
MaxShrubHeight: 3.5             # [m] Maximum shrub height
ShorefaceToe: 0.0               # [m] Start location of shoreface toe


Output Files
------------

There are three main sets of output files. These are writen to the 
*output* folder as the model is running.
*  *output/elevation*: elevations of the entire model grid.
*  *output/profile*: elevations along the river profile
*  *output/river*: x, and y coordinates of the river profile

Each of these files is a CSV formatted text file. To create a plot
of one of these output files, use the *plot* subcommand. For example::

  $ rafem plot elevation

will plot the final elevations for the simulation in the current directory.
Use *rafem plot --help* to see further options.

Examples
--------

To run a simulation using the sample input files described above, you first
need to create a set of sample files. This can be done by hand or by running
`rafem setup` to get a default set of parameters that you can then edit.
For example::

  $ mkdir example
  $ cd example
  $ rafem setup

This command has created a new file, *rafem.yaml*, that you can edit for your
simulation.  To run *rafem* using this file::

  $ rafem run

This will have create a new folder, *output*, that contains the output files.
You can look at some of the output with the *plot* subcommand. For example,
the following will create a plot the final elevations::

  $ rafem plot elevation

Use the *--help* option to get help about other command line options.
