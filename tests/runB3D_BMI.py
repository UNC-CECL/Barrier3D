# Example run script for the BMI version of Barrier3D
#
# Written by K. Anarde
#
# NOTES to users:
#       - if using Barrier3D for the first time, remember to $ pip install -e .

# import 1) the new class, Barrier3D, which is the model
#        2) load_input: readers for input files, which include the parameter file as well as data files for elev.,
#        storms, etc. Uses the config class to define imported input parameters. NOTE: this is not needed if you use the
#        "from_path" function

from barrier3d import Barrier3dBmi
from tests import Barrier3D_Plotting_Functions as B3Dfunc
import time

# create an instance of the new BMI class, which is the model
barrier3d = Barrier3dBmi()

# specify data directory with initial conditions
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/Barrier3D/tests/test_params/barrier3d-parameters.yaml"
barrier3d.initialize(datadir)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d._model._TMAX):
    barrier3d.update()  # update the model by a time step

    # Print time step to screen
    print("\r", 'Time Step: ', time_step, end="")

SimDuration = time.time() - Time
print()
print('Elapsed Time: ', SimDuration, 'sec')  # Print elapsed time of simulation

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(barrier3d._DuneDomain, barrier3d._Dmax)
# Plot 2: Elevation Domain For Last Time Step
B3Dfunc.plot_ElevTMAX(barrier3d._TMAX, barrier3d._time_index, barrier3d._DuneDomain, barrier3d._DomainTS)
# 3: Elevation Domain Frames
B3Dfunc.plot_ElevFrames(barrier3d._TMAX, barrier3d._DomainTS)
# 4: Animation Frames of Barrier and Dune Elevation
B3Dfunc.plot_ElevAnimation(
    barrier3d._InteriorWidth_AvgTS,
    barrier3d._ShorelineChange,
    barrier3d._DomainTS,
    barrier3d._DuneDomain,
    barrier3d._SL,
    barrier3d._x_s_TS,
    barrier3d._Shrub_ON,
    barrier3d._PercentCoverTS,
    barrier3d._TMAX,
    barrier3d._DeadPercentCoverTS
)
# 5: Cross-Shore Transect Every 100 m Alongshore For Last Time Step
B3Dfunc.plot_XShoreTransects(
    barrier3d._InteriorDomain,
    barrier3d._DuneDomain,
    barrier3d._SL,
    barrier3d._TMAX
)
# 6: Shoreline Positions Over Time
B3Dfunc.plot_ShorelinePositions(barrier3d._x_s_TS, barrier3d._x_b_TS)
# 7: Shoreline Change Rate Over Time
B3Dfunc.plot_ShorelineChangeRate(barrier3d._x_s_TS)
# 8: Run-up vs Inundation count
B3Dfunc.plot_RuInCount(barrier3d._RunUpCount, barrier3d._InundationCount)
# 9: Shoreface LTA14 transects over time
B3Dfunc.plot_LTATransects(barrier3d._SL, barrier3d._TMAX, barrier3d._x_b_TS, barrier3d._x_t_TS, barrier3d._x_s_TS)
# 10: Average Island Elevation Over Time
B3Dfunc.plot_AvgIslandElev(barrier3d._h_b_TS)
# 11: Shoreface Slope Over Time
B3Dfunc.plot_ShorefaceSlope(barrier3d._s_sf_TS)
# 12: Average Interior Width Over Time
B3Dfunc.plot_AvgInteriorWidth(barrier3d._InteriorWidth_AvgTS)
# 13: Shoreface Overwash Flux Over Time
B3Dfunc.plot_OverwashFlux(barrier3d._QowTS)
# 14: Width, Berm Elevation, SF Slope, Shoreline Change, and Overwash Flux Over Time (all in one)
B3Dfunc.plot_StatsSummary(
    barrier3d._s_sf_TS,
    barrier3d._x_s_TS,
    barrier3d._TMAX,
    barrier3d._InteriorWidth_AvgTS,
    barrier3d._QowTS,
    barrier3d._QsfTS,
    barrier3d._Hd_AverageTS
)
# 15: 3D Plot of Island Domain For Last Time Step
B3Dfunc.plot_3DElevTMAX(barrier3d._TMAX, barrier3d._time_index, barrier3d._SL, barrier3d._DuneDomain, barrier3d._DomainTS)
# 16: 3D Animation Frames of Island Elevation and Shrubs (no translation)
B3Dfunc.plot_3DElevFrames(barrier3d._DomainTS, barrier3d._SL, barrier3d._TMAX, barrier3d._DuneDomain)
# 17: 3D Animation Frames of Island Evolution (with barrier translation)
B3Dfunc.plot_3DElevAnimation(
    barrier3d._DomainTS,
    barrier3d._SL,
    barrier3d._TMAX,
    barrier3d._DuneDomain,
    barrier3d._DomainWidth,
    barrier3d._x_s_TS,
    barrier3d._ShorelineChange
)
