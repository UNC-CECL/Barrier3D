#
# K. Anarde's first attempt at the BMI!
#
# NOTE to Ian: to get the command-line interface to work you need to % make install
# NOTE to Eric:  RuntimeWarning: invalid value encountered in double_scalars using command line interface
# Q for Eric: in the command-line, where does the information come from for the parameter yaml file? config?

# import 1) the new class, Barrier3D, which is the model
#        2) load_input: readers for input files, which include the parameter file as well as data files for elev.,
#        storms, etc. Uses the config class to define imported input parameters. NOTE: this is not needed if you use the
#        "from_path" function
from barrier3d import Barrier3d, load_inputs
import Barrier3D_Functions as B3Dfunc
import time

# specify data directory with initial conditions and use load_inputs function: processes input parameters for a barrier3d simulation
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/Barrier3D/tests/test_params"
fmt = "yaml"  # alternatively "py", NOTE to Eric: I can't get "py" to work
params = load_inputs(datadir, prefix="barrier3d", fmt=fmt)
# NOTE to Eric: the above throws and error: FutureWarning: Support for multi-dimensional indexing (e.g. `obj[:, None]`)
# # is deprecated and will be removed in a future version.  Convert to a numpy array before indexing instead.

# initialize model, which currently has five properties: QowTS, QsfTS, x_b_TS, x_s_TS, x_t_TS - will add more later
barrier3d = Barrier3d(**params)

# or I can skip the params step with the "from_path" function...cool!
#barrier3d = Barrier3d.from_path(datadir, fmt=fmt)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d._TMAX):
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
