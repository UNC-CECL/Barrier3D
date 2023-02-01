"""
    Example run script for the BMI version of Barrier3D (Version 2.0)

    NOTE to users:
        - if using Barrier3D for the first time, remember to $ pip install -e .
"""

from barrier3d import Barrier3dBmi
from barrier3d.tools import plot as B3Dfunc
import time

# create an instance of the new BMI class, which is the model
barrier3d = Barrier3dBmi()

# specify data directory with initial conditions
datadir = "tests/test_params/"
input_file = "barrier3d-default-parameters.yaml"
barrier3d.initialize(datadir + input_file)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d._model._TMAX):
    barrier3d.update()  # update the model by a time step

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_dune_height(barrier3d._model._DuneDomain, barrier3d._model._Dmax)
B3Dfunc.plot_shrub_percent_cover_tmax(
    barrier3d._model._PercentCoverTS,
    barrier3d._model._TMAX,
    barrier3d._model._DeadPercentCoverTS,
)
B3Dfunc.plot_ElevTMAX(
    barrier3d._model._TMAX,
    barrier3d._model._DuneDomain,
    barrier3d._model._DomainTS,
    barrier3d._model._BermEl,
    barrier3d._model._Shrub_ON,
    barrier3d._model._PercentCoverTS,
    barrier3d._model._DeadPercentCoverTS,
    barrier3d._model._DuneWidth,
)
