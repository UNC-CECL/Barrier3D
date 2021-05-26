"""
    Test that the BMI version of Barrier 3D (Version 2.0) is equivalent to Ian's original code in Version 1.0. Note that
    the input files for Version 1.0 are located in `V1_NoBMI/Parameter` and for Version 2.0 in 'tests/test_params/'.
    Here we test the model with SHRUBS ON.

    NOTE to users:
      - if using Barrier3D for the first time, remember to $ pip install -e .
"""

import time
import matplotlib.pyplot as plt
from barrier3d import Barrier3dBmi
from barrier3d.tools import plot as B3Dfunc

# specify data directories with initial conditions
datadir_V1 = "V1_NoBMI/"
datadir_V2 = "tests/test_params/"

# Version 2.0 ------------------------------
# create an instance of the new BMI class, which is the model
barrier3d = Barrier3dBmi()
input_file = "barrier3d-parameters.yaml"
barrier3d.initialize(datadir_V2 + input_file)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d._model._TMAX):
    barrier3d.update()
    print("\r", "Time Step: ", time_step, end="")
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

# Version 1.0 ------------------------------

# starts running immediately and ends with plots!
execfile(datadir_V1 + "Barrier3D.py")

# ----------------------------- #

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_dune_height(barrier3d._model._DuneDomain, barrier3d._model._Dmax)
B3Dfunc.plot_dune_height(DuneDomain, Dmax)

# check that shoreline change is the same between the two and that shoreface slope starts in equilibrium
plt.figure()
plt.plot(x_s_TS, "b")
plt.plot(barrier3d._model.x_s_TS, "g")

plt.figure()
plt.plot(s_sf_TS, "b")
plt.plot(barrier3d._model._s_sf_TS, "g")
