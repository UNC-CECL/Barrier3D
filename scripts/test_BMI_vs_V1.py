"""
    Test output from the BMI version versus Ian's original code. Both models source the input storm, dune, and elevation
    files from: '/Users/KatherineAnardeWheels/PycharmProjects/Barrier3d/tests/test_params/'
"""

import os
import time

import matplotlib.pyplot as plt

from barrier3d import Barrier3dBmi
from scripts import Barrier3D_Plotting_Functions as B3Dfunc

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
    print("\r", "Time Step: ", time_step, end="")

SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

# This is currently the only plotting function that works with the BMI version
# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(barrier3d._model._DuneDomain, barrier3d._model._Dmax)

# ----------------------------- #

os.chdir("/Users/KatherineAnardeWheels/PycharmProjects/Barrier3D/V1_NoBMI/")

# (starts running immediately)
execfile("Barrier3D.py")

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(DuneDomain, Dmax)

# I also want to check that shoreline change is the same between the two and that shoreface slope starts in equilibrium
plt.figure()
plt.plot(x_s_TS, "b")
plt.plot(barrier3d._model.x_s_TS, "g")

plt.figure()
plt.plot(s_sf_TS, "b")
plt.plot(barrier3d._model._s_sf_TS, "g")
