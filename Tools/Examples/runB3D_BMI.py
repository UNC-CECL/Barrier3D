# Example run script for the BMI version of Barrier3D
#
# Written by K. Anarde
#
# NOTES to users:
#       - if using Barrier3D for the first time, remember to $ pip install -e .

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

# This is currently the only plotting function that works with the BMI version
# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(barrier3d._model._DuneDomain, barrier3d._model._Dmax)

