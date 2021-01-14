# Example run script for version 1 of Barrier3D (no BMI)
#
# Written by K. Anarde
#
import os
from V1_NoBMI import Barrier3D_Functions as B3Dfunc

os.chdir('/Users/KatherineAnardeWheels/PycharmProjects/Barrier3D/V1_NoBMI')

# (starts running immediately)
execfile('Barrier3D.py')

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(DuneDomain, Dmax)

# Plot 2: Elevation Domain For Last Time Step
B3Dfunc.plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS)

# 3: Elevation Domain Frames
B3Dfunc.plot_ElevFrames(TMAX, DomainTS)
