"""
    Example run script for Version 1.0 of Barrier3D (NoBMI)
"""

import V1_NoBMI.Barrier3D_Functions as B3Dfunc

datadir = "V1_NoBMI/"

# starts running immediately and ends with plots!
execfile(datadir + "Barrier3D.py")

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(DuneDomain, Dmax)

# Plot 2: Elevation Domain For Last Time Step
B3Dfunc.plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS)

# 3: Elevation Domain Frames
B3Dfunc.plot_ElevFrames(TMAX, DomainTS)
