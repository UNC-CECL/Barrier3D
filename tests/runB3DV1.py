# Now Ian's model (starts running immediately)
execfile('Barrier3D.py')

import Barrier3D_Functions as B3Dfunc

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(DuneDomain, Dmax)

# Plot 2: Elevation Domain For Last Time Step
B3Dfunc.plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS)

# 3: Elevation Domain Frames
B3Dfunc.plot_ElevFrames(TMAX, DomainTS)
