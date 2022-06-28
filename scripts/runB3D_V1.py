"""
    Example run script for Version 1.0 of Barrier3D (NoBMI). Note, you must be in the top Barrier3D working directory.
"""
import version1_local_copy.Barrier3D_Functions as B3Dfunc

datadir = "version1_local_copy/"

TMAX = []
t = []
DuneDomain = []
DomainTS = []
Dmax = []

# starts running immediately and ends with plots!
exec(open(datadir + "Barrier3D.py").read())

# Plot 1: Dune Height Over Time (input in decameter)
B3Dfunc.plot_DuneHeight(DuneDomain, Dmax)

