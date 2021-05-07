# Reads WIS data, calculates TWL time series, and counts storms

# Data needed: one .onlns WIS file and one .csv tide gauge file

# IRBR 9Jun20


import numpy as np
import csv
import math
import matplotlib.pyplot as plt


################################################################
### Load WIS data 1980-2014
### KA: from Wahl et al., 2016 - need Hs, Tp, and direction theta (they used offshore wave buoy in 28m water depth, 1980-2013)

# pathpre = ('C:\Barrier3D\Parameterization\Storms\WIS-63183')
pathpre = "/Users/KatherineAnardeWheels/Research/BARis/UNC/VCR/SyntheticStorms"
name = "/ST63183_v03.onlns"


filename = pathpre + name

data = np.loadtxt(filename)
DT = data[:, 0]
Hs = data[:, 9]
Tp = data[:, 11]
WavD = data[:, 15]

# Date Time
year = []
month = []
day = []
hour = []

for n in range(len(DT)):
    DTstring = str(DT[n])
    year.append(int(DTstring[0:4]))
    month.append(int(DTstring[4:6]))
    day.append(int(DTstring[6:8]))
    hour.append(int(DTstring[8:10]))


################################################################
### Load Water Level from CSV
### KA: from Wahl et al., 2016 - need one hour time records of water level (they used tidal gauge)
### KA: what are the decadal trends?
filename2 = pathpre + "/Tide-8631044.csv"
with open(filename2, newline="") as csvfile:
    sl = list(csv.reader(csvfile))
sl = [float(i) for i in sl[0]]  # KA: is this in m NAVD88?
sl = sl[1:-1]  # Trim off first and last readings from tide gauge to match waves
# KA: are these in one hour time steps? Assuming so (other option for tidal gauges is 6 min)
SL = np.asarray(sl)


################################################################
### Running means (i.e., Pre-processing - remove non-stationarity)
### KA: Wahl et al., 2016 - 30 day running medians (shifted by 1 h each time step)

N = 24 * 30  # For 30 day running mean

SL_rm = np.convolve(SL, np.ones((N,)) / N, mode="same")

Hs_rm = np.convolve(Hs, np.ones((N,)) / N, mode="same")

Tp_rm = np.convolve(Tp, np.ones((N,)) / N, mode="same")

# KA: need to check if wave direction is non-stationary over decades by plotting here

# REMOVE THE RUNNING MEDIAN

# WHAT IS THE MEDIAN OVER THE LAST THREE YEARS?

# ADD THIS MEDIAN TO THE NEW CORRECTED TIME SERIES SO REPRESENTATIVE OF RECENT CLIMATE

################################################################
### Storm surge (non-tidal residuals) & astronomical tide

# KA: From the corrected time series, use t-tide python to perform a year-by-year tidal analysis

# Get from Matlab t_Tide for:

# NTR - nontidal residual
# AT - tidal amplitude


################################################################
### R2%, TWL, Rlow
beta = 0.04  # Beach slope, Hog Island
L0 = (9.8 * Tp ** 2) / (2 * math.pi)  # Wavelength
R2 = []
Rlow = []
TWL = []  # KA: this is observed + R2 (not corrected + R2)

# KA: is this Stockdon 2006 broken down into components?
for n in range(len(L0)):
    # Setup
    setup = 0.35 * beta * math.sqrt(Hs[n] * L0[n])

    # Incident band swash
    Sin = 0.75 * beta * math.sqrt(Hs[n] * L0[n])

    # Infragravity band swash
    Sig = 0.06 * math.sqrt(Hs[n] * L0[n])

    # Swash
    swash = math.sqrt((Sin ** 2) + (Sig ** 2))

    # R2%
    r2 = 1.1 * (setup + (swash / 2))
    R2.append(r2)

    # TWL & Rlow
    twl = (
        SL[n] + r2
    )  # KA: not clear from Wahl if this is observations and not corrected....try both?
    TWL.append(twl)
    Rlow.append(twl - (swash / 2))

TWL = np.asarray(TWL)
Rlow = np.asarray(Rlow)
R2 = np.asarray(R2)


################################################################
### Plot
plt.figure()
fig = plt.gcf()
fig.set_size_inches(14, 18)
plt.rcParams.update({"font.size": 14})

x = []  # TEMP: this is not exact!
for t in range(len(year)):
    yr = 1980 + t / (365 * 24)
    x.append(yr)

# SL
plt.subplot(6, 1, 1)
plt.plot(x, SL, color="silver")
plt.plot(x, SL_rm, color="teal")
plt.ylabel("Sea level [mNAVD88]")

# Hs
plt.subplot(6, 1, 2)
plt.plot(x, Hs, color="silver")
plt.plot(x, Hs_rm, color="teal")
plt.ylabel("Hs [m]")

# Tp
plt.subplot(6, 1, 3)
plt.plot(x, Tp, color="silver")
plt.plot(x, Tp_rm, color="teal")
plt.ylabel("Tp [m]")

# WavD
plt.subplot(6, 1, 4)
plt.plot(x, WavD, color="dimgray")
plt.ylabel("Wave Direction [degree]")

# R2
plt.subplot(6, 1, 5)
plt.plot(x, R2, color="dimgray")
plt.ylabel("R2% [m]")

# TWL
plt.subplot(6, 1, 6)
plt.plot(x, TWL, color="dimgray")
plt.ylabel("TWL [mNAVD88]")

plt.show()


################################################################
### Find Storms

# for each year, find when the TWL exceeds an erosion threshold
# From Wahl: they used the 5th percentile of dune toe heights (free parameter for dune erosion)
BermEl = 1.7

# find the annual average TWL from all threshold exceedances from a given year
# calculate annual averages of MSL (here the 30 day running mean), tidal amplitude, residual,
# and R2% during the TWL exceedances

Hs_over_yearly = []
for y in range(35):
    start = 365 * 24 * y
    stop = 365 * 24 * (y + 1)
    hh = Hs[start:stop]  # KA: update to include tidal and nontidal residual and MWL, R2
    tt = TWL[start:stop]
    Hs_over = hh[tt > BermEl]
    Hs_over_yearly.append(np.mean(Hs_over))

# KA: need to confirm that Hs really is driving the large TWL (likely)
# KA: also check that we only consider events where the Hs thresholds were exceeded and the
# simultaneous surge was positive (see comments below)
Hs_min = min(Hs_over_yearly)
Hs_threshold = (
    math.floor(Hs_min / 0.05) * 0.05
)  # Threshold Hs needed to qualify as storm event, rounded to nearst 0.05 m

# KA: separate the seasons (June through November)?

# Plot yearly means
plt.figure()
x = np.linspace(1980, 2014, 35)
plt.plot(x, Hs_over_yearly)
fig = plt.gcf()
fig.set_size_inches(14, 5)
plt.xlabel("Year")
plt.ylabel("Hs [m]")
plt.hlines(Hs_threshold, x[0], x[-1], colors="black", linestyles="dashed")
plt.show()

# Find storms
Storms = np.zeros(
    [0, 12]
)  # Storm start, storm stop, datenum start, datenum stop, Hs, dur, TWL, SL, Tp, R2, Rlow, year

t = 0
while t < len(hour):
    if (
        Hs[t] >= Hs_threshold
    ):  # Future improvement: discard storms where simultaneous surge is negative (but need to separate surge from tide first)
        stormStart = t
        height = Hs[t]
        dur = 1
        t += 1
        while (
            len([x for x in Hs[t : t + 25] if x > Hs_threshold]) > 0
        ):  # If Hs drops below Hs_threshold for only 24 hrs or less, exceedence is assumed part of same weather system (Wahl et al., 2016; Li et al., 2014)
            if Hs[t] > Hs_threshold:
                dur += 1
                t += 1
            else:
                t += 1
        if (
            dur > 8
        ):  # KA: this seems arbitrary, isn't this hours? need to check and see if he did this right (number of storms seem low)
            stormStop = t
            datenumStart = DT[stormStart]
            datenumStop = DT[stormStop]
            storm = [
                stormStart,
                stormStop,
                datenumStart,
                datenumStop,
                max(Hs[stormStart : stormStop + 1]),
                dur,
                max(TWL[stormStart : stormStop + 1]),
                max(SL[stormStart : stormStop + 1]),
                max(Tp[stormStart : stormStop + 1]),
                max(R2[stormStart : stormStop + 1]),
                max(Rlow[stormStart : stormStop + 1]),
                year[stormStart],
            ]
            Storms = np.vstack([Storms, storm])
        t += 1
    else:
        t += 1

# Plot storm TWL histogram
plt.figure()
Bin = np.linspace(0.2, 4.7, 46)
stormtwl = Storms[:, 6]
plt.hist(stormtwl, bins=Bin)
plt.xlabel("Storm TWL [mNAVD88]")
plt.title("1980 - 2014")

print("Max TWL:          ", max(Storms[:, 6]))
print("Max Rexcess:      ", (max(Storms[:, 6]) - BermEl))

print("% Rhigh > BermEl: ", ((Storms[:, 6]) > BermEl).sum() / len(Storms) * 100)
print("% Rlow  > BermEl: ", ((Storms[:, 10]) > BermEl).sum() / len(Storms) * 100)


################################################################
