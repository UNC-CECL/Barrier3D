# Makes & saves the list of storms from spreadsheet data for

# ~ Barrier3D ~

# Ian R.B. Reeves

# Version Number: 1
# Updated: 21 Sep 2020


import numpy as np
import xlrd


# Simulated Storms from WIS 63183 & Wachapreague
length = 9986  # Number of storms in list

loc = "C:/Barrier3D/Parameterization/VCR_SyntheticStorms.xlsx"


wb = xlrd.open_workbook(loc)  # Open Excel storms file
s1 = wb.sheet_by_index(1)
StormList = np.zeros([length, 7])
for r in range(length):

    StormList[r, 0] = s1.cell_value(r + 1, 0)  # Hs
    StormList[r, 1] = s1.cell_value(r + 1, 1)  # Dur
    StormList[r, 2] = s1.cell_value(r + 1, 2)  # TWL
    StormList[r, 3] = s1.cell_value(r + 1, 3)  # NTR
    StormList[r, 4] = s1.cell_value(r + 1, 4)  # Tp
    StormList[r, 5] = s1.cell_value(r + 1, 5)  # Tide
    StormList[r, 6] = s1.cell_value(r + 1, 6)  # Rlow


np.save("C:/Barrier3D/Parameterization/StormList.npy", StormList)
