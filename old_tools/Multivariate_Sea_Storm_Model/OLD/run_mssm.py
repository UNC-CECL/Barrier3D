from oct2py import Oct2Py
import os

os.chdir(
    "/Users/KatherineAnardeWheels/PyCharmProjects/Barrier3D/Tools/Multivariate_Sea_Storm_Model"
)
# octave.addpath(
#     "/Users/KatherineAnardeWheels/PyCharmProjects/Barrier3D/Tools/Multivariate_Sea_Storm_Model"
# )
sCopula = "c-vine"
sWIS_filename = "../ST63183_v03.onlns"
sWaterLevel_filename = "../Tide-8631044-Combined.txt"
fBeta = 0.04
fBermEl = 1.9
nSimStorm = 10000
bPlot = 0

oc = Oct2Py()
[stStorms, stSimStorms] = oc.multivariateSeaStorm(
    sCopula, sWIS_filename, sWaterLevel_filename, fBeta, fBermEl, nSimStorm, bPlot
)
