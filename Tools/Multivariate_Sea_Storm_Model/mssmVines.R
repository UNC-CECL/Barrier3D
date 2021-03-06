# mssmVines <- function(U) {
# -----------------------------mssmVines----------------------------
# Purpose: This function utilizes the CDVine package to fit a Vine 
# copula model to U, the d-dimensional copula data set in [0,1] space 
#
# Inputs:
#       - U: the uniform marginals     
#
# SEE ALSO: VineCopula Package, MultiHazard toolbox
#
# Record of revisions:
#       Date            Programmer          Description of Change
#       =========================================================
#       08/10/20        K. Anarde           original code
#
# ------------------------------user input------------------------------

# set working directory and load R packages
setwd("/Users/KatherineAnardeWheels/Research/BARis/UNC/VCR/SyntheticStorms") 

#library(VineCopula)
library(CDVine)
library(R.matlab)

# user inputs
U <- read.table("U_mssmVCR.txt")   # N x d data matrix (with uniform margins)
nSimNum = 10000;

# ------------------------------model-------------------------------

# from MultiHazard toolbox:
# Standard trivariate copulas lack flexibility to model joint distributions where heterogeneous dependencies exist between the variable pairs. Pair copula constructions construct multivariate distribution using a cascade of bivariate copulas (some of which are conditional). As the dimensionality of the problem increases the number of mathematically equally valid decompositions quickly becomes large. Bedford and Cooke (2001,2002) introduced the regular vine, a graphical model which helps to organize the possible decompositions. The Canonical (C-) and D- vine are two commonly utilized sub-categories of regular vines, in the trivariate case a vine copula is simultaneously a C- and D-vine. Lets fit a regular vine copula model

# copula model parameters
type  <- "DVine";          # options are "CVine" or "DVine"

# returns pair-copula families composing the C- or D-vine copula ($family), 
# its parameters ($par and $par2) 
Model <- CDVineCopSelect(U, familyset = NA, type = type, selectioncrit="AIC", indeptest=FALSE, level=0.05)

# generate new data Usim from the copula
uSim  <- CDVineSim(nSimNum, family=Model$family, par=Model$par, Model$par2, type=type)

# write outputs to a mat file
writeMat('Usim_mssmVCR-Dvine.mat', uSim = uSim, fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)