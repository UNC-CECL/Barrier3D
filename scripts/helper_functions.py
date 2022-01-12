from matplotlib import pyplot as plt
from yaml import full_load, dump
from barrier3d import Barrier3d
import numpy as np



# this is a function that opens a yaml file and changes the variable
def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)

    """
    Args: 
        var_name: name in the yaml file
        new_vals: name in the set_barrier3d_yaml function below
        file_name: location of the yaml file

    Returns: nothing, uploads yaml file/links old and new variable names
    """


# yaml variables that will be changed from this code (all of them included)
def set_barrier3d_yaml(
        datadir="tests/test_params/",
        parameter_file="barrier3d-parameters.yaml",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
        storm_file="StormTimeSeries_1000yr.npy",
        shrub_on=0,
        tmax=100,
        storm_start=2,
        barrier_len=500.0,
        dune_width=20.0,
        L_shoreface=500.0,
        D_shoreface=10.0,
        bay_depth=3.0,
        MHW=0.46,
        dune_param_start=True,
        D_start=0.5,
        berm_el=1.9,
        growth_param_start=True,
        rmin=0.35,
        rmax=0.85,
        Hd_diffu=0.75,
        D_max_el=3.4,
        C1=8.8,
        C2=4.6,
        dune_restart=0.075,
        Rat=0.0,
        SLR_constant=True,
        SLR_const=0.004,
        beta=0.04,
        # storm_series=[],
        nn=0.5,
        mm=2.0,
        Rin_r=2.0,
        Rin_i=0.1,
        Qs_min=1.0,
        max_up_slope=0.25,
        threshold_in=0.25,
        kr=0.000075,
        ki=0.0000075,
        Cbb_r=0.7,
        Cbb_i=0.85,
        Qs_bb_min=1,
        Cx=10.0,
        OWss_i=2,
        OWss_r=1,
        k_sf=5000.0,
        s_sf_eq=0.02,
        seed_min=100.0,
        seed_max=1000.0,
        disp_mu=-0.721891,
        disp_sigma=1.5,
        D_shrub=2.75,
        germ_rate=0.6,
        time_fruit=5.0,
        female=0.5,
        shrub_el_min=1.2,
        shrub_el_max=2.3,
        spray_dist=170,
        burial_limit=0.75,
        uproot_limit=-0.2,
        saline_limit=5,
        Qshrub_max=0.15,
        max_shrub_height=5.3,
        shoreface_toe=0,
        seeded_RNG=True,
):
    """

    Args:
        datadir:                Location of the parameter files
        parameter_file:         yaml file that contains the variables
        elevation_file:         File that contains the initial elevation data
        dune_file:              File that contains the dune data
        storm_file:             File that contains storm data
        shrub_on:               Turns shrubs on and off (0 off, 1 on); make sure they are off
        tmax:                   [y] Duration of simulation
        storm_start:            [y] Year when storm can start occurring (NOTE: if changed, need new storm time series)
        barrier_len:            [m] Static length (alongshore) of island segment (comprised of 10x10 cells)
        dune_width:             [m] Width (cross-shore) of island dune field
        L_shoreface:            [m] Initial length of shoreface
        D_shoreface:            [m] Height of shoreface
        bay_depth:              [m] Depth of bay behind island segment
        MHW:                    [m NAVD88] Elevation of Mean High Water (NOTE: if changed, need new storm time series)
        dune_param_start:       Dune height will come from external file (True/False)
        D_start:                [m] Initial height of dune domain above berm elevation
        berm_el:                (NOTE: if BermEl is changed, the MSSM storm list and storm time series needs to be remade)
        growth_param_start:     Dune growth parameter will come from external file (True/False)
        rmin:                   Minimum growth rate for logistic dune growth
        rmax:                   Maximum growth rate for logistic dune growth
        Hd_diffu:               Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
        D_max_el:               [m NAVD88] Maximum elevation of dunes
        C1:                     [m] Empirical dune erosion parameter
        C2:                     [m] Empirical dune erosion parameter
        dune_restart:           [m] Restart height for dunes lowered to essentially zero
        Rat:                    [m / y] Rate of shoreline retreat attributed to alongshore transport; (-) = erosion, (+) = accretion
        SLR_constant:           Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series (True/False)
        SLR_const:              [m / y] Relative sea-level rise rate
        beta:                   Beach slope for runup calculations
        storm_series:           Time series of storms
        nn:                     Flow routing constant
        mm:                     Exponent constant for sediment transport
        Rin_r:                  Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
        Rin_i:                  Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
        Qs_min:                 [m^3 / hr] Minimum discharge needed for sediment transport
        max_up_slope:           Maximum slope water can flow upward
        threshold_in:           [m^3 / hr] Threshold to determine if in inundation regime
        kr:                     Sediment flux constant, run-up regime
        ki:                     Sediment flux constant, inundation regime
        Cbb_r:                  Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime
        Cbb_i:                  Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime
        Qs_bb_min:              [m^3 / hr] Minimum sediment flux in back-barrier bay (below which sediment won't flux)
        Cx:                     Multiplier with the average slope of the interior for constant "C" in inundation transport rule
        OWss_i:                 Overwash substep
        OWss_r:                 Overwash substep
        k_sf:                   [m^3 / m / y] Shoreface flux rate constant
        s_sf_eq:                Equilibrium shoreface slope
        seed_min:               [1 / yr] Seeds produced per shrub per year (fecundity)
        seed_max:               [1 / yr] Seeds produced per shrub per year (fecundity)
        disp_mu:                For lognormal probability distribution of seed dispersal distance
        disp_sigma:             For lognormal probability distribution of seed dispersal distance
        D_shrub:                [m] Minimum elevation of fronting dune for shrub growth
        germ_rate:              Germination rate
        time_fruit:             [yr] Age shrubs need to be before they start fruiting
        female:                 Percentage of shrubs that are female
        shrub_el_min:           [m NAVD88] Elevation range for shrub growth, minimum bound
        shrub_el_max:           [m NAVD88] Elevation range for shrub growth, maximum bound
        spray_dist:             [m] Distance from ocean shoreline that shrubs can establish
        burial_limit:           [m] Maximum percentage of height that a shrub can be buried up to before dying
        uproot_limit:           [m] Shrubs eroded beyond this limit killed
        saline_limit:           [m^3 / hr] Dishcharge limit to determine shrub mortality via saline flooding
        Qshrub_max:             Maximum percentage of overwash reduction through a shrub cell with full percent cover
        max_shrub_height:       [m] Maximum shrub height
        shoreface_toe:          [m] Start location of shoreface toe
        seeded_RNG:             Use seeded random number generator for reproducibility (True/False)

    Returns:
        nothing -- but modifies the yaml Barrier3D input file
    """

    fid = datadir + parameter_file

    # update yaml file
    set_yaml("Shrub_ON", shrub_on, fid)
    set_yaml("TMAX", tmax, fid)
    set_yaml("StormStart", storm_start, fid)
    set_yaml("BarrierLength", barrier_len, fid)
    set_yaml("DuneWidth", dune_width, fid)
    set_yaml("DShoreface", D_shoreface, fid)
    set_yaml("LShoreface", L_shoreface, fid)
    set_yaml("BayDepth", bay_depth, fid)
    set_yaml("MHW", MHW, fid)
    set_yaml("DuneParamStart", dune_param_start, fid)
    set_yaml("DStart", D_start, fid)
    set_yaml("BermEl", berm_el, fid)
    set_yaml(
        "GrowthParamStart", growth_param_start, fid
    )  # Dune growth parameter WILL NOT come from external file
    set_yaml("HdDiffu", Hd_diffu, fid)
    set_yaml("Dmaxel", D_max_el, fid)
    set_yaml("C1", C1, fid)
    set_yaml("C2", C2, fid)
    set_yaml("DuneRestart", dune_restart, fid)
    set_yaml("Rat", Rat, fid)
    set_yaml("RSLR_Constant", SLR_constant, fid)
    set_yaml("RSLR_const", SLR_const, fid)
    set_yaml("beta", beta, fid)
    set_yaml("StormSeries", [], fid)  # lets just always set this to empty
    set_yaml("nn", nn, fid)
    set_yaml("mm", mm, fid)
    set_yaml("Rin_r", Rin_r, fid)
    set_yaml("Rin_i", Rin_i, fid)
    set_yaml("Qs_min", Qs_min, fid)
    set_yaml("MaxUpSlope", max_up_slope, fid)
    set_yaml("threshold_in", threshold_in, fid)
    set_yaml("Kr", kr, fid)
    set_yaml("Ki", ki, fid)
    set_yaml("Cbb_r", Cbb_r, fid)
    set_yaml("Cbb_i", Cbb_i, fid)
    set_yaml("Qs_bb_min", Qs_bb_min, fid)
    set_yaml("Cx", Cx, fid)
    set_yaml("OWss_i", OWss_i, fid)
    set_yaml("OWss_r", OWss_r, fid)
    set_yaml("k_sf", k_sf, fid)
    set_yaml("s_sf_eq", s_sf_eq, fid)
    set_yaml("Seedmin", seed_min, fid)
    set_yaml("Seedmax", seed_max, fid)
    set_yaml("Seedmin", seed_min, fid)
    set_yaml("disp_mu", disp_mu, fid)
    set_yaml("disp_sigma", disp_sigma, fid)
    set_yaml("DShrub", D_shrub, fid)
    set_yaml("GermRate", germ_rate, fid)
    set_yaml("TimeFruit", time_fruit, fid)
    set_yaml("Female", female, fid)
    set_yaml("ShrubEl_min", shrub_el_min, fid)
    set_yaml("ShrubEl_max", shrub_el_max, fid)
    set_yaml("SprayDist", spray_dist, fid)
    set_yaml("BurialLimit", burial_limit, fid)
    set_yaml("UprootLimit", uproot_limit, fid)
    set_yaml("SalineLimit", saline_limit, fid)
    set_yaml("Qshrub_max", Qshrub_max, fid)
    set_yaml("MaxShrubHeight", max_shrub_height, fid)
    set_yaml("ShorefaceToe", shoreface_toe, fid)
    set_yaml("SeededRNG", seeded_RNG, fid)
    set_yaml("rmin", rmin, fid)
    set_yaml("rmax", rmax, fid)

    # external file names used for initialization
    set_yaml("storm_file", storm_file, fid)
    set_yaml("dune_file", dune_file, fid)
    set_yaml("elevation_file", elevation_file, fid)


# example run file (changing 1 parameter in a Barrier 3D run)
def run_changing_dune_params(datadir, parameter_file, rmin, rmax):
    # note to Lexi: since these were out of order, you needed to specify the variable name they correspond with!
    set_barrier3d_yaml(
        rmin=rmin, rmax=rmax, datadir=datadir, parameter_file=parameter_file
    )
    a_barrier3d_model = Barrier3d.from_yaml(datadir)

    # increase time step
    for time_step in range(1, a_barrier3d_model._TMAX):
        a_barrier3d_model.update()  # update the model's main time loop
        a_barrier3d_model.update_dune_domain()  # now update the dune domain and increase time by one year

        # Print time step to screen
        print("\r", "Time Step: ", time_step, end="")

    return a_barrier3d_model

# and example call to the above functions
# datadir = "tests/test_params/"
# parameter_file = "barrier3d-parameters.yaml"
# a_barrier3d_model = run_changing_dune_params(
#     datadir=datadir,
#     parameter_file=parameter_file,
#     rmin=0.45,
#     rmax=0.85,
# )
# nn influences Qi
def change_nn(datadir, parameter_file, nn):
    set_barrier3d_yaml(
        nn=nn, datadir=datadir, parameter_file=parameter_file
    )
    a_barrier3d_model = Barrier3d.from_yaml(datadir)
    # increase time step
    for time_step in range(1, a_barrier3d_model._TMAX):
        a_barrier3d_model.update()  # update the model's main time loop
        a_barrier3d_model.update_dune_domain()  # now update the dune domain and increase time by one year
        # Print time step to screen
        print("\r", "Time Step: ", time_step, end="")

    qow = a_barrier3d_model.QowTS

    return a_barrier3d_model, nn, qow

def plot_nn(nn, qow, a_barrier3d_model):
    x = np.linspace(1,a_barrier3d_model._TMAX, len(qow))
    plt.scatter(x, qow)
    plt.title('nn = ' + str(nn))
    plt.xlabel('time step')
    plt.ylabel('Overwash (m3/m)')

# Kr, Ki, mm influence Qsi
def change_overwash_params(datadir, parameter_file, Kr, Ki, mm):
    set_barrier3d_yaml(
        mm=mm, datadir=datadir, parameter_file=parameter_file, kr=Kr, ki=Ki
    )
    a_barrier3d_model = Barrier3d.from_yaml(datadir)
    # increase time step
    for time_step in range(1, a_barrier3d_model._TMAX):
        a_barrier3d_model.update()  # update the model's main time loop
        a_barrier3d_model.update_dune_domain()  # now update the dune domain and increase time by one year
        # Print time step to screen
        print("\r", "Time Step: ", time_step, end="")

    qow = a_barrier3d_model.QowTS

    return a_barrier3d_model, mm, Kr, Ki, qow

def plot_mm(qow, mm, a_barrier3d_model):
    x = np.linspace(1,a_barrier3d_model._TMAX, len(qow))
    plt.figure(1)
    plt.scatter(x, qow, c='tab:green')
    plt.title('mm = ' + str(mm))
    plt.xlabel('time step')
    plt.ylabel('Overwash (m3/m)')
def plot_overwash_params(Kr, Ki, qow, a_barrier3d_model):
    x = np.linspace(1,a_barrier3d_model._TMAX, len(qow))
    plt.figure(2)
    plt.scatter(x, qow, c='tab:purple')
    plt.title('Kr = ' + str(Kr) + ' Ki = ' + str(Ki))
    plt.xlabel('time step')
    plt.ylabel('Overwash (m3/m)')