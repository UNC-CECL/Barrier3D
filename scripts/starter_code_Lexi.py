from yaml import full_load, dump
from barrier3d import Barrier3d
​
# this is a function that opens a yaml file and changes the variable
def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)
​
​
# yaml variables that will be changed in from this code
def set_barrier3d_yaml(
    datadir="tests/test_params/",
    parameter_file="barrier3d-parameters.yaml",
    elevation_file="InitElevHog.npy",
    dune_file="DuneStart_1000dam.npy",
    storm_file="StormTimeSeries_1000yr.npy",  # File that contains storm data
    shrub_on=0,  # 0 off, 1 on
    tmax=100,  # [y] Duration of simulation
    storm_start=2,  # [y] Year when storm can start occurring (NOTE: if changed, need new storm time series)
    barrier_len=500.0,  # [m] Static length (alongshore) of island segment
    dune_width=20.0,  # [m] Width (cross-shore) of island dune field
    L_shoreface=500.0,  # [m] Initial length of shoreface
    D_shoreface=10.0,  # [m] Height of shoreface
    bay_depth=3.0,  # [m] Depth of bay behind island segment
    MHW=0.46,  # [m NAVD88] Elevation of Mean High Water (NOTE: if changed, need new storm time series)
    dune_param_start=True,  # Dune height will come from external file
    D_start=0.5,  # [m] Initial height of dune domain above berm elevation
    berm_el=1.9,
    growth_param_start=True,  # Dune growth parameter will come from external file
    rmin=0.35,  # Minimum growth rate for logistic dune growth
    rmax=0.85,  # Maximum growth rate for logistic dune growth
    Hd_diffu=0.75,  # [m] Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
    D_max_el=3.4,  # [m NAVD88] Maximum elevation of dunes
    C1=8.8,  # [m] Empirical dune erosion parameter
    C2=4.6,  # [m] Empirical dune erosion parameter
    dune_restart=0.075,  # [m] Restart height for dunes lowered to essentially zero
    Rat=0.0,  # [m / y] Rate of shoreline retreat attributed to alongshore transport; (-) = erosion, (+) = accretion
    SLR_constant=True,
    # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
    SLR_const=0.004,  # [m / y] Relative sea-level rise rate
    beta=0.04,  # Beach slope for runup calculations
    # storm_series=[],  # Time series of storms
    nn=0.5,  # Flow routing constant
    mm=2.0,  # Exponent constant for sediment transport
    Rin_r=2.0,  # Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
    Rin_i=0.1,  # Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
    Qs_min=1.0,  # [m^3 / hr] Minimum discharge needed for sediment transport
    max_up_slope=0.25,  # Maximum slope water can flow upward
    threshold_in=0.25,  # [m^3 / hr] Threshold to determine if in inundation regime
    kr=0.000075,  # Sediment flux constant, run-up regime
    ki=0.0000075,  # Sediment flux constant, inundation regime
    Cbb_r=0.7,  # Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime
    Cbb_i=0.85,  # Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime
    Qs_bb_min=1,  # [m^3 / hr] Minimum sediment flux in back-barrier bay (below which sediment won't flux)
    Cx=10.0,  # Multiplier with the average slope of the interior for constant "C" in inundation transport rule
    OWss_i=2,  # Overwash substep
    OWss_r=1,  # Overwash substep
    k_sf=5000.0,  # [m^3 / m / y] Shoreface flux rate constant
    s_sf_eq=0.02,  # Equilibrium shoreface slope
    seed_min=100.0,  # [1 / yr] Seeds produced per shrub per year (fecundity)
    seed_max=1000.0,  # [1 / yr] Seeds produced per shrub per year (fecundity)
    disp_mu=-0.721891,  # For lognormal probability distribution of seed dispersal distance
    disp_sigma=1.5,  # For lognormal probability distribution of seed dispersal distance
    D_shrub=2.75,  # [m] Minimum elevation of fronting dune for shrub growth
    germ_rate=0.6,  # Germination rate
    time_fruit=5.0,  # [yr] Age shrubs need to be before they start fruiting
    female=0.5,  # Percentage of shrubs that are female
    shrub_el_min=1.2,  # [m NAVD88] Elevation range for shrub growth, minimum bound
    shrub_el_max=2.3,  # [m NAVD88] Elevation range for shrub growth, maximum bound
    spray_dist=170,  # [m] Distance from ocean shoreline that shrubs can establish
    burial_limit=0.75,  # [m] Maximum percentage of height that a shrub can be buried up to before dying
    uproot_limit=-0.2,  # [m] Shrubs eroded beyond this limit killed
    saline_limit=5,  # [m^3 / hr] Dishcharge limit to determine shrub mortality via saline flooding
    Qshrub_max=0.15,  # Maximum percentage of overwash reduction through a shrub cell with full percent cover
    max_shrub_height=5.3,  # [m] Maximum shrub height
    shoreface_toe=0,  # [m] Start location of shoreface toe
    seeded_RNG=True,
):
​
    """
​
    Args:
        elevation_file: File that contains initial elevations [m MHW]
        dune_file: File that contains initial dune height values [m]
        storm_file:
        shrub_on:
        tmax: float,
        storm_start:
        barrier_len:
        dune_width:
        L_shoreface:
        D_shoreface:
        bay_depth:
        MHW:
        dune_param_start:
        D_start:
        berm_el:
        growth_param_start:
        rmin:
        rmax:
        Hd_diffu:
        D_max_el:
        C1:
        C2:
        dune_restart:
        Rat:
        SLR_constant:
        SLR_const:
        beta:
        storm_series:
        nn:
        mm:
        Rin_r:
        Rin_i:
        Qs_min:
        max_up_slope:
        threshold_in:
        kr:
        ki:
        Cbb_r:
        Cbb_i:
        Qs_bb_min:
        Cx:
        OWss_i:
        OWss_r:
        k_sf:
        s_sf_eq:
        seed_min:
        seed_max:
        disp_mu:
        disp_sigma:
        D_shrub:
        germ_rate:
        time_fruit:
        female:
        shrub_el_min:
        shrub_el_max:
        spray_dist:
        burial_limit:
        uproot_limit:
        saline_limit:
        Qshrub_max:
        max_shrub_height:
        shoreface_toe:
        seeded_RNG:
        parameter_file:
        datadir:
​
    Returns:
        nothing -- but modifies the yaml Barrier3D input file
    """
​
    fid = datadir + parameter_file
​
    # update yaml file
    set_yaml("Shrub_ON", shrub_on, fid)  # make sure that shrubs are turned off
    set_yaml("TMAX", tmax, fid)
    set_yaml("StormStart", storm_start, fid)
    set_yaml(
        "BarrierLength", barrier_len, fid
    )  # [m] Static length of island segment (comprised of 10x10 cells)
    set_yaml("DuneWidth", dune_width, fid)
    set_yaml("DShoreface", D_shoreface, fid)
    set_yaml("LShoreface", L_shoreface, fid)
    set_yaml("BayDepth", bay_depth, fid)
    set_yaml("MHW", MHW, fid)
    # (NOTE: if MHW is changed, the storm time series needs to be remade)
    set_yaml(
        "DuneParamStart", dune_param_start, fid
    )  # Dune height will come from external file
    set_yaml("DStart", D_start, fid)
    set_yaml("BermEl", berm_el, fid)
    # (NOTE: if BermEl is changed, the MSSM storm list and storm time series needs to be remade)
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
    set_yaml("rmin", rmin, fid)  # Minimum growth rate for logistic dune growth
    set_yaml("rmax", rmax, fid)  # Maximum growth rate for logistic dune growth
​
    # external file names used for initialization
    set_yaml("storm_file", storm_file, fid)
    set_yaml("dune_file", dune_file, fid)
    set_yaml("elevation_file", elevation_file, fid)
​
​
# example run file (changing 1 parameter in a Barrier 3D run)
# I am confused on this
def run_changing_dune_params(datadir, parameter_file, rmin, rmax):
    # note to Lexi: since these were out of order, you needed to specify the variable name they correspond with!
    set_barrier3d_yaml(
        rmin=rmin, rmax=rmax, datadir=datadir, parameter_file=parameter_file
    )
    a_barrier3d_model = Barrier3d.from_yaml(datadir)
​
    # increase time step
    for time_step in range(1, a_barrier3d_model._TMAX):
​
        # note to Lexi: when using the BMI version of Barrier3D, these commands are grouped behind the scenes in bmi.py
        # take a look and make sure this makes sense to you
        a_barrier3d_model.update()  # update the model's main time loop
        a_barrier3d_model.update_dune_domain()  # now update the dune domain and increase time by one year
​
        # Print time step to screen
        print("\r", "Time Step: ", time_step, end="")
​
    return a_barrier3d_model
​
​
# and example call to the above functions
datadir = "tests/test_params/"
parameter_file = "barrier3d-parameters.yaml"
a_barrier3d_model = run_changing_dune_params(
    datadir=datadir,
    parameter_file=parameter_file,
    rmin=0.45,
    rmax=0.85,
)