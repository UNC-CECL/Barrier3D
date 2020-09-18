# K. Anarde's first attempt at the BMI!

# NOTE to Ian: to get the command-line interface to work you need to % make install
# NOTE to Eric:  RuntimeWarning: invalid value encountered in double_scalars using command line interface
# Q for Eric: in the command-line, where does the information come from for the parameter yaml file? config?

# import 1) the new class, Barrier3D, which is the model
#        2) load_input: readers for input files, which include the parameter file as well as data files for elev.,
#        storms, etc. Uses the config class to define imported input parameters. NOTE: this is not needed if you use the
#        "from_path" function
from barrier3d import Barrier3d, load_inputs

# specify data directory with initial conditions and use load_inputs function: processes input parameters for a barrier3d simulation
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/Barrier3D/tests/test_params"
fmt = "yaml"  # alternatively "py", NOTE to Eric: I can't get "py" to work
params = load_inputs(datadir, prefix="barrier3d", fmt=fmt)
# NOTE to Eric: the above throws and error: FutureWarning: Support for multi-dimensional indexing (e.g. `obj[:, None]`)
# # is deprecated and will be removed in a future version.  Convert to a numpy array before indexing instead.

# initialize model, which currently has five properties: QowTS, QsfTS, x_b_TS, x_s_TS, x_t_TS - will add more later
barrier3d = Barrier3d(**params)

# or I can skip the params step with the "from_path" function...cool!
barrier3d = Barrier3d.from_path(datadir, fmt=fmt)

# increase time step
barrier3d.update()