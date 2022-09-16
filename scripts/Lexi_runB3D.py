from scripts import helper_functions as hlp
from barrier3d import Barrier3d
from matplotlib import pyplot as plt

b3d = Barrier3d.from_yaml("tests/test_params/")

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    b3d.InteriorDomain, cmap="terrain"
)
fig1.colorbar(mat)
ax1.set_title("Initial Elevation $(dam)$")
ax1.set_ylabel("barrier width (dam)")
ax1.set_xlabel("barrier length (dam)")
plt.show()


# datadir = "tests/test_params/"
# parameter_file = "barrier3d-parameters.yaml"
#
# # a_barrier3d_model, nn, qow = hlp.change_nn(datadir, parameter_file, nn=0.7)
# # hlp.plot_nn(nn, qow, a_barrier3d_model)
#
# a_barrier3d_model, mm, Kr, Ki, qow = hlp.change_overwash_params(datadir, parameter_file, Kr=7.5e-05, Ki=7.5e-06, mm=2)
# # hlp.plot_mm(qow=qow, mm=mm, a_barrier3d_model=a_barrier3d_model)
# hlp.plot_overwash_params(Kr=Kr, Ki=Ki, qow=qow, a_barrier3d_model=a_barrier3d_model)
#
