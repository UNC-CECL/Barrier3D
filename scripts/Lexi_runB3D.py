from scripts import helper_functions as hlp
from barrier3d import Barrier3d
from matplotlib import pyplot as plt
import numpy as np
import os
# from cascade.outwasher_reorg_b3d import plot_ElevAnimation


# path = "C:/Users/Lexi/Documents/Research/Outwasher/Output/"
# runID = "b3d_domains"
# newpath = path + runID + "/"
# newpath = "D:/NC State/Outwasher/Output/b3d_domains/"

b3d = Barrier3d.from_yaml("tests/test_params/")

# print("Cx = ", b3d._Cx)
# AvgSlope = b3d._BermEl / 20
# print("B3D avg slope of the interior = ", AvgSlope)

directory = "D:/NC State/Outwasher/Output/b3d_domains/"
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)
# hours = b3d._TMAX
# elevs = np.zeros([b3d._TMAX, b3d._DomainWidth, b3d._BarrierLength])
# for time_step in range(1, b3d._TMAX):
for time_step in range(1, 125):
    b3d.update()  # update the model's main time loop
    b3d.update_dune_domain()  # now update the dune domain and increase time by one year
    int_domain = np.flip(b3d.InteriorDomain)
    os.chdir(directory)
    newpath = "D:/NC State/Outwasher/Output/b3d_domains/Elevation_data/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)
    name = "domain_" + str(time_step)
    np.save(newpath + name, int_domain)
    # AnimateDomain = int_domain
    #
    # # Plot and save
    # plt.rcParams.update({"font.size": 15})
    # elevFig1 = plt.figure(figsize=(15, 7))
    # ax = elevFig1.add_subplot(111)
    # cax = ax.matshow(
    #     AnimateDomain,
    #     # origin="upper",
    #     cmap="terrain",
    #     vmin=-0.3, vmax=0.4
    # )  # , interpolation='gaussian') # analysis:ignore
    # ax.xaxis.set_ticks_position("bottom")
    # elevFig1.colorbar(cax)
    # plt.xlabel("Alongshore Distance (dam)")
    # plt.ylabel("Cross-Shore Distance (dam)")
    # plt.title("Elevation (dam)")
    # plt.tight_layout()
    # timestr = "Time = " + str(time_step) + " hrs"
    # plt.text(1, 1, timestr)
    # name = "elev_" + str(time_step)
    # elevFig1.savefig(name, facecolor='w')  # dpi=200
    # plt.close(elevFig1)
    # # Print time step to screen
    print("\r", "Time Step: ", time_step, end="\n")

print("\r", "Qow: ", b3d._QowTS, end="\n")
# for time_step in range(1, 10):  # b3d._TMAX):
#     b3d.update()  # update the model's main time loop
#     b3d.update_dune_domain()  # now update the dune domain and increase time by one year
#     int_domain = np.flip(b3d.InteriorDomain)
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     # fig1, (ax1, ax3) = plt.subplots(1, 2, sharey=True)
#     mat = ax.matshow(
#         int_domain,
#         cmap="terrain",
#         vmin=-0.25, vmax=0.25,
#     )
#     fig.colorbar(mat)
#     ax.set_title("Initial Elevation $(dam)$")
#     ax.set_ylabel("barrier width (dam)")
#     ax.set_xlabel("barrier length (dam)")
#     plt.gca().xaxis.tick_bottom()
#     # Print time step to screen
#     print("\r", "Time Step: ", time_step, end="\n")

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     b3d.InteriorDomain, cmap="terrain",
#     vmin=-0.3, vmax=0.4
# )
# fig2.colorbar(mat)
# ax2.set_title("Initial Elevation $(dam)$")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.show()

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
