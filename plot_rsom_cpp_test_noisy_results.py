import os

import numpy as np



rsom_rs_results_path = "results/escape_spectra/"
folders = os.listdir(rsom_rs_results_path)
tuples_rs = []

for f in folders:

    print("f")
    print(f)

    folder = rsom_rs_results_path + f

    print("folder")
    print(folder)

    testid = folder[8 + len("escape_spectra/"): 27 +len("escape_spectra/")]

    print("testid")
    print(testid)

    testid_split = testid.split("_")
    n = int(testid_split[0][1:])
    mindeg = int(testid_split[1][6:])
    sigma = float(testid_split[2][5:]) / 10 #!! /10
    print("n")
    print(n)
    print("mindeg")
    print(mindeg)
    print("sigma")
    print(sigma)

    #rot errs
    rot_errs = np.array([], dtype=np.float64)
    if (testid[-1] != "_"):
        testid = testid + "_" #temporary fix
    file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
    print("file_rot_errs")
    print(file_rot_errs)
    with open(file_rot_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                rot_errs = np.append(rot_errs, float(line))

    rot_errs_mean_rs = np.average(rot_errs)
    print("rot_errs_mean_rs")
    print(rot_errs_mean_rs)

    #transl errs
    transl_errs = np.array([], dtype=np.float64)
    file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
    print("file_transl_errs")
    print(file_transl_errs)
    with open(file_transl_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                transl_errs = np.append(transl_errs, float(line))

    transl_errs_mean_rs = np.average(transl_errs)
    print("transl_errs_mean_rs")
    print(transl_errs_mean_rs)

    #exec time
    exec_times = np.array([], dtype=np.float64)
    file_exec_times = folder + "/" + testid + "exec_times.txt"
    print("file_exec_times")
    print(file_exec_times)
    with open(file_exec_times, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                exec_times = np.append(exec_times, float(line))

    exec_times_mean_rs = np.average(exec_times)
    print("exec_times_mean_rs")
    print(exec_times_mean_rs)

    tuple_rs = (n, mindeg, sigma, rot_errs_mean_rs, transl_errs_mean_rs, exec_times_mean_rs)
    tuples_rs.append(tuple_rs)


############################

rsom_icp_results_path = "results_icp/"
folders = os.listdir(rsom_icp_results_path)

tuples_icp = []

for f in folders:

    print("f")
    print(f)

    folder = rsom_icp_results_path + f

    print("folder")
    print(folder)

    testid = folder[12: 31]

    print("testid")
    print(testid)

    testid_split = testid.split("_")
    n = int(testid_split[0][1:])
    mindeg = int(testid_split[1][6:])
    sigma = float(testid_split[2][5:])
    print("n")
    print(n)
    print("mindeg")
    print(mindeg)
    print("sigma")
    print(sigma)

    #rot errs
    rot_errs = np.array([], dtype=np.float64)
    if (testid[-1] != "_"):
        testid = testid + "_" #temporary fix
    file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
    print("file_rot_errs")
    print(file_rot_errs)
    with open(file_rot_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                rot_errs = np.append(rot_errs, float(line))

    rot_errs_mean_icp = np.average(rot_errs)
    print("rot_errs_mean_icp")
    print(rot_errs_mean_icp)

    #transl errs
    transl_errs = np.array([], dtype=np.float64)
    file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
    print("file_transl_errs")
    print(file_transl_errs)
    with open(file_transl_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                transl_errs = np.append(transl_errs, float(line))

    transl_errs_mean_icp = np.average(transl_errs)
    print("transl_errs_mean_icp")
    print(transl_errs_mean_icp)

    #exec time
    exec_times = np.array([], dtype=np.float64)
    file_exec_times = folder + "/" + testid + "exec_times.txt"
    print("file_exec_times")
    print(file_exec_times)
    with open(file_exec_times, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                exec_times = np.append(exec_times, float(line))

    exec_times_mean_icp = np.average(exec_times)
    print("exec_times_mean_icp")
    print(exec_times_mean_icp)

    tuple_icp = (n, mindeg, sigma, rot_errs_mean_icp, transl_errs_mean_icp, exec_times_mean_icp)
    tuples_icp.append(tuple_icp)

############################


rsom_procrustes_results_path = "results_procrustes/"
folders = os.listdir(rsom_procrustes_results_path)

tuples_procrustes = []

for f in folders:

    print("f")
    print(f)

    folder = rsom_procrustes_results_path + f

    print("folder")
    print(folder)

    testid = folder[19: 38]

    print("testid")
    print(testid)

    testid_split = testid.split("_")
    n = int(testid_split[0][1:])
    mindeg = int(testid_split[1][6:])
    sigma = float(testid_split[2][5:])
    print("n")
    print(n)
    print("mindeg")
    print(mindeg)
    print("sigma")
    print(sigma)

    #rot errs
    rot_errs = np.array([], dtype=np.float64)
    if (testid[-1] != "_"):
        testid = testid + "_" #temporary fix
    file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
    print("file_rot_errs")
    print(file_rot_errs)
    with open(file_rot_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                rot_errs = np.append(rot_errs, float(line))

    rot_errs_mean_procrustes = np.average(rot_errs)
    print("rot_errs_mean_procrustes")
    print(rot_errs_mean_procrustes)

    #transl errs
    transl_errs = np.array([], dtype=np.float64)
    file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
    print("file_transl_errs")
    print(file_transl_errs)
    with open(file_transl_errs, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                transl_errs = np.append(transl_errs, float(line))

    transl_errs_mean_procrustes = np.average(transl_errs)
    print("transl_errs_mean_procrustes")
    print(transl_errs_mean_procrustes)

    #exec time
    exec_times = np.array([], dtype=np.float64)
    file_exec_times = folder + "/" + testid + "exec_times.txt"
    print("file_exec_times")
    print(file_exec_times)
    with open(file_exec_times, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                exec_times = np.append(exec_times, float(line))

    exec_times_mean_procrustes = np.average(exec_times)
    print("exec_times_mean_procrustes")
    print(exec_times_mean_procrustes)

    tuple_procrustes = (n, mindeg, sigma, rot_errs_mean_procrustes, transl_errs_mean_procrustes, exec_times_mean_procrustes)
    tuples_procrustes.append(tuple_procrustes)

############################ PLOT ############################

import matplotlib.pyplot as plt

mindeg_to_plot = 2

xpoints = np.array([0, 0.1, 0.2, 0.5, 1, 2])

############################ PLOT RS ############################


rot_errs_rs = np.zeros_like(xpoints, dtype=np.float64)
transl_errs_rs = np.zeros_like(xpoints, dtype=np.float64)
exec_times_rs = np.zeros_like(xpoints, dtype=np.float64)

for t in tuples_rs:
    t_n = t[0]
    t_mindeg = t[1]
    if (t_mindeg != mindeg_to_plot):
        continue
    t_sigma = t[2]
    t_rot_errs_mean_rs = t[3]
    t_transl_errs_mean_rs = t[4]
    t_exec_times_mean_rs = t[5]

    sigma_index = np.where(xpoints == t_sigma)
    rot_errs_rs[sigma_index] = t_rot_errs_mean_rs
    transl_errs_rs[sigma_index] = t_transl_errs_mean_rs
    exec_times_rs[sigma_index] = t_exec_times_mean_rs

############################ PLOT ICP ############################

rot_errs_icp = np.zeros_like(xpoints, dtype=np.float64)
transl_errs_icp = np.zeros_like(xpoints, dtype=np.float64)
exec_times_icp = np.zeros_like(xpoints, dtype=np.float64)

for t in tuples_icp:
    t_n = t[0]
    t_mindeg = t[1]
    if (t_mindeg != mindeg_to_plot):
        continue
    t_sigma = t[2]
    t_rot_errs_mean_icp = t[3]
    t_transl_errs_mean_icp = t[4]
    t_exec_times_mean_icp = t[5]

    sigma_index = np.where(xpoints == t_sigma)
    rot_errs_icp[sigma_index] = t_rot_errs_mean_icp
    transl_errs_icp[sigma_index] = t_transl_errs_mean_icp
    exec_times_icp[sigma_index] = t_exec_times_mean_icp

############################ PLOT PROCRUSTES ############################

rot_errs_procrustes = np.zeros_like(xpoints, dtype=np.float64)
transl_errs_procrustes = np.zeros_like(xpoints, dtype=np.float64)
exec_times_procrustes = np.zeros_like(xpoints, dtype=np.float64)

for t in tuples_procrustes:
    t_n = t[0]
    t_mindeg = t[1]
    if (t_mindeg != mindeg_to_plot):
        continue
    t_sigma = t[2]
    t_rot_errs_mean_procrustes = t[3]
    t_transl_errs_mean_procrustes = t[4]
    t_exec_times_mean_procrustes = t[5]

    sigma_index = np.where(xpoints == t_sigma)
    rot_errs_procrustes[sigma_index] = t_rot_errs_mean_procrustes
    transl_errs_procrustes[sigma_index] = t_transl_errs_mean_procrustes
    exec_times_procrustes[sigma_index] = t_exec_times_mean_procrustes


fig = plt.figure()
gs = fig.add_gridspec(3, hspace=1)
axs = gs.subplots(sharex=True)
fig.suptitle('mindeg ' + str(mindeg_to_plot))
axs[0].plot(xpoints, rot_errs_rs, 'o', label="TPGO-RS", color='green')
axs[0].set_title('Rotation errors')
# axs[0].set_ylim(-0.1, 3.14)
axs[1].plot(xpoints, transl_errs_rs, 'o', label="TPGO-RS", color='green')
axs[1].set_title('Translation errors')

axs[2].plot(xpoints, exec_times_rs, 'o', label="TPGO-RS", color='green')
axs[2].set_title('Execution times [ms]')
# axs[2].set_ylim()

axs[0].plot(xpoints, rot_errs_icp, 'o', label="TPGO-ICP", color='red')
axs[0].set_title('Rotation errors')
axs[1].plot(xpoints, transl_errs_icp, 'o', label="TPGO-ICP", color='red')
axs[1].set_title('Translation errors')

axs[2].plot(xpoints, exec_times_icp, 'o', label="TPGO-ICP", color='red')
axs[2].set_title('Execution times [ms]')

axs[0].plot(xpoints, rot_errs_procrustes, 'o', label="TPGO-PROCR", color='blue')
axs[0].set_title('Rotation errors')
axs[1].plot(xpoints, transl_errs_procrustes, 'o', label="TPGO-PROCR", color='blue')
axs[1].set_title('Translation errors')
axs[2].plot(xpoints, exec_times_procrustes, 'o', label="TPGO-PROCR", color='blue')
axs[2].set_title('Execution times [ms]')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
    ax.legend()

plt.show()

