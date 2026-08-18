import os

import numpy as np
from matplotlib.colors import ListedColormap


rsom_rs_results_path = "~/Downloads/cpp"
folders = os.listdir(rsom_rs_results_path)
tuples_rs = []
tuples_lambdas_acceptable = []

for f in folders:
    print("f")
    print(f)

    folder = ssom_rs_results_path + "/" + f

    print("folder")
    print(folder)

    len_timestamp = 17
    testid = folder[len(ssom_rs_results_path) + 1 : len(folder) - len_timestamp]

    print("testid")
    print(testid)

    testid_split = testid.split("_")
    n = int(testid_split[0][1:])
    mindeg = int(testid_split[1][6:])
    sigma = float(testid_split[2][5:]) / 100  #!! /100
    print("n")
    print(n)
    print("mindeg")
    print(mindeg)
    print("sigma")
    print(sigma)

    # rot errs
    rot_errs = np.array([], dtype=np.float64)
    if testid[-1] != "_":
        testid = testid + "_"  # temporary fix
    file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
    print("file_rot_errs")
    print(file_rot_errs)
    with open(file_rot_errs, "r") as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                rot_errs = np.append(rot_errs, float(line))

    rot_errs_mean_rs = np.average(rot_errs)
    print("rot_errs_mean_rs")
    print(rot_errs_mean_rs)

    # transl errs
    transl_errs = np.array([], dtype=np.float64)
    file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
    print("file_transl_errs")
    print(file_transl_errs)
    with open(file_transl_errs, "r") as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                transl_errs = np.append(transl_errs, float(line))

    transl_errs_mean_rs = np.average(transl_errs)
    print("transl_errs_mean_rs")
    print(transl_errs_mean_rs)

    # scale errs
    scale_errs = np.array([], dtype=np.float64)
    file_scale_errs = folder + "/" + testid + "lambda_errors_mean.txt"
    print("file_scale_errs")
    print(file_scale_errs)
    with open(file_scale_errs, "r") as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                scale_errs = np.append(scale_errs, float(line))

    scale_errs_mean_rs = np.average(scale_errs)
    print("scale_errs_mean_rs")
    print(scale_errs_mean_rs)

    # exec time
    exec_times = np.array([], dtype=np.float64)
    file_exec_times = folder + "/" + testid + "exec_times.txt"
    print("file_exec_times")
    print(file_exec_times)
    with open(file_exec_times, "r") as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                exec_times = np.append(exec_times, float(line))

    exec_times_mean_rs = np.average(exec_times)
    print("exec_times_mean_rs")
    print(exec_times_mean_rs)

    tuple_rs = (
        n,
        mindeg,
        sigma,
        rot_errs_mean_rs,
        transl_errs_mean_rs,
        scale_errs_mean_rs,
        exec_times_mean_rs,
    )
    tuples_rs.append(tuple_rs)

    ### plot acceptable lambdas
    lambdas_acceptable = np.array([], dtype=np.float64)
    file_lambdas_acceptable = folder + "/" + testid + "lambdas_acceptable.txt"
    print("file_lambdas_acceptable")
    print(file_lambdas_acceptable)
    with open(file_lambdas_acceptable, "r") as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                # print(line)
                lambdas_acceptable = np.append(lambdas_acceptable, float(line))
    print("lambdas_acceptable")
    print(lambdas_acceptable)

    # Store the acceptable lambdas
    tuple_lambdas_acceptable = (n, mindeg, sigma, lambdas_acceptable)
    tuples_lambdas_acceptable.append(tuple_lambdas_acceptable)

############################

ssom_icp_results_path = "results_icp/"
# folders = os.listdir(ssom_icp_results_path)

# tuples_icp = []

# for f in folders:

#     print("f")
#     print(f)

#     folder = ssom_icp_results_path + f

#     print("folder")
#     print(folder)

#     testid = folder[12: 31]

#     print("testid")
#     print(testid)

#     testid_split = testid.split("_")
#     n = int(testid_split[0][1:])
#     mindeg = int(testid_split[1][6:])
#     sigma = float(testid_split[2][5:])
#     print("n")
#     print(n)
#     print("mindeg")
#     print(mindeg)
#     print("sigma")
#     print(sigma)

#     #rot errs
#     rot_errs = np.array([], dtype=np.float64)
#     if (testid[-1] != "_"):
#         testid = testid + "_" #temporary fix
#     file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
#     print("file_rot_errs")
#     print(file_rot_errs)
#     with open(file_rot_errs, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 rot_errs = np.append(rot_errs, float(line))

#     rot_errs_mean_icp = np.average(rot_errs)
#     print("rot_errs_mean_icp")
#     print(rot_errs_mean_icp)

#     #transl errs
#     transl_errs = np.array([], dtype=np.float64)
#     file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
#     print("file_transl_errs")
#     print(file_transl_errs)
#     with open(file_transl_errs, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 transl_errs = np.append(transl_errs, float(line))

#     transl_errs_mean_icp = np.average(transl_errs)
#     print("transl_errs_mean_icp")
#     print(transl_errs_mean_icp)

#     #exec time
#     exec_times = np.array([], dtype=np.float64)
#     file_exec_times = folder + "/" + testid + "exec_times.txt"
#     print("file_exec_times")
#     print(file_exec_times)
#     with open(file_exec_times, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 exec_times = np.append(exec_times, float(line))

#     exec_times_mean_icp = np.average(exec_times)
#     print("exec_times_mean_icp")
#     print(exec_times_mean_icp)

#     tuple_icp = (n, mindeg, sigma, rot_errs_mean_icp, transl_errs_mean_icp, exec_times_mean_icp)
#     tuples_icp.append(tuple_icp)

############################


ssom_procrustes_results_path = "results_procrustes/"
# folders = os.listdir(ssom_procrustes_results_path)

# tuples_procrustes = []

# for f in folders:

#     print("f")
#     print(f)

#     folder = ssom_procrustes_results_path + f

#     print("folder")
#     print(folder)

#     testid = folder[19: 38]

#     print("testid")
#     print(testid)

#     testid_split = testid.split("_")
#     n = int(testid_split[0][1:])
#     mindeg = int(testid_split[1][6:])
#     sigma = float(testid_split[2][5:])
#     print("n")
#     print(n)
#     print("mindeg")
#     print(mindeg)
#     print("sigma")
#     print(sigma)

#     #rot errs
#     rot_errs = np.array([], dtype=np.float64)
#     if (testid[-1] != "_"):
#         testid = testid + "_" #temporary fix
#     file_rot_errs = folder + "/" + testid + "rot_errors_mean.txt"
#     print("file_rot_errs")
#     print(file_rot_errs)
#     with open(file_rot_errs, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 rot_errs = np.append(rot_errs, float(line))

#     rot_errs_mean_procrustes = np.average(rot_errs)
#     print("rot_errs_mean_procrustes")
#     print(rot_errs_mean_procrustes)

#     #transl errs
#     transl_errs = np.array([], dtype=np.float64)
#     file_transl_errs = folder + "/" + testid + "transl_errors_mean.txt"
#     print("file_transl_errs")
#     print(file_transl_errs)
#     with open(file_transl_errs, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 transl_errs = np.append(transl_errs, float(line))

#     transl_errs_mean_procrustes = np.average(transl_errs)
#     print("transl_errs_mean_procrustes")
#     print(transl_errs_mean_procrustes)

#     #exec time
#     exec_times = np.array([], dtype=np.float64)
#     file_exec_times = folder + "/" + testid + "exec_times.txt"
#     print("file_exec_times")
#     print(file_exec_times)
#     with open(file_exec_times, 'r') as f:
#         for count, line in enumerate(f, start=1):
#             if count % 2 == 0:
#                 # print(line)
#                 exec_times = np.append(exec_times, float(line))

#     exec_times_mean_procrustes = np.average(exec_times)
#     print("exec_times_mean_procrustes")
#     print(exec_times_mean_procrustes)

#     tuple_procrustes = (n, mindeg, sigma, rot_errs_mean_procrustes, transl_errs_mean_procrustes, exec_times_mean_procrustes)
#     tuples_procrustes.append(tuple_procrustes)

############################ PLOT ############################

import matplotlib.pyplot as plt

mindeg_to_plot = 3

xpoints = np.array([0.0, 0.01, 0.02, 0.05,0.01,0.02,0.05,0.1,0.2,0.5,1,2])

############################ PLOT RS ############################


rot_errs_rs = np.zeros_like(xpoints, dtype=np.float64)
transl_errs_rs = np.zeros_like(xpoints, dtype=np.float64)
scale_errs_rs = np.zeros_like(xpoints, dtype=np.float64)
exec_times_rs = np.zeros_like(xpoints, dtype=np.float64)

num_tests_per_instance = 5  # must match what was used in the C++ tests
print("xpoints.shape[0]")
print(xpoints.shape[0])
lambdas_acceptable = np.zeros([xpoints.shape[0], num_tests_per_instance], dtype=np.float64)  # will be binary (0 or 1) for plotting

for t in tuples_rs:
    t_n = t[0]
    t_mindeg = t[1]
    if t_mindeg != mindeg_to_plot:
        continue
    t_sigma = t[2]
    t_rot_errs_mean_rs = t[3]
    t_transl_errs_mean_rs = t[4]
    t_scale_errs_mean_rs = t[5]
    t_exec_times_mean_rs = t[6]

    sigma_index = np.where(xpoints == t_sigma)
    rot_errs_rs[sigma_index] = t_rot_errs_mean_rs
    transl_errs_rs[sigma_index] = t_transl_errs_mean_rs
    scale_errs_rs[sigma_index] = t_scale_errs_mean_rs
    exec_times_rs[sigma_index] = t_exec_times_mean_rs
    

n_tmp = 5
for t in enumerate(tuples_lambdas_acceptable):
    t_n = t[0]
    t_mindeg = t[1]
    if t_n != n_tmp:
        continue
    if t_mindeg != mindeg_to_plot:
        continue
    t_sigma = t[2]

    sigma_index = np.where(xpoints == t_sigma)
    lambdas_acceptable[sigma_index, :] = t[3]  # Assuming this is already binary (0 or 1) for each sigma   
    
print("lambdas_acceptable")
print(lambdas_acceptable)

# Build a binary matrix for acceptable lambdas and plot as red/green squares.
rows = []
for idx, t in enumerate(tuples_lambdas_acceptable):
    if t[1] == mindeg_to_plot:
        rows.append((t[0], idx))  # (n, original_index)
rows.sort(key=lambda x: x[0])

# Build binary matrix from tuples_lambdas_acceptable filtered & sorted in `rows`
num_rows = len(rows)
if num_rows == 0:
    print("No acceptable-lambda rows to plot for mindeg", mindeg_to_plot)
else:
    mat = np.zeros((num_rows, num_tests_per_instance), dtype=int)
    ylabels = []
    for i, (n_val, orig_idx) in enumerate(rows):
        t = tuples_lambdas_acceptable[orig_idx]
        arr = np.asarray(t[3], dtype=int)
        # Ensure length matches expected number of tests
        if arr.size != num_tests_per_instance:
            tmp = np.zeros(num_tests_per_instance, dtype=int)
            tmp[: min(arr.size, num_tests_per_instance)] = arr[: min(arr.size, num_tests_per_instance)]
            arr = tmp
        mat[i, :] = arr
        ylabels.append(f"n={n_val}, σ={t[2]}")

    # Plot binary matrix: 0 -> red, 1 -> green
    fig_table, ax_table = plt.subplots(figsize=(max(6, num_tests_per_instance), max(3, num_rows * 0.3)))
    cmap = ListedColormap(["red", "green"])
    im = ax_table.imshow(mat, cmap=cmap, aspect="auto", interpolation="nearest", vmin=0, vmax=1)

    ax_table.set_yticks(np.arange(num_rows))
    ax_table.set_yticklabels(ylabels, fontsize=8)
    ax_table.set_xticks(np.arange(num_tests_per_instance))
    ax_table.set_xticklabels([f"run {i+1}" for i in range(num_tests_per_instance)])
    ax_table.set_title("Acceptable lambdas (red=0, green=1)")
    plt.tight_layout()

############################ PLOT ICP ############################

# rot_errs_icp = np.zeros_like(xpoints, dtype=np.float64)
# transl_errs_icp = np.zeros_like(xpoints, dtype=np.float64)
# exec_times_icp = np.zeros_like(xpoints, dtype=np.float64)

# for t in tuples_icp:
#     t_n = t[0]
#     t_mindeg = t[1]
#     if (t_mindeg != mindeg_to_plot):
#         continue
#     t_sigma = t[2]
#     t_rot_errs_mean_icp = t[3]
#     t_transl_errs_mean_icp = t[4]
#     t_exec_times_mean_icp = t[5]

#     sigma_index = np.where(xpoints == t_sigma)
#     rot_errs_icp[sigma_index] = t_rot_errs_mean_icp
#     transl_errs_icp[sigma_index] = t_transl_errs_mean_icp
#     exec_times_icp[sigma_index] = t_exec_times_mean_icp

############################ PLOT PROCRUSTES ############################

# rot_errs_procrustes = np.zeros_like(xpoints, dtype=np.float64)
# transl_errs_procrustes = np.zeros_like(xpoints, dtype=np.float64)
# exec_times_procrustes = np.zeros_like(xpoints, dtype=np.float64)

# for t in tuples_procrustes:
#     t_n = t[0]
#     t_mindeg = t[1]
#     if (t_mindeg != mindeg_to_plot):
#         continue
#     t_sigma = t[2]
#     t_rot_errs_mean_procrustes = t[3]
#     t_transl_errs_mean_procrustes = t[4]
#     t_exec_times_mean_procrustes = t[5]

#     sigma_index = np.where(xpoints == t_sigma)
#     rot_errs_procrustes[sigma_index] = t_rot_errs_mean_procrustes
#     transl_errs_procrustes[sigma_index] = t_transl_errs_mean_procrustes
#     exec_times_procrustes[sigma_index] = t_exec_times_mean_procrustes


fig = plt.figure()
gs = fig.add_gridspec(3, hspace=1)
axs = gs.subplots(sharex=True)
fig.suptitle("mindeg " + str(mindeg_to_plot))
axs[0].plot(xpoints, rot_errs_rs, "o", label="TPGO-RS", color="green")
axs[0].set_title("Rotation errors")
# axs[0].set_ylim(-0.1, 3.14)
axs[1].plot(xpoints, transl_errs_rs, "o", label="TPGO-RS", color="green")
axs[1].set_title("Translation errors")

axs[2].plot(xpoints, exec_times_rs, "o", label="TPGO-RS", color="green")
axs[2].set_title("Execution times [ms]")
# axs[2].set_ylim()

# axs[0].plot(xpoints, rot_errs_icp, 'o', label="TPGO-ICP", color='red')
# axs[0].set_title('Rotation errors')
# axs[1].plot(xpoints, transl_errs_icp, 'o', label="TPGO-ICP", color='red')
# axs[1].set_title('Translation errors')

# axs[2].plot(xpoints, exec_times_icp, 'o', label="TPGO-ICP", color='red')
# axs[2].set_title('Execution times [ms]')

# axs[0].plot(xpoints, rot_errs_procrustes, 'o', label="TPGO-PROCR", color='blue')
# axs[0].set_title('Rotation errors')
# axs[1].plot(xpoints, transl_errs_procrustes, 'o', label="TPGO-PROCR", color='blue')
# axs[1].set_title('Translation errors')
# axs[2].plot(xpoints, exec_times_procrustes, 'o', label="TPGO-PROCR", color='blue')
# axs[2].set_title('Execution times [ms]')



# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
    ax.legend()

    

plt.show()
