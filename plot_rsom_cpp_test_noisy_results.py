import numpy as np

folder = "results/n5_mindeg2_sigma00_20250205_1707_59/"

testid = folder[8: 27]

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
file_rot_errs = folder + testid + "rot_errors_mean.txt"
print("file_rot_errs")
print(file_rot_errs)
with open(file_rot_errs, 'r') as f:
    for count, line in enumerate(f, start=1):
        if count % 2 == 0:
            # print(line)
            rot_errs = np.append(rot_errs, float(line))

rot_errs_mean = np.average(rot_errs)
print("rot_errs_mean")
print(rot_errs_mean)

#transl errs
transl_errs = np.array([], dtype=np.float64)
file_transl_errs = folder + testid + "transl_errors_mean.txt"
print("file_transl_errs")
print(file_transl_errs)
with open(file_transl_errs, 'r') as f:
    for count, line in enumerate(f, start=1):
        if count % 2 == 0:
            # print(line)
            transl_errs = np.append(transl_errs, float(line))

transl_errs_mean = np.average(transl_errs)
print("transl_errs_mean")
print(transl_errs_mean)

#exec time
exec_times = np.array([], dtype=np.float64)
file_exec_times = folder + testid + "exec_times.txt"
print("file_exec_times")
print(file_exec_times)
with open(file_exec_times, 'r') as f:
    for count, line in enumerate(f, start=1):
        if count % 2 == 0:
            # print(line)
            exec_times = np.append(exec_times, float(line))

exec_times_mean = np.average(exec_times)
print("exec_times_mean")
print(exec_times_mean)