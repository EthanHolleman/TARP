import matplotlib.pyplot as plt
import os



TOP = '/media/ethan/EH_DATA/TARP_Runs/Test_Min_Dist_log'

dirs = [os.path.join(TOP, subdir) for subdir in os.listdir(TOP)]
result_dirs = [os.path.join(subdir, 'Results') for subdir in dirs]
log_files = []

for r in result_dirs:
    result_files = os.listdir(r)
    for f in result_files:
        if 'log' in f:
            log_files.append(os.path.join(r, f))
data = []
for log in log_files:
    with open(log) as l:
        data.append(l.readline().strip().split('\t'))

for d in data:
    if len(d) > 2:
        print('o')
