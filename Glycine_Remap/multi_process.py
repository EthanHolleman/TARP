import subprocess
import sys
import os
import multiprocessing as mp
from multiprocessing import Pool

import tqdm


def run_process(process):
    subprocess.call(process)

def mute():
    sys.stdout = open(os.devnull, 'w')

def run_pool(processes, cpu_count=mp.cpu_count()):
    # Process list should be a list of rady to go subprocess commands
    # so should be a list of lists
    pool = Pool(processes=cpu_count, initializer=mute)
    print('{} remaps to run'.format(len(processes)))
    print('Starting run with {} workers'.format(cpu_count))
    for _ in tqdm.tqdm(pool.imap_unordered(run_process, processes), total=len(processes)):
        pass
