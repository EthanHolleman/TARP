import subprocess
from multiprocessing import Pool


def run_process(process):
    subprocess.call(process)

def run_pool(processes):
    # Process list should be a list of rady to go subprocess commands
    # so should be a list of lists
    pool = Pool(processes=len(processes))
    pool.map(run_process, processes)
