import subprocess
import sys
import os
import multiprocessing as mp
from multiprocessing import Pool
import tqdm

from Transposer.remap import transpose

def mute():
    sys.stdout = open(os.devnull, 'w')

def run_pool(processes, cpu_count=mp.cpu_count()):
    # Process list should be a list of rady to go subprocess commands
    # so should be a list of lists
    pool = Pool(processes=cpu_count, initializer=mute)
    print('{} remaps to run'.format(len(processes)))
    print('Starting run with {} workers'.format(cpu_count))

    for _ in tqdm.tqdm(pool.starmap(transpose, processes), total=len(processes)):
        pass
# example of how this should work
# processes shuold be a list of lists where instead of being subprocess
# formated will be ordered parameters that are passed into the run process function
# basic layout
# gather arguement keywords iteratively for all files
# process those into list or tuple
# unpack and pass those into remap through the run_process function
'''
transpose(BT_index='/media/ethan/Vault/Trans_Gypsy_Files/Gmax2.1_Bowtie/Gmax2.1_Bowtie_Index',
          solo_con='/home/ethan/Documents/Final_GM_Cons/consensus_Gmr364SOLO.fna',
          intact_con='/home/ethan/Documents/Final_GM_Cons/consensus_Gmr364INTACT.fna',
          cur_acc='/media/ethan/Vault/Trans_Gypsy_Files/Chr2Acc/chr2acc_Gmax2.1',
          old_acc='/media/ethan/Vault/Trans_Gypsy_Files/Chr2Acc/chr2acc_Gmax1.0',
          output='/home/ethan/Documents/transpose_test', name='GMR 364',
          cur_BDB='/media/ethan/Vault/Trans_Gypsy_Files/Gmax2.1_Blast_DB/Gmax2.1_Blast_DB',
          old_BDB='/media/ethan/Vault/Trans_Gypsy_Files/Gmax1.0_Blast_DB/Gmax1.0_Blast_DB',
          old_elements='/media/ethan/Vault/Trans_Gypsy_Files/Gypsy_Fams/Gmr364.fna',
          allowance=150, verbose=True)
'''
