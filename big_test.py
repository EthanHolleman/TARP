DIR = '/media/ethan/EH_DATA/Gypsy_Seperated'
BTI = '/media/ethan/EH_DATA/Gmax2.1_Bowtie/Gmax_2.1_BTI'
BDB_OLD, BDB_NEW = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/GMAX_1.0_BDB', '/media/ethan/EH_DATA/GMAX_2.1_BDB_parsed/GM_2.1_BD'
ACC_OLD, ACC_NEW = '/media/ethan/EH_DATA/GMax1.1_assembly/chr2acc', '/media/ethan/EH_DATA/GMAX_2.1_BDB_parsed/chr2acc.txt'
OUT = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_2'

import os
import os

def get_family_groups():
    files = os.listdir(DIR)
    fam_dict = {}
    for f in files:
        temp = f[3:]
        if 'INTACT' in f:
            temp = temp.split('INTACT')[0]
        else:
            temp = temp.split('SOLO')[0]
        if temp in fam_dict:
            fam_dict[temp].append(f)
        else:
            fam_dict[temp] = [f]

    return fam_dict


def make_cmd(run_name, solo, intact):
    cmd = ['python3', 'TARP.py', '-I', intact, '-S', solo, '-P', BDB_OLD, '-C',
            BDB_NEW, '-acc_c', ACC_NEW, '-acc_o', ACC_OLD, '-name', run_name,
            '-B', BTI, '-O', OUT]
    return ' '.join(cmd)

def exec_cmd(cmd):
    call = os.system(cmd)
    if call == 0:
        return True
    else:
        return False

def run_all():
    families = get_family_groups()
    for k, v in families.items():
        rn = 'GMR_' + str(k)
        if len(v) == 1:
            continue
        if 'SOLO' in v[0]:
            solo, intact = v[0], v[1]

        else:
            solo, intact = v[1], v[0]
        solo, intact = os.path.join(DIR, solo), os.path.join(DIR, intact)

        cmd = make_cmd(rn, solo, intact)
        print(cmd)
        exec_cmd(cmd)

run_all()
