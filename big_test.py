DIR = '/media/ethan/EH_DATA/Gypsy_Seperated'
BTI = '/media/ethan/EH_DATA/Gmax2.1_Bowtie/GM_2.1_BTI'
BDB_OLD, BDB_NEW = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/GMAX_1.0_BDB', '/media/ethan/EH_DATA/GMAX_2.1_BDB_parsed/GM_2.1_BDB'
ACC_OLD, ACC_NEW = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/chr2acc', '/media/ethan/EH_DATA/GMAX_2.1_BDB_parsed/chr2acc'
OUT = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_12-26'
LOG = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_12-26/log.txt'
import os
import time

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


def make_cmd(run_name, intact, solo=False):
    if solo:
        cmd = ['python3', 'TARP', '-I', intact, '-S', solo, '-P', BDB_OLD, '-C',
                BDB_NEW, '-acc_c', ACC_NEW, '-acc_o', ACC_OLD, '-name', run_name,
                '-B', BTI, '-O', OUT]
    else:
        cmd = ['python3', 'TARP', '-I', intact, '-P', BDB_OLD, '-C',
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
    with open(LOG, 'w') as log:
        for k, v in families.items():
            rn = 'GMR_' + str(k)
            if len(v) > 1:
                if 'SOLO' in v[0]:
                    solo, intact = v[0], v[1]
                elif 'INTACT' in v[0]:
                    solo, intact = v[1], v[0]
            elif len(v) == 1:
                if 'INTACT' in v[0]:
                    solo, intact = None, v[0]
            else:
                continue
        
            if solo:
                solo = os.path.join(DIR, solo)
            intact =  os.path.join(DIR, intact)
            cmd = make_cmd(rn, intact, solo)
            log.write(cmd + '\n')
            exec_cmd(cmd)
start = time.time()
run_all()
end = time.time()
total_time = end - start
with open(LOG, 'a+') as log:
    log.write('Start: {}\nEnd: {}\nTotal: {}'.format(str(start), str(end), str(total_time)))
