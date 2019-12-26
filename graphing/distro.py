import csv
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

TOP = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_3'
SUMM = '/media/ethan/EH_DATA/Element_Summaries/IILTRGypsy(1).txt'

def old_el_reader(summary):
    # read soybase summary file for superfam and return fam dict
    # fam dict keys = family names values = lists [0] == number of solo
    # elements list[1] == number of intact elements
    fam_dict = {}
    with open(summary) as s:
        reader = csv.reader(s, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            try:
                fam = row[6].upper()
                if fam not in fam_dict:
                    fam_dict[fam] = [0, 0]
                if row[7] == 'INTACT':
                    fam_dict[fam][1] += 1
                elif row[7] == 'SOLO':
                    fam_dict[fam][0] += 1
            except IndexError as e:
                continue

    return fam_dict

print(len(old_el_reader(SUMM)))


def get_remap_files(top_dir):
    fams = [os.path.join(top_dir, fam) for fam in os.listdir(top_dir)]
    for f in fams:
        results = os.path.join(f, 'Results')
        results_files = os.listdir(results)
        for r in results_files:
            if r[-3:] == 'csv':
                yield os.path.join(results, r)


def family_sums(remap_files):
    fam_dict = {}
    for remap_file in remap_files:
        with open(remap_file) as f:
            reader = csv.reader(f)
            fam = os.path.basename(remap_file)[:-4]
            fam_dict[fam.replace('_', '')] = [0, 0]  # solo, intact
            for row in reader:
                if row[5] == 'S':
                    fam_dict[fam.replace('_', '')][0] += 1
                elif row[5] == 'I':
                    fam_dict[fam.replace('_', '')][1] += 1
    return fam_dict


def diff_dict(remap_dict, old_dict):
    diff_dict = {}
    for fam in old_dict:
        try:
            s, i = remap_dict[fam][0] - old_dict[fam][0], remap_dict[fam][1] - old_dict[fam][1]
            diff_dict[fam] = [s, i]
        except KeyError as e:
            print('key', fam)
            continue
    return diff_dict


old_dict = old_el_reader(SUMM)
remap_dict = family_sums(get_remap_files(TOP))
print(len(old_dict), len(remap_dict))
#diff = diff_dict(remap_dict, old_dict)



def stack_bars(f_dict, solo_col='skyblue', intact_col='tan', width=1):
    SOLO_COL = solo_col
    INTACT_COL = intact_col

    groups, N = [k for k in f_dict], len(f_dict)
    solo_values = [f_dict[f][0] for f in groups]
    intact_values = [f_dict[f][1] for f in groups]

    solo_patch = mpatches.Patch(color=SOLO_COL, label='Solo')
    intact_patch = mpatches.Patch(color=INTACT_COL, label='Intact')
    plt.legend(handles=[solo_patch, intact_patch], loc='upper right')


    ind = np.arange(N)
    plt.xticks(ind, groups, fontsize=5, rotation=90)
    p_solo = plt.bar(ind, solo_values, width, color='tan')
    p_intact = plt.bar(ind, intact_values, width, bottom=solo_values, color='skyblue')
    plt.ylabel('Total Elements')
    plt.title('Remapped Element Totals by Family')

    plt.show()

#stack_bars(diff)
