import csv
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

TOP = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_12-30_AMB'
SUMM = '/media/ethan/EH_DATA/Element_Summaries/IILTRGypsy.txt'

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
                fam = row[6].upper().replace('/', '')
                if fam not in fam_dict:
                    fam_dict[fam] = [0, 0]  # solo and intact count
                if row[7] == 'INTACT':
                    fam_dict[fam][1] += 1
                elif row[7] == 'SOLO':
                    fam_dict[fam][0] += 1
            except IndexError as e:
                continue

    return fam_dict



def get_remap_files(top_dir):
    csv_files = []
    fams = [os.path.join(top_dir, fam) for fam in os.listdir(top_dir)]
    j = 0
    for i, f in enumerate(fams):
        if 'txt' not in f:
            results = os.path.join(f, 'Results')
            results_files = os.listdir(results)
            for r in results_files:
                if r[-3:] == 'csv' and 'matches' not in r:
                    csv_files.append(os.path.join(results, r))
    return csv_files


def family_sums(remap_files):
    fam_dict = {}
    for i, remap_file in enumerate(remap_files):
        with open(remap_file) as f:
            reader = csv.reader(f)
            fam = os.path.basename(remap_file)[:-4]
            fam_dict[fam.replace('_', '').upper()] = [0, 0]  # solo, intact
            for row in reader:
                if row[5] == 'S':
                    fam_dict[fam.replace('_', '').upper()][0] += 1
                elif row[5] == 'I':
                    fam_dict[fam.replace('_', '').upper()][1] += 1
    return fam_dict


def diff_dict(remap_dict, old_dict):
    diff_dict = {}
    for fam in old_dict:
        try:
            s, i = remap_dict[fam][0] - old_dict[fam][0], remap_dict[fam][1] - old_dict[fam][1]
            diff_dict[fam] = [s, i]
        except KeyError as e:
            continue
    return diff_dict


old_dict = old_el_reader(SUMM)
old_dict_2 = {}
remap_dict = family_sums(get_remap_files(TOP))

for f in old_dict:
    if f in remap_dict:
        old_dict_2[f] = old_dict[f]

def bin_items(l, n=50):
    l, c, t, f, b = sorted(l, reverse=True), 0, n, [0], []
    count = l.pop()
    while l:
        if count < t:
            f[c] += 1
            count = l.pop()
            b.append(t)
        else:
            t += n
            c += 1
            f.append(0)

    return f, b

def make_labels(b1, b2, ):
    labels = []
    i = 0
    for i in b:
        labels.append('{} - {}'.format(i - n, i))

    print(labels)
    return labels

import numpy as np


def set_breaks(l1, l2, n=50):
    b1 = np.bincount(l1)
    b2 = np.bincount(l2)
    print(sorted(l1))


sum_list_1 = [sum(l) for l in remap_dict.values()]
sum_list_2 = [sum(l) for l in old_dict_2.values()]



def element_distrabutions(dict_1, dict_2):
    sum_1 = [sum(l) for l in dict_1.values()]
    sum_2 = [sum(l) for l in dict_2.values()]
    fig, (ax1, ax2) = plt.subplots(1, 2, num='seaborn')

    fig.suptitle('Gypsy Superfamily Transposable Element Distrabution')

    ax1.hist(sum_1, color='skyblue')
    ax2.hist(sum_2, color='tan')
    ax1.set_yscale('log')
    ax1.text(1500, 50, 'Total Elements\n{}'.format(sum(sum_1)))
    ax2.text(1600, 50, 'Total Elements\n{}'.format(sum(sum_2)))
    ax2.set_yscale('log')
    ax1.set_title('W82 2.1 (Remapped)')
    ax2.set_title('W82 1.1')

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax1.set_xlabel('Number of Elements')
    ax2.set_xlabel('Number of Elements')

    ax1.set_ylabel('Number of Families')
    ax2.set_ylabel('Number of Families')

    plt.yscale('log')
    plt.show()

element_distrabutions(remap_dict, old_dict_2)

def dist_by_num_elements(dict_1, dict_2, col1='skyblue', col2='tan', solo_index=0,
                         intact_index=1, n=50):
    sum_list_1 = [sum(l) for l in dict_1.values()]
    sum_list_2 = [sum(l) for l in dict_2.values()]
    bin_list_1, breaks_1 = bin_items(sum_list_1, n)  # n = 50
    bin_list_2, breaks_2 = bin_items(sum_list_2, n)
    y_pos_1 = np.arange(len(bin_list_1))
    y_pos_2 = np.arange(len(bin_list_2))

    plt.figure()
    plt.subplot(211)
    plt.title('Williams 82 2.1 (Remapped)')
    plt.xlabel('Number of Elements')
    plt.ylabel('Number of Families')
    plt.bar(y_pos_1, bin_list_1, color=col1)
    plt.xticks(y_pos_1, make_labels(breaks_1), rotation='vertical', fontsize=5)
    plt.yscale('log')

    plt.subplot(212)
    plt.title('Williams 82 1.1 (Current Library)')
    plt.xlabel('Number of Elements')
    plt.ylabel('Number of Families')
    plt.bar(y_pos_2, bin_list_2, color=col2)
    plt.xticks(y_pos_1, make_labels(breaks_2), rotation='vertical', fontsize=5)
    plt.yscale('log')
    plt.subplots_adjust(hspace=0.8)




    plt.show()

import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
#b = dist_by_num_elements(remap_dict, old_dict_2)
#print(max(dist_by_num_elements(old_dict_2)))

def linear_compare(dict_1, dict_2):
    x, y = [], []
    x1, y1 = [], []
    for k in dict_1:
        x.append(dict_1[k][0])
        y.append(dict_2[k][0])
        x1.append(dict_1[k][1])
        y1.append(dict_2[k][1])
    #plt.scatter(x, y)
    print(x)
    print(y)
    fig, ax = plt.subplots()

    ax.scatter(x, y, alpha=0.8, color='skyblue')
    ax.scatter(x1, y1, color='orange', alpha=0.8)
    ax.plot(np.arange(0, 5000), color='grey')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    transform = ax.transAxes
    plt.title('Solo and Intact Element Differences by Family')
    plt.ylabel('Remapped Family Total')
    plt.xlabel('Original Family Total')
    solo_patch = mpatches.Patch(color='orange', label='Solo')
    intact_patch = mpatches.Patch(color='skyblue', label='Intact')
    plt.legend(handles=[solo_patch, intact_patch], loc='bottom right')
    plt.show()


linear_compare(remap_dict, old_dict_2)



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
    plt.xticks(ind, groups, fontsize=4, rotation=90)
    p_solo = plt.bar(ind, solo_values, width, color=solo_col)
    p_intact = plt.bar(ind, intact_values, width, bottom=solo_values, color=intact_col)
    plt.ylabel('Total Elements')
    plt.title('Remapped Element Totals by Family')
    plt.show()

#stack_bars(diff_dict(old_dict_2, remap_dict))
