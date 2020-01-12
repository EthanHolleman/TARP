import os
import csv
TOP = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_12-30_AMB'
OLD = '/media/ethan/EH_DATA/Element_Summaries/IILTRGypsy.txt'

def get_old_elements():
    elements = []
    with open(OLD) as old:
        reader = csv.reader(old, delimiter='\t')
        next(reader)
        for row in reader:
            if len(row) > 11:
                e = extract_old_element(row)
                if e:
                    elements.append(e)
    return elements

def extract_old_element(old_el):
    chr, type, start, end,= old_el[8], old_el[7], int(old_el[10]), int(old_el[11])
    if 'Gm' in chr:
        chr = int(chr.split('Gm')[-1].lstrip('0'))
    else:
        return None
    if type == 'SOLO':
        type = 'S'
    elif type == 'INTACT':
        type = 'I'
    return 0, 0, chr, start, end-start, type

def get_result_csvs():
    csvs = []
    families = [os.path.join(TOP, fam) for fam in os.listdir(TOP)]
    results = [os.path.join(fam, 'Results') for fam in families]
    for r in results:
        result_files = os.listdir(r)
        for f in result_files:
            if f[-3:] == 'csv' and 'matches' not in f:
                    csvs.append(os.path.join(r, f))
    return csvs

def get_elements(csvs, delimiter=','):
    elements = []
    for c in csvs:
        with open(c) as C:
            reader = csv.reader(C, delimiter=delimiter)
            next(reader)
            for row in reader:
                elements.append(row[:-1])
    return elements

def extract_element(row):
    # extracts graphing info from element csv entry
    chr = int(row[2])
    start = int(row[3])
    end = int(start) + int(row[4])
    type = row[5]

    return chr, start, end, type

def sort_elements(element_list, chrs=20, chr_index=2, start_index=3):
    bins = [[] for j in range(0, chrs)]
    for el in element_list:
        print(int(el[chr_index])-1)
        bins[int(el[chr_index])-1].append(el)
    for i, b in enumerate(bins):
        t = sorted(b, key=lambda x: x[int(start_index)])
        b[i] = t

    return bins

def sort_new_elements(element_list, chrs=20):
    bins = [[] for j in range(0, chrs)]
    new_bins = []
    for el in element_list:
        bins[int(el[2])-1].append(el)
    for i, b in enumerate(bins):
        t = sorted(b, key=lambda x: int(x[3]))
        new_bins.append(t)
    return new_bins

def format_sorted_elements(sorted_elements):
    for i in range(0, len(sorted_elements)):
        for j in range(0, len(sorted_elements[i])):
            sorted_elements[i][j] = extract_element(sorted_elements[i][j])

    return sorted_elements
