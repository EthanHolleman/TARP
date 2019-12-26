import csv
import os


CHRS = [56831395, 48577505, 45779434, 52388848, 42234498, 51416163, 44630203,
        47837022, 50189283, 51566898, 34766867, 40091314, 45874162, 49042192,
        51756343, 37887014, 41641078, 58018742, 50746504, 47900577]



TOP = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_3'
WRITE = '/media/ethan/EH_DATA/PAG/test.txt'

def read_results(top_dir):
    fams = [os.path.join(top_dir, fam) for fam in os.listdir(top_dir)]
    csv_results = []
    for f in fams:
        results = os.path.join(f, 'Results')
        results_files = os.listdir(results)
        for r in results_files:
            if r[-3:] == 'csv':
                csv_results.append(os.path.join(results,r))
    return csv_results


def extract_element(row):
    # extracts graphing info from element csv entry
    chr = row[2]
    start = row[3]
    end = int(start) + int(row[4])
    type = row[5]

    return chr, start, end, type

def write_row(chr, start, end, name, stain):
    return 'chr{}\t{}\t{}\t{}\t{}\n'.format(chr, start, end, name, stain)



def bin_elements(elements, num_chr=20):
    bins = [[] for i in range(1, 21)]
    print(bins)
    for el in elements:
        chr = int(el[0])-1
        bins[chr].append(el)
    return bins


def sort_bin_elements(bins):
    for i, b in enumerate(bins):
        bins[i] = sorted(b, key=lambda x: x[1])

    return [item for sublist in bins for item in sublist]


def ideogram(csv_results, output):
    scale = 100000  # change to increase element length to show up better on graph
    elements = []
    stain_dict = {0: 'gneg', 'S': 'gpos25', 'I':'gpos50'}
    with open(output, 'w') as out:
        out.write('#chr\tchromStart\tchromEnd\tname\tgieStain\n')
        for csv_file in csv_results:
            with open(csv_file) as cf:
                reader = csv.reader(cf)
                reader.next()
                for row in reader:
                    elements.append(extract_element(row))

        sorted_elements = sort_bin_elements(bin_elements(elements))
        cc, i, le = None, 0, None
        while sorted_elements:
            ce = sorted_elements.pop(0)
            if int(ce[0]) != cc:
                if cc != 0 and cc != None:
                    out.write(write_row(cc, i, CHRS[int(cc) -1], 'NA', stain_dict[0]))
                    # write end of chromosome if not first chr
                cc, i = int(ce[0]), 0
                out.write(write_row(cc, 0, ce[1], 'NA', stain_dict[0]))
                out.write(write_row(cc, ce[1], int(ce[2]) + scale, 'NA', stain_dict[ce[3]]))
                i = int(ce[2]) + scale

                # write start new chr and first element
            else:
                out.write(write_row(cc, i, ce[1], 'NA', stain_dict[0]))
                out.write(write_row(cc, ce[1], ce[2] + scale, 'NA', stain_dict[ce[3]]))
                i = int(ce[2]) + scale

id = ideogram(read_results(TOP), WRITE)
