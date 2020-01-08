import csv
import os
from data_gathering import *


CHRS = [56831395, 48577505, 45779434, 52388848, 42234498, 51416163, 44630203,
        47837022, 50189283, 51566898, 34766867, 40091314, 45874162, 49042192,
        51756343, 37887014, 41641078, 58018742, 50746504, 47900577]



TOP = '/media/ethan/EH_DATA/TARP_Runs/Big_Run_12-30_AMB'
WRITE = '/media/ethan/EH_DATA/PAG/remapped_ideo.txt'



def extract_element(row):
    # extracts graphing info from element csv entry
    chr = int(row[2])
    start = int(row[3])
    end = int(start) + int(row[4])
    type = row[5]

    return chr, start, end, type




def write_row(chr, start, end, name, stain):
    return 'chr{}\t{}\t{}\t{}\t{}\n'.format(chr, start, end, name, stain)

#elements = format_sorted_elements((get_old_elements()))
#elements = format_sorted_elements(sort_new_elements(get_elements(get_result_csvs())))



def ideogram(sorted_elements, output):
    scale = 0  # change to increase element length to show up better on graph
    stain_dict = {0: 'gneg', 'S': 'gpos25', 'I':'gpos50'}
    with open(output, 'w') as out:
        out.write('chrom\tstart\tend\tname\tgieStain\n')
        i = 0
        for cc, chr in enumerate(sorted_elements):
            print(len(chr), 'l')
            cp =  0
            while chr:
                i+=1
                ce = chr.pop(0)  # pull element from chr list
                #if ce[1] > cp:
                out.write(write_row(cc+1, cp, ce[1], 'CHR', stain_dict[0]))
                cp = ce[1]  # set current position to start of element
                out.write(write_row(cc+1, ce[1], ce[2]+scale, ce[3], stain_dict[ce[3]]))  # color element based on type
                cp = ce[2] + scale  # set current position to end of element
                for el in chr:
                    el = (el[0], el[1] + scale, el[2] + scale, el[3])
                CHRS[cc] += scale
                #else:
            # at this point written all elements
            # need to write end of chromsome
            out.write(write_row(cc+1, cp, CHRS[cc], 'CHR', stain_dict[0]))


#ideogram(elements, WRITE)















'''







        while sorted_elements:
            ce = extract_element(sorted_elements.pop(0))
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
                '''
