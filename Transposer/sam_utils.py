import csv
import os
import itertools
import subprocess

from operator import or_
from functools import reduce
from collections import deque



def sort_elements(elements):
    '''
    Takes in a list of elements, sorts by chromosome and then by start location.
    '''
    chr_dict = {}
    sorted_list = []
    for e in elements:
        if e.chr in chr_dict:
            chr_dict[e.chr].append(e)
        else:
            chr_dict[e.chr] = [e]

    for chr in sorted(chr_dict):
        sorted_list += sorted(chr_dict[chr], key=lambda e: (e.startLocation))

    return sorted_list


def sort_sams(search_list):
    '''
    Takes in run object and returns sorted list of the solo and intact element
    for that object. Sort by chromosomes and sort by position in the
    chromosomes.
    Each sam can sort itself but want to make sure you do all processing first
    before anyting else happens
    This would be called in the run and given a list of all the search jobs after
    they have been run. Then this function can remove duplicates in each sam
    and the duplicates between sams the sort and return a sorted list of
    elements which then can be written to a fasta file.
    '''
    els_list = [s.element_set for s in search_list]
    # get all element sets from search objects
    all_els = list(itertools.chain.from_iterable(els_list))
    # chain all lists together
    sort_els = sort_elements(all_els)
    # sort the elements as if one large list

    return sort_els


def pruner(seq, f=75):
    '''
    Removes elements by their start location. If start location of two elements
    is within f distance they are assumed to be duplicates and one is removed.
    If multible elements form a "chain" where each is within f distance of the
    last all these elements will be removed except the first in the chain.
    '''
    seq = deque(seq)
    run = [seq.popleft()]
    while seq:
        n = seq.popleft()
        if run[-1].chr != n.chr:
            yield run[0]
            run = [n]
        elif n.startLocation - run[-1].startLocation < f:
            run.append(n)
        else:
            yield run[0]
            run = [n]
    if run:
        yield run[0]


def write_results(sels, path):
    '''
    Writes both fasta and csv results of completed TARP run.
    '''
    fasta = path + '.fa'
    csv_file = path + '.csv'
    fasta, csv_file = open(fasta, 'w'), open(csv_file, 'w')

    writer = csv.writer(csv_file)  # make csv writer
    write_csv_header(writer)  # write header row

    for el in sels:
        t = el
        fasta.write(t.get_header() + '\n' + t.seq + '\n')
        writer.writerow(t.get_row())


def write_csv_header(writer):
    '''
    Writes Header for backmap csv results. Takes in a csv writer object.
    '''
    writer.writerow(['Name', 'Accession', 'Chr', 'Start', 'Length',
                     'Status', 'Seq', 'Left Flank', 'Right Flank'])


def rename_elements(sels):  # deal with generator stuff for now
    '''
    Renames elements based on their relative chromosomal location to other
    elements on thar chromosome. Format is [el name]:[chr number]-[order]
    '''
    cur_chr, i = None, 1
    for el in sels:
        if cur_chr != el.chr:
            cur_chr = el.chr
            i = 1  # reset i at each new chromosome
        el.name = '{}-{}'.format(cur_chr, i)
        i += 1
        yield el


def cigarParser(cigar):
    '''
    reads the cigar info from an Bowtie allignment and translates
    into the length on the alligned reference
    '''
    length = 0
    temp = ""
    for char in cigar:
        temp = temp + "" + char
        if char == "M" or char == "D":
            temp = temp[:-1]
            length += int(temp)
            temp = ""
        elif char == "I" or char == "H":
            temp = ""

    return length


def process_rname(r_name):
    '''
    formats r_name from same file.
    '''
    return r_name.split('|')[3]


def get_chr_num(acc_path, acc):
    '''
    Given an acc_2_chr file path and a specific acc, returns the chr number
    linked to that accession in the given acc_2_chr file.
    '''
    try:
        lines = []
        with open(acc_path) as names:
            for line in names:
                l = line.strip().split('\t')
            if l[1] == acc:
                return int(l[0])
    except FileNotFoundError as e:
        return -1


def make_acc_dict(acc_path):
    '''
    Given an chr_2_acc type file returns a dictionary with acc number mapped
    to chromosome numbers.
    '''
    acc_dict = {}
    with open(acc_path) as names:
        for line in names:
            l = line.strip().split('\t')
            acc_dict[l[1]] = l[0]
    return acc_dict


def search_BDB(BDB, start, end, entry, r_seq=True):
    '''
    search the blast db for a seq in an entry
    '''
    cmd = ['blastdbcmd', '-db', BDB, '-dbtype', 'nucl',
           '-range', str(start) + '-' + str(end), '-entry', entry]
    cmd = ' '.join(cmd)
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        return ''
    if r_seq:
        return ''.join(str(output).split('\\n')[1:-1])
        # returns string of just the seq
    else:
        return output  # else return all output


def get_seq_flanks(start, end, entry, BDB, flank_len=20):
    '''
    Uses search BDB to get the seq and right and left flanks.
    Returns as a tuple (seq, left, right)
    '''
    if entry == None:
        return ('', '', '')
    else:
        seq = search_BDB(BDB, int(start) - flank_len,
                         int(end) + flank_len, entry)
        left, right = seq[:flank_len], seq[len(seq) - flank_len:]
        seq = seq[flank_len:-flank_len]
        return (seq, left, right)
