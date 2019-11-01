#!/usr/bin/env python3
import csv
import os
import itertools
from operator import or_
from functools import reduce

from Transposer.search import sort_elements

# move this to the run object class since each sam is saved in the jobs attribute


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

def prune(e, n=75):
    '''
    Takes a sorted list of elements and removes elements that are within
    n number of bases (start position) from the last element. This prevents
    similar consensus sequences from hitting elements at only few base
    pair positional difference from making it into the final output.
    '''
    chunks = []
    cur_chunk = []
    for i in range(1, len(e)):
        d = e[i].startLocation - e[i-1].startLocation
        if d > n:  # no nearby elements
            if len(cur_chunk) > 0:
                if len(cur_chunk) >= 1:
                    yield cur_chunk[0]
                cur_chunk = []
            yield e[i]
        else:
            cur_chunk.append(e[i])

def pruner(seq, n=75):
    cur_chunk = []
    i = 0
    for i in range(0, len(seq)):
        if i == len(seq) -1:
            break
        d = seq[i+1].startLocation - seq[i].startLocation
        if d > n:
            yield seq[i]


def write_results(sels, path):
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
    writer.writerow(['Name', 'Accession', 'Chr', 'Start', 'Length',
                     'Status', 'Seq', 'Left Flank', 'Right Flank'])


def rename_elements(sels): # deal with generator stuff for now
    i = 1
    cur_chr = None
    for el in sels:
        if cur_chr != el.chr:
            cur_chr = el.chr
            i = 1
        el.name = '{}:{}-{}'.format(el.name, cur_chr, i)
        i += 1
        yield el
