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


def write_fasta(sorted_elements, output):
    with open(output, 'w') as out:
        for el in sorted_elements:
            out.write(el.get_header() + '\n' + el.seq + '\n')


def write_csv(sorted_elements, output):
    with open(output, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(['Name', 'Accession', 'Chr', 'Start', 'Length',
                         'Seq', 'Left Flank', 'Right Flank'])
        for el in sorted_elements:
            writer.writerow(el.get_row())


def rename_elements(sorted_elements):
    sorted_elements = list(sorted_elements)  # deal with generator stuff for now
    i = 1
    cur_chr = None
    for el in sorted_elements:
        if cur_chr != el.chr:
            cur_chr = el.chr
            i = 1
        el.name = '{}:{}-{}'.format(el.name, cur_chr, i)
        i += 1

    return sorted_elements
