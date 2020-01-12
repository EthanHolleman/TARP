#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
from Transposer.sam_utils import *
from Transposer.element import Element
from fasta_tools import check_formating



class Search():


    def __init__(self, BTI, con_file, out_file, num_old_els, acc, BDB, intact_len, type="I"):
        self.BTI = BTI
        self.BDB = BDB
        self.acc = acc
        self.acc_dict = make_acc_dict(acc)
        self.con_file = con_file
        self.intact_len = intact_len
        self.out_file = out_file  # path to a specific file
        self.num_old_els = num_old_els
        self.type = type  # from a solo or intact consensus
        self.element_set = None
        # old elments must be in fasta for this to be correct 1 seq per line

    def remove_dups(self):
        '''
        Removes elements from the element set that are duplicates within an individual
        sam file. Additional processing required at the run level to seperate
        elements that are duplicate between sam files.
        '''
        id_set = set([])
        nr_list = set([])
        for e in self.element_set:
            id = str(e.chr) + str(e.startLocation)
            if id not in id_set:
                nr_list.add(e)
                id_set.add(id)

        self.element_set = nr_list

    def type_elements(self, allowance=25):
        '''
        If Sam is solo type then sort solo elements and remove those that are
        actually the ends of an LTR. If sam is an intact type then change the
        status of all elements in the set to intact.
        '''
        if self.type == 'S':
            solo_set = set([])
            sort_e = sort_elements(self.element_set)
            # sort solos by start location
            i = 0
            n = len(sort_e) - 1
            if n == 0:  # if only one element still gets added
                sort_e[n].status = 'S'
                solo_set.add(sort_e[n])
            else:
                while i < n:  # if only one element fails
                    current = sort_e[i]
                    next = sort_e[i + 1]  # watch out of bounds error
                    if current.chr != next.chr:
                        solo_set.add(current)
                        i += 1
                    elif next.startLocation - self.intact_len + allowance <= current.endLocation:
                        # elements are close enough to be intact
                        i += 2
                    else:
                        current.status = 'S'
                        solo_set.add(current)
                        i += 1
                    if i == n:  # last element if only one must be a solo
                        sort_e[n].status = 'S'
                        solo_set.add(sort_e[n])
                        break
            self.element_set = solo_set
        else:
            for e in self.element_set:
                e.status = 'I'

    # need to put get seq and flanks in a more general position so that
    # the old elements can use it

    def make_element_set(self):
        '''
        Takes in a sam file from the bowtie2 search and reads the allignments
        into a list of Element objects. Removes duplicate elements.
        '''
        flank_len = 20
        try:
            element_set = set([])
            with open(self.out_file) as path:
                reader = csv.reader(path, delimiter='\t')
                for i,row in enumerate(reader):
                    if row[0][0] != '@':
                        name = row[0]
                        acc = row[2]
                        if row[2] == '*':
                            continue
                        chr = self.acc_dict[str(acc)]
                        start = int(row[3])  # 1 based start
                        length = cigarParser(row[5])
                        end = start + length
                        seq, left, right = get_seq_flanks(start, end, acc, self.BDB)
                        # get flanking sequences for backmapping to avoid additional searches
                        element_set.add(
                            Element(name, acc, chr, start, end, length, self.type, seq, left, right))

            self.element_set = element_set
            self.remove_dups()  # remove duplicate elements
            self.type_elements()  # identify false solos (actually LTRs)

            return 0
        except FileNotFoundError as e:
            return set([])


    def search_BTI(self, defualt=True, custom=None, k_fuct=10,
                   preset='--sensitive', threads=8):
        '''
        Searches the BTI index using a consensus sequence. By defualt will use
        the sensitive setting for the search and look for a max of 10 times
        the number of old elements linked to the searched consensus sequence.
        '''
        defualt_cmd = ['bowtie2', '-x', self.BTI, '-f', self.con_file, '-k',
                       round(self.num_old_els * k_fuct + self.num_old_els), preset, '-S', self.out_file, '--n-ceil', 'L,0,0.20.']
        cmd = [str(c) for c in defualt_cmd]
        try:
            FNULL = open(os.devnull, 'w')
            retcode = subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
            if retcode != 0:
                print('Bowtie2 exited with code {}'.format(retcode))
                sys.exit(retcode)
            self.make_element_set()

        except subprocess.CalledProcessError as e:
            return e
