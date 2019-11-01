#!/usr/bin/env python3
import subprocess
import os
import csv

from Transposer.element import Element
from fasta_tools import check_formating

# each search object needs to include or need somekind of way for it
# to access the intact consensus file when the solo elements are
# being processed

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
    return r_name.split('|')[3]


def get_chr_num(acc_path, acc):
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
    acc_dict = {}
    with open(acc_path) as names:
        for line in names:
            l = line.strip().split('\t')
            acc_dict[l[1]] = l[0]
    return acc_dict


def sort_elements(elements):
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


def search_BDB(BDB, start, end, entry, r_seq=True):
    '''
    search the blast db for a seq in an entry
    '''
    cmd = ['blastdbcmd', '-db', BDB, '-dbtype', 'nucl', '-range', str(start) + '-' + str(end), '-entry', entry]
    cmd = ' '.join(cmd)
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        return ''
    if r_seq:
        return ''.join(str(output).split('\\n')[1:])
        # returns string of just the seq
    else:
        return output  # else return all output


def get_seq_flanks(start, end, entry, BDB, flank_len):
    '''
    Uses search BDB to get the seq and right and left flanks.
    Returns as a tuple (seq, left, right)
    '''
    seq = self.search_BDB(BDB, start-flank_len, end+flank_len, entry)
    left, right = seq[:flank_len], seq[len(seq)-flank_len:]
    seq = seq[flank_len:-flank_len]
    return tuple([seq, left, right])


class Search():
    '''
    # TODO: going to need to pass the BDB to the search objects now
    becuase it is creating the element sets this also means that
    going to need the the accession files
    '''

    def __init__(self, BTI, con_file, out_file, num_old_els, acc, BDB, intact_len, type="I"):
        self.BTI = BTI
        self.BDB = BDB
        self.acc = acc
        self.acc_dict = make_acc_dict(acc)
        self.con_file = con_file
        self.intact_len = intact_len
        # Con file length needs to be the length of the
        # intact element always even if the type is solo
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
            while i < n:
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

    def search_BDB(self, start, end, entry, r_seq=True):
        '''
        search the blast db for a seq in an entry
        '''
        cmd = ['blastdbcmd', '-db', self.BDB, '-dbtype', 'nucl', '-range', str(start) + '-' + str(end), '-entry', entry]
        cmd = ' '.join(cmd)
        try:
            output = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            return ''
        if r_seq:
            return ''.join(str(output).split('\\n')[1:])
            # returns string of just the seq
        else:
            return output  # else return all output

    # need to put get seq and flanks in a more general position so that
    # the old elements can use it

    def make_element_set(self):
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
                        seq = self.search_BDB(start-flank_len, end+flank_len, acc)
                        left, right = seq[:flank_len], seq[len(seq)-flank_len:-1]  # ' at end of sequence
                        seq = seq[flank_len:-flank_len]
                        # get flanking sequences for backmapping to avoid additional searches
                        element_set.add(
                            Element(name, acc, chr, start, end, length, self.type, seq, left, right))

            self.element_set = element_set
            self.remove_dups()
            self.type_elements()

            return 0
        except FileNotFoundError as e:
            return set([])


    def search_BTI(self, defualt=True, custom=None, k_fuct=3, preset='--sensitive', threads=8):
        defualt_cmd = ['bowtie2', '-x', self.BTI, '-f', self.con_file, '-k',
                       round(self.num_old_els * k_fuct), preset, '-S', self.out_file, '--n-ceil', 'L,0,0.20.']
        cmd = [str(c) for c in defualt_cmd]
        try:
            FNULL = open(os.devnull, 'w')
            retcode = subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
            self.make_element_set()

        except subprocess.CalledProcessError as e:
            return e


# moving functionalty away from the sam file object its not really doing anything unique
