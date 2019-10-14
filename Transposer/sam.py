import csv
import os
import subprocess
from Transposer.blast_BD import Blast_DB
from Transposer.element import Element
# want this to work with eleemnt objects


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

def sort_elements(elements):
    chr_dict = {}
    sorted_list = []
    for e in elements:
        if e.chr in chr_dict:
            chr_dict[e.chr].append(e)
        else:
            chr_dict[e.chr] = [e]

    for chr in sorted(chr_dict):
        sorted_list += sorted(chr_dict[chr], key= lambda e: (e.startLocation))

    return sorted_list


class Sam():
    def __init__(self, path, bdb, accession_path, type, element_set=set([])):
        self.path = path
        self.acc_path = accession_path
        self.name = os.path.basename(path)
        self.element_set = element_set
        self.dbd = Blast_DB(bdb, accession_path)
        self.type = type

    def make_element_set(self):
        try:
            element_set = set([])
            with open(self.path) as path:
                reader = csv.reader(path, delimiter='\t')
                for row in reader:
                    if row[0][0] != '@':
                        name = row[0]
                        acc = process_rname(row[2])
                        chr = get_chr_num(self.acc_path, acc)
                        start = int(row[3])  # 1 based start
                        length = cigarParser(row[5])
                        end = start + length
                        sequence = self.dbd.search(start, end, acc)
                        element_set.add(
                            Element(name, acc, chr, start, end, length, type, sequence))
            self.element_set = element_set
            return 0
        except FileNotFoundError as e:
            return set([])

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


    def type_elements(self, len_intact, allowance=25):
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
            n = len(sort_e)-1
            print(n, 'number of element')
            while i < n:
                print(i)
                current = sort_e[i]
                next = sort_e[i+1]  # watch out of bounds error
                if current.chr != next.chr:
                    solo_set.add(current)
                    i+=1
                elif next.startLocation - len_intact - allowance <= current.endLocation:
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
