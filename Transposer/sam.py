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


class Sam():

    def __init__(self, path, bdb, accession_path):
        self.path = path
        self.name = os.path.basename(path)
        self.element_set = set([])
        self.dbd = Blast_DB(bdb, accession_path)

    def make_element_set(self):
        try:
            element_set = set([])
            with open(self.path) as path:
                reader = csv.reader(path, delimiter='\t')
                for row in reader:
                    if row[0][0] != '@':
                        name = row[0]
                        acc = process_rname(row[2])
                        start = int(row[3])  # 1 based start
                        length = cigarParser(row[5])
                        end = start + length
                        sequence = self.dbd.search(start, end, acc)
                        element_set.add(
                            Element(name, acc, start, end, length, None, sequence))
            self.element_set = element_set
            return 0
        except FileNotFoundError as e:
            return set([])
