import subprocess
import os

from Transposer.blast_BD import Blast_DB
from Transposer.sam import Sam


class Search():

    def __init__(self, BTI, con_file, out_file, num_old_els):
        self.BTI = BTI
        self.sam = None
        self.con_file = con_file
        self.out_file = out_file  # path to a specific file
        self.num_old_els = num_old_els
        # old elments must be in fasta for this to be correct 1 seq per line

    def search_BTI(self, bdb, acc_path, defualt=True, custom=None, k_fuct=1.20, preset='--very-fast'):
        defualt_cmd = ['bowtie2', '-x', self.BTI, '-f', self.con_file, '-k',
                       round(self.num_old_els * k_fuct), preset, '-S', self.out_file]
        try:
            subprocess.call(defualt_cmd, shell=True)
            self.sam = Sam(path=self.out_file, bdb=bdb,
                           accession_path=acc_path)
            return 0
        except subprocess.CalledProcessError as e:
            return e
