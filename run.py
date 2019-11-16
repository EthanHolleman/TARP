#!/usr/bin/env python3
import os
import subprocess

from progress.bar import Bar
from progress.bar import ChargingBar

from Transposer.search import Search
from par.multi_process import run_pool
from Clustering.ClstrFile import ClstrFile
from Clustering.big_hit import run_cd_hit


from fasta_tools import make_consensus
from fasta_tools import check_formating
# return the length of the sequence which will
# be in the second line


def get_intact_length(con_file_path):
    check_formating(con_file_path)
    length = []
    with open(con_file_path) as con:
        length = con.readlines()

    return len(length[1])


def make_clstr(fasta, clstr_dir):
    clstr_file_name = make_clstr_name(fasta, clstr_dir)
    print('Clustering', fasta)
    run_cd_hit(clstr_file_name, fasta)
    return clstr_file_name + '.clstr'


def make_clstr_name(path, clstr_dir):
    base = os.path.basename(path).split('.')[0]
    return os.path.join(clstr_dir, base)


def make_dirs(output, run_name):
    '''
    Create a run dir to store all files and then three dirs within it for
    cluster files, sam files and result files. Set the write_dirs variable
    to a list of the three low level dirs.
    '''
    run_dir = os.path.join(output, run_name)
    if os.path.isdir(run_dir) is False:
        os.mkdir(run_dir)
    clstr_dir = os.path.join(run_dir, 'Clusters')
    sam_dir = os.path.join(run_dir, 'Sam_Files')
    data_dir = os.path.join(run_dir, 'Results')
    write_dirs = [clstr_dir, sam_dir, data_dir]
    for dir in write_dirs:
        if os.path.isdir(dir) is False:
            os.mkdir(dir)

    return tuple(write_dirs)


class Run():

    def __init__(self, cur_BDB, cur_acc, old_BDB, old_acc, BTI, cie, run_name, output, csi=None, min_els=2):
        self.write_dirs = make_dirs(output, run_name)
        self.min_els = min_els
        self.old_BDB = old_BDB  # old_acc
        self.cur_BDB = cur_BDB  # cur_acc
        self.BTI = BTI  # path to bowtie index
        self.cie = cie  # path to file
        self.csi = csi  # path to file
        self.old_acc = old_acc  # path to acc2chr file for old assembly
        self.cur_acc = cur_acc  # path to acc2chr file for new assembly
        self.cie_clstrs = ClstrFile(
            path=make_clstr(self.cie, self.write_dirs[0]))
        self.num_cie = num_cie = sum(1 for line in open(self.cie)) / 2
        if csi is not None:
            self.csi_clstrs = ClstrFile(
                path=make_clstr(self.csi, self.write_dirs[0]))
            self.num_csi = num_cie = sum(1 for line in open(self.csi)) / 2
        else:
            self.csi_clstrs = None
            self.num_csi = 0
        self.cie_fastas = None
        self.csi_fastas = None
        self.cie_cons = None
        self.csi_cons = None
        self.jobs = None

    def select_clusters(self, min_elements=3):
        '''
        Remove clusters that do not meet the min element
        threshold. May want to calculate this differently later
        or make so you can pass in a function that allows for
        calculating min value.
        '''
        self.cie_clstrs.trim_clusters(min_elements)
        if self.csi_clstrs is not None:
            self.csi_clstrs.trim_clusters(min_elements)

    def make_clstr_fastas(self):
        '''
        Write fasta files from the clusters. Should trim the clusters before
        calling this method pretty much always.
        '''
        self.cie_cons = self.cie_clstrs.write_cluster_fastas(
            self.cie, self.write_dirs[0], new_dir=False)
        if self.csi_clstrs is not None:
            self.csi_cons = self.csi_clstrs.write_cluster_fastas(
                self.csi, self.write_dirs[0], new_dir=False)

    def make_consensensi_teo(self, min_elements=2, n=21):
        bar = Bar('Making intact consensuses', max=len(self.cie_cons))
        for clstr in self.cie_clstrs.clusters_set:
            con_name = clstr.fasta + '_intact_con'
            bar.next()
            make_consensus(clstr.fasta, con_name, min_elements, n=n)
            clstr.consensus = con_name

        bar_2 = Bar('Making solos consensuses',
                    max=len(self.csi_cons))
        if self.csi_clstrs is not None:
            for clstr in self.csi_clstrs.clusters_set:
                con_name = clstr.fasta + '_solo_con'
                bar_2.next()
                make_consensus(clstr.fasta, con_name, min_elements, n=n)
                clstr.consensus = con_name

    def make_consensensi(self, min_elements=2, n=21):
        '''
        Need a way to get the paths to the fasta files and then process
        them with clustal omega and that kind of thing. Need to store the
        paths to consensus sequence in the self.con variables in run.
        '''
        # reduce output here
        cie_cons, csi_cons = [], []
        bar = Bar('Making intact consensuses', max=len(self.cie_cons))
        for fasta in self.cie_cons:
            bar.next()
            con_name = fasta + '_intact_con'
            make_consensus(fasta, con_name, min_elements, n=n)
            cie_cons.append(con_name)
        self.cie_cons = cie_cons

        if self.csi_cons is not None:
            print('\n')
            bar_2 = ChargingBar('Making solos consensuses',
                                max=len(self.cie_cons))
            for fasta in self.csi_cons:
                bar_2.next()
                con_name = fasta + 'solo_con'
                make_consensus(fasta, con_name, min_elements, n=n)
                csi_cons.append(con_name)
            self.csi_cons = csi_cons

    def make_jobs_two(self):
        sam_dir = self.write_dirs[1]
        jobs = []
        ave_len = 0
        for clstr in self.cie_clstrs.clusters_set:
            sam_name = 'clstr_' + clstr.num + '_intact.sam'
            sam_file = os.path.join(sam_dir, sam_name)
            l = get_intact_length(clstr.consensus)
            jobs.append(Search(BTI=self.BTI, con_file=clstr.consensus,
                               out_file=sam_file, num_old_els=clstr.num_elements,
                               type='I', acc=self.cur_acc, BDB=self.cur_BDB,
                               intact_len=l))
            ave_len += l

        for clstr in self.csi_clstrs.clusters_set:
            sam_name = 'clstr_' + clstr.num + '_solo.sam'
            sam_file = os.path.join(sam_dir, sam_name)
            jobs.append(Search(BTI=self.BTI, con_file=clstr.consensus,
                               out_file=sam_file, num_old_els=clstr.num_elements,
                               type='S', acc=self.cur_acc, BDB=self.cur_BDB,
                               intact_len=ave_len))
        self.jobs = jobs

    def run_jobs(self, threads=1):
        '''
        Runs the jobs and stores sam file objects. using methods from sam class
        removes duplicates and types the elements for each sam object. After this
        the elements are ready to be placed into a final order and then written to
        a fasta file.
        '''
        print('\nRunning {} bowtie2 jobs'.format(len(self.jobs)))
        sam_dir = self.write_dirs[1]  # stored at 1 index always
        for job in self.jobs:
            job.search_BTI()

    def write_meta(self, run_name):
        '''
        Writes to output dir metadata about the TARP run. Will include info
        such as number of solo and intact clusters created, number of jobs run,
        total run time, the command used, output directory paths etc.
        '''
        pass
