import os
import subprocess

from progress.bar import Bar
from progress.bar import ChargingBar

from Transposer.search import Search
from Clustering.ClstrFile import ClstrFile
from Clustering.big_hit import run_cd_hit

from Run.run_utils import *

from fasta_tools import make_consensus
from fasta_tools import check_formating
# return the length of the sequence which will
# be in the second line


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

    def __repr__(self):
        return (f'{self.__class__.__name__}('
           f'{self.cie}, {self.csi})')

    def select_clusters(self, min_elements=3):
        '''
        Remove clusters that do not meet the min element
        threshold. May want to calculate this differently later
        or make so you can pass in a function that allows for
        calculating min value.
        '''
        def clstr_count(clstr_set):
            '''
            Returns clstrs that are above min number of elements. If no clusters
            are above min elements returns 2 if no clstr has more than 1 element
            returns 1.
            '''
            sum = 0
            clstrs = []
            for clstr in clstr_set:
                if clstr.num_elements >= min_elements:
                    clstrs.append(clstr)
                else:
                    sum += clstr.num_elements
            if clstrs:
                return clstrs
            elif sum == len(clstr_set):
                return clstr_set
            else:
                return None

        self.cie_clstrs.clusters_set = clstr_count(self.cie_clstrs.clusters_set)
        self.csi_clstrs.clusters_set = clstr_count(self.csi_clstrs.clusters_set)


    def make_clstr_fastas(self):
        '''
        Write fasta files from the clusters. Should trim the clusters before
        calling this method pretty much always. At this point the clusters_set
        should be a tuple with first value being 0-2. 0 = clusters were trimmed
        and contain clusters beyond min element threshold. 1 = all clusters in
        set have 1 element and 2 = no clusters have beyond min value but some
        have more than 1 element.
        '''

        def fasta_decision(clstr_file, og_fasta):
            if clstr_file.clusters_set != None:
                w = self.write_dirs[0]
                return clstr_file.write_cluster_fastas(og_fasta, w)
            else:
                return []

        self.cie_cons = fasta_decision(self.cie_clstrs, self.cie)
        if self.csi_clstrs is not None:
            self.cie_cons = fasta_decision(self.csi_clstrs, self.csi)


    def make_consensensi_teo(self, min_elements=2, n=21):

        def make_con(clstr_file, con_str):
            if clstr_file.clusters_set != []:
                for clstr in clstr_file.clusters_set:
                    con_name = clstr.fasta + con_str
                    if clstr.num_elements > 1:
                        make_consensus(clstr.fasta, con_name, min_elements, n=n)
                        clstr.consensus = con_name
                    else:
                        clstr.consensus = clstr.fasta
        print('Making Intact consensuses')
        make_con(self.cie_clstrs, '_intact_con')
        if self.csi_clstrs is not None:
            print('Making Solo consensuses')
            self.cie_cons = make_con(self.csi_clstrs, '_solo_con')


    def make_jobs_two(self):
        '''
        Creates a list of jobs which are composed of search objects. Solo and
        intact jobs are grouped together under run.jobs but each search obejct
        is labeled with type to distinguish later on when solo searches must be
        processed to remove LTRs that are apart of intact elements.
        '''
        sam_dir = self.write_dirs[1]
        jobs = []
        ave_len = 0
        clstrs, sam_suffic = None, None

        def jobs(el_el_type):
            if el_type == 'I':
                clstrs = self.cie_clstrs.clusters_set
                sam_suffix = '_intact.sam'
            else:
                clstrs = self.csi_clstrs.clusters_set
                sam_suffix = '_solo.sam'

            for clstr in clstrs:
                sam_name = 'clstr_' + clstr.num + sam_suffix
                sam_file = os.path.join(sam_dir, sam_name)
                if el_type == 'I':
                    l = get_intact_length(clstr.consensus)
                    ave_len += l
                else:
                    l = ave_len
                jobs.append(Search(self.BTI, clstr.consensus, sam_file,
                            clstr.num_old_els, el_type, self.cur_acc, self.cur_BDB,
                            intact_len=l))
        jobs('I')  # make intact jobs
        ave_len = ave_len // len(self.cie_clstrs.clusters_set)
        jobs('S')  # make solo jobs

        self.jobs = jobs

    def run_jobs(self, threads=8):
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
        Writes to output dir metadata about the TARPs run. Will include info
        such as number of solo and intact clusters created, number of jobs run,
        total run time, the command used, output directory paths etc.
        '''
        pass
