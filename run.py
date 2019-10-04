import os

from Transposer.sam import Sam
from Transposer.search import Search
from par.multi_process import run_pool
from Clustering.ClstrFile import ClstrFile
from Transposer.blast_BD import Blast_DB


class Run():

    def __init__(self, cur_BDB, cur_acc, old_BDB, old_acc, BTI, cie, run_name, output, csi=None):
        self.old_BDB = Blast_DB(old_BDB, old_acc)
        self.cur_BDB = Blast_DB(cur_BDB, cur_acc)
        self.BTI = BTI  # path to bowtie index
        self.cie = cie  # path to file
        self.csi = csi  # path to file
        self.cie_clstrs = ClstrFile(path=cie)
        self.num_cie = num_cie = sum(1 for line in open(self.cie))
        if csi is not None:
            self.csi_clstrs = ClstrFile(path=csi)
            self.num_csi = num_cie = sum(1 for line in open(self.csi))
        else:
            self.csi_clstrs = None
            self.num_csi = 0
        self.cie_fastas = None
        self.csi_fastas = None
        self.cie_cons = None
        self.csi_cons = None
        self.jobs = None
        self.write_dirs = None

    def make_dirs(self):
        '''
        Create a run dir to store all files and then three dirs within it for
        cluster files, sam files and result files. Set the write_dirs variable
        to a list of the three low level dirs.
        '''
        run_dir = os.path.join(self.output, self.run_name)
        if os.isdir(run_dir) is False:
            os.mkdir(run_dir)
        clstr_dir = os.path.join(run_dir, 'Clusters')
        sam_dir = os.path.join(run_dir, 'Sam_Files')
        data_dir = os.path.join(run_dir, 'Results')
        write_dirs = [clstr_dir, sam_dir, data_dir]
        for dir in write_dirs:
            if os.isdir(dir) is False:
                os.mkdir(dir)

        self.write_dirs = tuple(write_dirs)

    def make_clstr_fastas(self):
        '''
        Make clustered fasta files for any input elements given. By defualt you
        must provide at least the path to cie (current intact elements) fasta
        file. If you also have solo elements you can provide the path for that
        fasta as well and clustering will be calculated for it as well.
        By defualt cluster fastas are stored in the current working directory
        but a new location can be specified in the run object init.
        '''

    def select_clusters(self, min_elements=10, rm=True):
        '''
        Iterates through cluster files and removes those with less than
        min_elements number of elements. By defualt those with not enough
        elements are deleted.
        select the clusters and make consensus sequences of them in new directory
        under the run directory. return the paths to the con self variables
        for intact and solo respectively.
        '''
        space = []
        if self.cie_clstrs is None:
            return -1
        elif self.csi_fastas is None:
            space = [self.cie_fastas]
        else:
            space = [self.cie_fastas, self.csi_fastas]

        for s in space:
            files = [os.path.join(self.clstr)]
            # check each file for lines / 2

    def make_jobs(self):
        '''
        generates a list of bowtie commands to be run based on the number of
        consensus clusters in the run object. Each consensus cluster will be
        bowtied and then needs to tripped of duplicate elements. Solo and intact
        will need to be seperated. Will be a list with first list being intact
        jobs and second list being solo jobs.
        '''
        sam_dir = self.write_dirs[1]
        all_cons = []
        jobs = []
        switch = -1

        if self.csi is not None:  # making sure use the right previous element
            # for each job being created
            switch = len(self.cie) - 1
            all_cons = self.cie_cons + csi_cons
        else:
            all_cons = self.cie_cons

        for i, con in enumerate(all_cons):
            # naming the sam file
            sam_name = ''
            if i == switch:
                num_old_els = self.num_csi
                sam_name = os.path.basename(con).split('.')[0] + '_solo.sam'
            else:
                num_old_els = self.num_cie
                sam_name = os.path.basename(con).split('.')[0] + '_intact.sam'

            sam_file = os.path.join(sam_dir, sam_name)
            # cie will always be first paths in the list so if csi does not
            # exist switch will never change from -1 and old element
            # number will always be from the cie
            jobs.append(Search(BTI=BTI, con_file=con,
                               out_file=sam_file, num_old_els=num_old_els))

            self.jobs = jobs  # jobs now stored in the run object

    def run_jobs(self, cur=True):
        '''
        Runs the jobs and stores sam file objects. Probably need a way to
        keep track of which sams are solo and which are intact.
        '''
        sam_dir = self.write_dirs[1]  # stored at 1 index always

        for job in self.jobs():
            if cur:  # search current databases used for remap
                job.search_BTI(self.cur_BDB, self.cur_acc)
            else:  # used for backmap
                job.search_BTI(self.old_BDB, self.old_acc)
        # this gives you a bunch of sam files but now dont know which are solo
        # and which are not done by adding intact or solo to file name above
