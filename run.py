import os

from Transposer.search import Search
from par.multi_process import run_pool
from Clustering.ClstrFile import ClstrFile
from Clustering.big_hit import run_cd_hit


from fasta_tools import make_consensus
    # return the length of the sequence which will
    # be in the second line


def make_clstr(fasta, clstr_dir):
    clstr_file_name = make_clstr_name(fasta, clstr_dir)
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

    def __init__(self, cur_BDB, cur_acc, old_BDB, old_acc, BTI, cie, run_name, output, csi=None):
        self.write_dirs = make_dirs(output, run_name)
        self.old_BDB = old_BDB, old_acc
        self.cur_BDB = cur_BDB, cur_acc
        self.BTI = BTI  # path to bowtie index
        self.cie = cie  # path to file
        self.csi = csi  # path to file
        self.old_acc = old_acc  # path to acc2chr file for old assembly
        self.cur_acc = cur_acc  # path to acc2chr file for new assembly
        self.cie_clstrs = ClstrFile(
            path=make_clstr(self.cie, self.write_dirs[0]))
        self.num_cie = num_cie = sum(1 for line in open(self.cie))/2
        if csi is not None:
            self.csi_clstrs = ClstrFile(
                path=make_clstr(self.csi, self.write_dirs[0]))
            self.num_csi = num_cie = sum(1 for line in open(self.csi))/2
        else:
            self.csi_clstrs = None
            self.num_csi = 0
        self.cie_fastas = None
        self.csi_fastas = None
        self.cie_cons = None
        self.csi_cons = None
        self.jobs = None

    def select_clusters(self, min_elements=10):
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

    def make_consensensi(self, min_elements=5, n=21):
        '''
        Need a way to get the paths to the fasta files and then process
        them with clustal omega and that kind of thing. Need to store the
        paths to consensus sequence in the self.con variables in run.
        '''
        cie_cons, csi_cons = [], []
        for fasta in self.cie_cons:
            con_name = fasta + '_consensus'
            make_consensus(fasta, con_name, min_elements, n=n)
            cie_cons.append(con_name)
        self.cie_cons = cie_cons

        if self.csi_cons is not None:
            for fasta in self.csi_cons:
                con_name = fasta + '_consensus'
                make_consensus(fasta, con_name, min_elements, n=n)
                csi_cons.append(con_name)
            self.csi_cons = csi_cons

    def make_jobs(self):
        '''
        generates a list of bowtie commands to be run based on the number of
        consensus clusters in the run object. Each consensus cluster will be
        bowtied and then needs to tripped of duplicate elements. Solo and intact
        will need to be seperated. Will be a list with first list being intact
        jobs and second list being solo jobs.
        '''
        sam_dir = self.write_dirs[1]
        jobs = []
        for c in self.cie_cons:
            sam_name = os.path.basename(c).split('.')[0] + '_intact.sam'
            sam_file = os.path.join(sam_dir, sam_name)
            jobs.append(Search(BTI=self.BTI, con_file=c,
                               out_file=sam_file, num_old_els=self.num_cie,
                               type='I', acc=self.cur_acc, BDB=self.cur_BDB))
        ## TODO: edit the search objects so taking in all required parameters
        if self.csi_cons is not None:
            for c in self.csi_cons:
                sam_name = os.path.basename(c).split('.')[0] + '_solo.sam'
                sam_file = os.path.join(sam_dir, sam_name)
                jobs.append(Search(BTI=self.BTI, con_file=c,
                                   out_file=sam_file, num_old_els=self.num_csi,
                                   type='S', acc=self.cur_acc, BDB=self.cur_BDB))
            # need numbers of both the intact and solo files
        self.jobs = jobs  # jobs now stored in the run object

    def run_jobs(self, cur=True):
        '''
        Runs the jobs and stores sam file objects. using methods from sam class
        removes duplicates and types the elements for each sam object. After this
        the elements are ready to be placed into a final order and then written to
        a fasta file.
        '''

        sam_dir = self.write_dirs[1]  # stored at 1 index always

        for job in self.jobs:
            if cur:  # search current databases used for remap
                job.search_BTI(self.cur_BDB, self.cur_acc)
                # search object will now hold the element list
            else:  # used for backmap
                job.search_BTI(self.old_BDB, self.old_acc)

            #job.sam.type_elements()  # type elements and remove false solos
            #job.sam.remove_dups()  # remove duplicate elements within the sams

        # this gives you a bunch of sam files but now dont know which are solo
        # and which are not done by adding intact or solo to file name abov

# temp test will be removed and made more unit test like once everything is working
# backmapping elements done by translating the old elements into element objects
# if information is in the the header given then that is one way if it is not
# then can bowtie consensus sequences back to original assembly but can work
# on that later



'''
cur_BDB = '/media/ethan/EH_DATA/GMAX_2.1_BDB_parsed/GM_2.1_BD'
cur_acc = '/media/ethan/EH_DATA/GMax2.1_assembly/chr2acc.txt'
old_acc = '/media/ethan/EH_DATA/GMax1.1_assembly/chr2acc'
old_BDB = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/GMAX_1.0_BDB'
BTI = '/media/ethan/EH_DATA/GMAX_1.1_BTI/GMAX_1.1_BTI'
cie = '/media/ethan/EH_DATA/Gypsy_Seperated/Gmr3INTACT.fna'
csi = '/media/ethan/EH_DATA/Gypsy_Seperated/Gmr3SOLO.fna'


a = Run(cur_BDB=cur_BDB, cur_acc=cur_acc, old_acc=old_acc, BTI=BTI, cie=cie,
        csi=csi, run_name='TEST_OLD_BTI', output='/media/ethan/EH_DATA', old_BDB=old_BDB)
print(a.csi_clstrs)
print(a.cie)

a.select_clusters()
a.make_clstr_fastas()
a.make_consensensi()
a.make_jobs()
print(a.jobs)
a.run_jobs()
'''
