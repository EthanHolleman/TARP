from Clustering.element import ClstrElement
from fasta_tools import read_as_tuples


class Cluster():

    def get_average_similarity(self):
        sum_sim, l = 0, 0
        if self.num_elements != 0:
            for el in self.elements:
                sum_sim += el.similarity
            return sum_sim / self.num_elements
        else:
            return 0

    def __init__(self, num=None, parent_file=None):
        self.num = num.strip()
        self.elements = set([])
        self.parent_file = parent_file
        self.num_elements = 0
        self.average_similarity = self.get_average_similarity()
        self.fasta = None
        self.consensus = None
        # num = cluster number, elements = elements in the cluster given by
        # cd-hit file info, parent_file = clstr file cluster came from
        # num elements = number elements in the cluster,
        # fasta is the fasta file containing all elements in the cluster
        # this will only be present if the write_cluster_fastas method is
        # called and cluster is in the cluster set of the clstr_file object

    def add_element(self, el):
        if isinstance(el, ClstrElement):
            self.elements.add(el)
            self.num_elements += 1
            return True
        else:
            return False

    def get_cluster_elements(self):
        # need to search for elements in the fasta file
        # going to assume that headers is always made to the first space
        fasta = self.parent_file.split('.')[0]
        # fasta duplicate always made by cd hit which has the same
        # file name but no extension as the clstr file
        try:
            lines = []
            fasta_elements = []
            search_dict = {}

            lines = read_as_tuples(self.parent_file.split(
                '.')[0])  # read as tuples opens the file

            for element_tuple in lines:
                search_dict[element_tuple[0].split(' ')[0]] = element_tuple
            # make the search dictionary

            # sets have lengths greater than one
            cluster_hits = 0
            for clstr_element in self.elements:
                if len(self.elements) > 1:
                    print(len(self.elements))

        except FileNotFoundError as e:
            return e
