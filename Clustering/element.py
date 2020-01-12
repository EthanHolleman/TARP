
class ClstrElement():

    def __init__(self, name, cluster_num, ID, similarity, nt, rep=False):
        self.cluster = cluster_num
        self.ID = ID
        self.similarity = similarity
        self.nt = nt
        self.rep = rep
        self.name = name

    def __repr__(self):
        (f'{self.__class__.__name__}('
           f'{self.ID}, {self.similarity})')
