
class Element:

    def __init__(self, name, accession, chr, startLocation, endLocation, length,
                 status, seq, left=None, right=None):
        self.name = name
        self.accession = accession
        self.chr = int(chr)
        self.startLocation = int(startLocation)
        self.endLocation = int(endLocation)
        self.length = length
        self.status = status
        self.seq = seq
        self.left = left
        self.right = right

    def __repr__(self):
        return (f'{self.__class__.__name__}('
                      f'{self.name}, {self.startLocation})')

    def get_header(self):
        '''
        Returns formated header for writing to fasta files.
        '''
        return '>{} {}: {}, {}'.format(self.name, self.status,
                                       self.startLocation, self.length)

    def get_row(self):
        '''
        Returns row information as list for csv writer obejcts.
        '''
        return [self.name, self.accession, self.chr, self.startLocation,
                self.length, self.status, self.seq, self.left, self.right]
