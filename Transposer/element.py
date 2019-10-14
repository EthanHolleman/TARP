
class Element:


    def __init__(self, name, accession, chr, startLocation, endLocation, length, status, seq,):
        self.name = name
        self.accession = accession
        self.chr = int(chr)
        self.startLocation = int(startLocation)
        self.endLocation = int(endLocation)
        self.length = length
        self.status = status
        self.seq = seq

    def toStringFasta(self):
        # sends element to 2 line string to be added to fasta file
        return ">{},{},{},{},{},{}".format(self.name,
                                           self.startLocation,
                                           self.endLocation,
                                           self.length,
                                           self.status,
                                           "\n" + self.seq)
