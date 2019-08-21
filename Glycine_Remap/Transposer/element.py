
class Element:

    def __init__(self, name, accession, startLocation, endLocation, length, status, seq,):
        self.name = name
        self.accession = accession
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

    def __eq__(self, other):
        if self.name == other.name
        and self.accession == other.accession
        and self.status == other.status:
            return True
        else:
            return False
