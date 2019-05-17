import subprocess


def createConsensusFile(element_list, fam):
    #creates the fasta file of elements to be passed to clustal omega
    status = element_list[0].status #all elements should be of same status
    file_name = '{}_{}_consensus.fasta'.format(fam,status)
    with open(file_name, 'w') as fasta:
        for element in element_list:
            fasta.write(element.toStringFasta())

    return file_name

def getRandomElements():
    #selects 7 random elements of the family to make fasta file of
    #list is passed to createConsensusFile
    pass

def clustalize():
    pass
    #runs the clustal omega command and returns the conensus sequence

def makeFamConsensus():
    pass
    #runs clustalize for both intact and solo elements
    #creates two consensus sequence fasta files
