import os
from interp import toSortedBam


def search(con, index, outputFile, verbose):
    '''
    uses a consensus sequence and bowtie2 index to preform a global allignment
    to the given index, returns results as a sam file
    '''

    try:
        commandSearch = "bowtie2 -x {} -r {} -a --non-deterministic -S {}".format(index, con, outputFile)
        samCommand = "Samtools view -S -b {} ".format(outputFile)
        os.system(commandSearch)
        outputFile = toSortedBam(outputFile, verbose)

        return outputFile
        # will run the command and return a samfile
    except FileNotFoundError:
        print("FileNotFoundError at search in search.py")


def makeIndex(name):
    pass
