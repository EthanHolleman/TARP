import os
import subprocess

import fasta_tools as ft
from Transposer.interp import toSortedBam


def search(con, index, outputFile, verbose=False):
    '''
    uses a consensus sequence and bowtie2 index to preform a global allignment
    to the given index, returns results as a sam file
    '''

    try:
        commandSearch = "bowtie2 -x {} -f {} -a --non-deterministic -S {}".format(index, con, outputFile)
        print(commandSearch)
        samCommand = "Samtools view -S -b {} ".format(outputFile)
        os.system(commandSearch)
        outputFile = toSortedBam(outputFile, verbose)

        return outputFile
        # will run the command and return a samfile
    except FileNotFoundError:
        print("FileNotFoundError at search in search.py")


def score_consensus(con, old_index, old_elements, verbose=False):
    TEMP_FILE = os.path.join(os.getcwd(), 'temp.fasta')
    '''
    This function will attempt to validate that a consensus sequence will return
    good results when looking for elements in the new assembly by searching
    the consensus against the old assembly and comparing the elements returned
    against the old elements. The assumption being that if it can yeild the
    old elements it should be able to find the new elements.
    '''
    bowtie_results = toSortedBam(search(con, old_index, TEMP_FILE), verbose)
    # need to change the toTxt becuase output will not be in a fasta format so
    # it will be better to read number of elements using a subprocess command
    allign_cmd = ['ssamtools', 'view', bowtie_results, '|', 'wc', '-l']
    num_alligns = subprocess.check_output(allign_cmd)
    num_old_elements = len(ft.read_as_tuples(old_elements))

    return len(elements / old_elements)  # acts as score

#def verify_consensus(score, cutoff):
    # could create tiers of warnings where beyond some point it
    # does not actually return anything
    #samtools view sample.sorted.bam | wc -l
    #print subprocess.check_output(["ping", "-c", "1", "8.8.8.8"])

    # could have two settings for this method, high intesity and low intensity.
    # one would check just the number of elements returned vs the number of
    # old elements and give a score just based on that. Have some kind of
    # cut off or warning based on the score. Other method would compare the
    # seq's % identity to which would be much more comp intesive
