import argparse
import sys

def set_args():
    parser = argparse.ArgumentParser(description='Gylcine_Remap Args')

    parser.add_argument('-I', type=str, help='Path to outdated intact type elements')
    parser.add_argument('-S', type=str, help='Path to outdated solo type elements')
    parser.add_argument('-P', type=str, help='Path to BLAST DB created from assembly old elements are mapped to')
    parser.add_argument('-C', type=str, help='Path to BLAST DB created from most recent assembly')
    parser.add_argument('-O', type=str, help='Path where results will be written')
    parser.add_argument('--T', type=int, defualt=mp.cpu_count(), help='Number of cores to run processes on, defualt is all')
    parser.add_argument('-B', type=str, help='Bowtie index created from most recent assembly')
    parser.add_argument('-acc_o', type=str, help='Path to accession file for outdated assembly')
    parser.add_argument('-acc_c', type=str, help='Path to accession file for current assembly')
    parser.add_argument('-name', type=str, help='Name of the run, all data will be written in dir under this name')

    args = parser.parse_args()
    exit = False
    if not args.I:
        print('Path to intact elements fasta file required (-I)')
        exit = True
    if not args.P or not args.C:
        print('Please supply path to BLAST database for both current \
        and outdated assemblies')
        exit = True
    if not args.O:
        print('Please supply output directory')
        exit = True
    if not args.B:
        print('Please supply bowtie2 index of current assembly')
        exit = True
    if not args.acc_c or not args.acc_c:
        print('Please provide paths to the acc2chr files for the current and \
        outdated assemblies. These can be found in the Genbank FTP folders for \
        a given assembly.')
    if exit:
        sys.exit()


# confiles can be made with one consensus to rule them all from fasta tools
# args that need to add
# run name, potentially need old bowtie index to do the cross validation
# of consensus sequences
# need to fix some of this crazy output from fasta tools
