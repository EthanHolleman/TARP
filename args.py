import argparse
import shutil
import sys
import multiprocessing as mp

DEPENDS = './depends.txt'
LOGO = './.logo.txt'


def print_logo():
    '''
    Prints out the logo from LOGO path and checks for
    depencendies.
    '''
    with open(LOGO) as logo:
        for l in logo.readlines():
            print(str(l.strip()))
    check_depends()


def check_depends():
    dps = []
    with open(DEPENDS) as dp:
        dps = dp.readlines()
    for dp in dps:
        if shutil.which(dp) is not None:
            print('Please install {} and rerun TARP'.format(dp))
            sys.exit

    print('All dependency checks passed\n')


def set_args():
    parser = argparse.ArgumentParser(description='Gylcine_Remap Args')

    parser.add_argument(
        '-I', type=str, help='Path to outdated intact type elements')
    parser.add_argument(
        '-S', type=str, help='Path to outdated solo type elements')
    parser.add_argument(
        '-P', type=str, help='Path to BLAST DB created from assembly old elements are mapped to')
    parser.add_argument(
        '-C', type=str, help='Path to BLAST DB created from most recent assembly')
    parser.add_argument(
        '-O', type=str, help='Path where results will be written')
    parser.add_argument('-T', type=int, default=mp.cpu_count(),
                        help='Number of cores to run processes on, default is all')
    parser.add_argument(
        '-B', type=str, help='Bowtie index created from most recent assembly')
    parser.add_argument('-acc_o', type=str,
                        help='Path to accession file for outdated assembly')
    parser.add_argument('-acc_c', type=str,
                        help='Path to accession file for current assembly')
    parser.add_argument(
        '-name', type=str, help='Name of the run, all data will be written in dir under this name')
    parser.add_argument(
        '-M', default=False, help='Tells TARP what kind of backmap to run if any. 1 = distance minimization, 2 = flanking, 3 = feature. 3 requires the arguement -F to have a valid path.')
    parser.add_argument(
        '-sum', default=False, help='Path to element summary file, currently required to be in soybase format. Must be given to run any backmap operation.')
    parser.add_argument('-E', default=2, help='Min number of elements required in a cluster for a consensus to be made and searched. Increasing this value will decrease runtime but also decrease the number of elements found.')
    parser.add_argument(
        '-F', default=False, help='Path to file containing genomic features backmaping. Required only for feature based backmap')
    args = parser.parse_args()
    print_logo()
    exit = False
    if not args.I:
        print('Path to intact elements fasta file required (-I)')
        exit = True
    if not args.P or not args.C:
        print(
            'Please supply path to BLAST database for both current and outdated assemblies')
        exit = True
    if not args.O:
        print('Please supply output directory')
        exit = True
    if not args.B:
        print('Please supply bowtie2 index of current assembly')
        exit = True
    if not args.acc_c or not args.acc_c:
        print('''Please provide paths to the acc2chr files for the current and
        outdated assemblies. These can be found in the Genbank FTP folders for
        a given assembly.''')
    if args.M:
        if args.sum == False:
            print('Please provide a soybase formated element summary file to run backmapping (-sum)')
            exit = True
        elif args.M == '3' and args.F == False:
            print('Please provide a genomic feature file (Genes or SNPs) to -F to run feature backmapping')
            exit = True
    if exit:
        print('Goodbye')
        sys.exit()
    else:
        return args
