import argparse
import multiprocessing as mp

def set_args():
    parser = argparse.ArgumentParser(description='Gylcine_Remap Args')

    parser.add_argument('-confiles', type=str, help='Path to the dir containg consensus files')
    parser.add_argument('-old', type=str, help='Path to outdated element files')
    parser.add_argument('-prev_BDB', type=str, help='Path to BLAST DB created from assembly old elements are mapped to')
    parser.add_argument('-cur_BDB', type=str, help='Path to BLAST DB created from most recent assembly')
    parser.add_argument('-out', type=str, help='Path where results will be written')
    parser.add_argument('-cores', type=int, defualt=mp.cpu_count(), help='Number of cores to run processes on, defualt is all')
    parser.add_argument('-BT_index', type=str, help='Bowtie index created from most recent assembly')
    parser.add_argument('-Old_acc', type=str, help='Path to accession file for outdated assembly')
    parser.add_argument('-Cur_acc', type=str, help='Path to accession file for current assembly')

    return parser
