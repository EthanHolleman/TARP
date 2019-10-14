#!/usr/bin/env python3

# once done writing this then change to no extension so can be run with
# just the command ./TARP parameters

from args import set_args
from run import Run
from Transposer.process_sams import write_fasta
# coult potentially include in the Run objects

def main():
    args = set_args() # read args and assign to args
    intact_old = args['I']
    solo_old = args['S']
    BDB_old = args['P']
    BDB_cur = args['C']
    output = args['O']
    BTI = args['B']
    acc_old = args['acc_o']
    acc_cur = args['acc_c']
    run_name = args['name']

    # assign passed arguements not really needed but done
    # to add some clarity later on

    # create the run object
    new_run = Run()

    # run operations refer to tests done here

    # then use process sam files to write the output
