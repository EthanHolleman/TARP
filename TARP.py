#!/usr/bin/env python3

# once done writing this then change to no extension so can be run with
# just the command ./TARP parameters
import time

from args import set_args
from run import Run
from Transposer.process_sams import write_fasta
from Transposer.process_sams import sort_elements
# coult potentially include in the Run objects

LOGO = './.logo.txt'

def print_logo():
    with open(LOGO) as logo:
        for l in logo.readlines():
            print(str(l.strip()))
    time.sleep(2)


def main():
    args = set_args()  # read args and assign to args
    intact_old = args.I
    solo_old = args.S
    BDB_old = args.P
    BDB_cur = args.C
    output = args.O
    BTI = args.B
    acc_old = args.acc_o
    acc_cur = args.acc_c
    run_name = args.name
    backmap = args.M

    print_logo()

    # assign passed arguements not really needed but done
    # to add some clarity later on

    # create the run object
    new_run = Run(cur_BDB=BDB_cur, old_BDB=BDB_old, BTI=BTI,
                  cur_acc=acc_cur, old_acc=acc_old, run_name=run_name,
                  output=output, cie=intact_old, csi=solo_old)

    # do everything up to including bowtieing all the consensus sequences
    new_run.select_clusters()
    new_run.make_clstr_fastas()  # convert clstr files to full fasta files
    new_run.make_consensensi()  # make cluster consensus sequences
    new_run.make_jobs()  # make bowtie jobs
    new_run.run_jobs()  # run bowtie jobs

    sorted_elements = sort_sams(run.jobs)
    write_fasta(sorted_elements, run.write_dirs[2])
    if backmap == True:
        pass
        # run all the backmapping stuff here
        # this is the flanking sequence stuff that needs to happen
        # create a search object form the old intact and solo elements
        # look at the flanker file

    # ========================
    # Start backmaping process here
    # should add a parameter that tells program to run the backmap portion
    # ========================

    # run.write_dirs[2] is the results file

    # at this point need to get all the output and
    # process those with the process_sams method
    # need to get list of sam files, jobs is a list of search objects
    # need to extract the sam files

    # need to see if process sams methods should be moved to search methods
    # and if they need to be adjusted
    # search objects now contain the element_set that is prunned and typed
    # after it is created
    # need easy way of getting the length of each consensus that was searched

    # needs to work with the different lengths of the consensus sequences
    # need to adjust sort_sams to do this


if __name__ == '__main__':
    main()

    # then use process sam files to write the outpt
