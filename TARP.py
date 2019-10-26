#!/usr/bin/env python3
import os
import time
from args import set_args
from args import check_depends
from run import Run
from Transposer.process_sams import write_fasta
from Transposer.process_sams import sort_elements
from Transposer.process_sams import sort_sams
from Transposer.process_sams import rename_elements
from Transposer.process_sams import prune
from Transposer.process_sams import write_csv
# coult potentially include in the Run objects
# # TODO: Need a more comprehensive duplicate removal system
# was getting hits that are less than 10 bases away from each other but
# are still being considered unique because it is not the exact positon
# quick and dirty way to do this would to be for each one added add
# keys that extend out x number of base pairs away ideally then keep the one
# that has fewest ns but basically going to be the same
# also need to experiment with using more sensitivity to remap stuff


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
    min = args.E
    # assign passed arguements not really needed but done
    # to add some clarity later on
    # create the run object
    new_run = Run(cur_BDB=BDB_cur, old_BDB=BDB_old, BTI=BTI,
                  cur_acc=acc_cur, old_acc=acc_old, run_name=run_name,
                  output=output, cie=intact_old, csi=solo_old, min_els=min)

    # do everything up to including bowtieing all the consensus sequences
    print('Made Run')
    new_run.select_clusters()
    print('Selected Clusters')
    new_run.make_clstr_fastas()  # convert clstr files to full fasta files
    print('Made Fastas')
    new_run.make_consensensi()  # make cluster consensus sequences
    print('\nMade Consensus Seqs')
    new_run.make_jobs()  # make bowtie jobs
    print('Made jobs')
    new_run.run_jobs()  # run bowtie jobs add threads arguement
    print('Ran jobs')

    print('Sorting Elements')
    sorted_elements = sort_sams(new_run.jobs)
    print('Pruning Remaining Duplicates')
    sorted_elements = prune(sorted_elements, n=100)
    print('Naming Elements')
    sorted_elements = rename_elements(sorted_elements)
    print('Writing to', new_run.write_dirs[2])

    out_base = os.path.join(new_run.write_dirs[2], run_name)

    write_fasta(sorted_elements, out_base + '.fa')
    write_csv(sorted_elements, out_base + '.csv')
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


if __name__ == '__main__':
    main()

    # then use process sam files to write the outpt
