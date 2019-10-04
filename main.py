# this is file that runs everything from it
# args file will also be on this outer level
# need to refactor glycine remap so it accepts normal paramters instead
# of arg parser arguements as args will be imported in this file and
# then used here
# use relative imports within the packages and absolute in main script

import os
from Clustering.ClstrFile import ClstrFile
from args import set_args
from Clustering.big_hit import get_element_files
from Clustering.big_hit import HIT_EM_WITH_IT
from Clustering.big_hit import make_consensus_dirs
args = set_args()

def main():

    full_out = os.path.join(args.out, args.name)
    if os.path.isdir(full_out) is True:
        print(args.name, 'already exists')
        # get input to overwrite where
    else:
        os.mkdir(full_out)

    element_files = get_element_files(args.old_elements)

    # make element dirs
    ef_paths = []  # storage for all element paths under full out
    for ef in element_files:
        ef_path = os.path.join(full_out, os.path.basename(ef))
        os.mkdir(ef_path)
        ef_paths.append(ef_path)

    for ef_path, ef_og = zip(ef_path, element_files):
        # iterate through all old elements with the location of where
        # their data should be written


    clstrfile_collection = [ClstrFile(path) for path in clstr_paths]

    for cf in clstrfile_collection:
        cd.write_cluster_fastas

    # need to write everything to correct dirs





    # add arg run name
    # create run name folder in args.output
    # for each collection of old elements create a folder based on the name
    # of those element files (should be given as a dir)
    # for each of those files run clstring create dir within the collection dir
    # contains clstr files.
    # make cf_files in another dir at the same level
    # make clstr_consensus seqs in another dir at same level
    # bowtie in par those clstr consensus
    # write output to another dir at same level
    # end structure should be
    # USER GIVEN DIR
    #|---Run dir (named after args.run_name)
    #   |---Collection A
    #   | -- Collection B  old elements
    #       | -- clstr results
    #       | -- clstr_fastas
    #       | -- clstr cons
    #       | -- bowtie results
