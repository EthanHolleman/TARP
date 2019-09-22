from multi_process import run_pool
from pairing import make_match_dict
from pairing import get_processes
from args import set_args

from Clustering.big_hit import HIT_EM_WITH_IT
from Clustering.big_hit import get_element_files
from Clustering.ClustrFile import ClstrFile
def main():
    args = set_args()

    # take path containing fastas of old elements
    # need arguement for where to store the clustered files
    # make clustered files and create consensus sequences
    clstr_paths = HIT_EM_WITH_IT(get_element_files(args.old_elements))  # old elements should be a directory
    # need to have a way of getting location to cd-hit program
    # probably do this in a config type file or run a setup script
    clstr_collection = [ClstrFile(path) for path in clstr_paths]
    # make collection of ClstrFile objects in a list
    for clstr_file in clstr_collection:
        clstr_file.write_cluster_fastas()











    match_dict = make_match_dict(args.confiles, args.old)
    args_iter = get_processes(match_dict, args.old_acc, args.cur_acc,
                              args.BT_index, args.out, args.cur_BDB,
                              args.old_BDB)
    run_pool(args_iter, cpu_count=args.cores)
