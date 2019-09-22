
import subprocess
from multiprocessing import Pool
import os

TEMPLATE_CMD = 'python3 Transposer/remap.py -i {} -s {} -l {} -k {} -a {} -o {} -name {} -c {} -p {} -m {} -e {}'
# # TODO: Transposer will not be called with subprocess and instead will
        # be a submodule so will only be a funciton call

def get_path(dir):
    '''
    Returns the path to a file / dir without the basename.
    Intended to be used for putting files in the directory as
    the file passed in using os.path.join.
    '''
    return '/'.join(dir.split('/')[1:-1])


def linker(old_element_file):
    consensus = 'consensus'
    extension = '.fna'
    o_e = old_element_file.split('.')[0]
    solo = '{}_{}SOLO{}'.format(consensus, o_e, extension)
    intact = '{}_{}INTACT{}'.format(consensus, o_e, extension)

    return tuple([solo, intact])


def make_match_dict(consensus_dir, old_elements_dir):
    '''
    Takes list of all files in the outdated elements dictionary and returns a
    dictionary where outdated element name is key and value is list of what
    the file names for the solo and intact consensus files for that element
    family should be assuming they were created using fasta_tools.
    Want to make the paths complete.
    '''
    match_dict = {}

    con_path = get_path(consensus_dir)  # get paths for both dirs
    old_el_path = get_path(old_elements_dir)

    con_files_set = set(os.listdir(consensus_dir))  # list files in dirs
    old_el_files = os.listdir(old_elements_dir)

    for old_el in old_el_files:
        # all items in dictionary are joined to any path leading up to the
        # actual consensus_dir and old_elements_dir directories. This way
        # the dictionary can be more easily be used for creating the list
        # of process commands to run
        key = os.path.join(old_el_path, old_el)
        match_dict[key] = [None, None]
        solo, intact = linker(old_el)

        if solo in con_files_set:
            match_dict[key][0] = (os.path.join(con_path, solo))
        if intact in con_files_set:
            match_dict[key][1] = (os.path.join(con_path, intact))
        # Solo element will always be at i 0 and intact will be at i 1
        # no matter if one or both are present for that family

    return match_dict


def get_processes(match_dict, old_acc, cur_acc, bowtie_index, output, cur_BDB, old_BDB):
    commands = []
    for file in match_dict:  # file is old elements paths
        # when passing the con arguments if only one is in the list need to
        # determine which it is and pass it to the correct arguement
        solo_con, intact_con = tuple(match_dict[file])
        commands.append(package_args(bowtie_index, cur_BDB, old_BDB, cur_acc,
                        old_acc, output, file, intact_con=intact_con,
                        solo_con=solo_con))
    return commands


def package_args(bowtie_index, cur_BDB, old_BDB, cur_acc, old_acc, output,
               old_elements, allowance=150, intact_con=None, solo_con=None,
               verbose=False):
               # defualt args are given here so when packaged all args needed
               # will have a value when passed into the remap functionality
               # intended to reduce issues with using starmap
    '''
    Returns a tuple of all process to be run by formating the TEMPLATE_CMD
    constant. Will need to make sure it can run with only one consensus.
    '''
    name = old_elements.split('.')[0]
    return tuple[bowtie_index, solo_con, intact_con, cur_acc, output, name,
                 cur_BDB, old_BDB, old_acc, old_elements, allowance, verbose]
