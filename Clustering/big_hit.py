import os
import subprocess
from fasta_tools import make_consensus
from fasta_tools import check_formating


def get_element_files():
    element_list = []
    els = os.listdir(ELEMENTS)
    for el in els:
        ext = str(el.split('.')[-1])
        if FILE_EXT == ext:
            element_list.append(os.path.join(ELEMENTS, el))
    print('Found', str(len(element_list)), 'files')
    return element_list


def run_cd_hit(output, input_file, identity=0.85):
    cmd = ['cd-hit-est', '-i', input_file, '-o', output,
           '-c', str(identity),  '-T', '6', '-d', '0']

    FNULL = open(os.devnull, 'w')
    try:
        hit = subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        # os.system(cmd)
    except subprocess.CalledProcessError as e:
        return e
    # cd hit for a single file


def make_consensus_clusters(clstr_fastas, output_dir, min_elements=4):
    # clstr_fastas is a dir containing the cluster fasta files
    min_elements *= 2  # number lines / 2 = number elements
    cf_paths = [os.path.join(clstr_fastas, file)
                for file in os.listdir(clstr_fastas)]
    # location of all target files
    con_canidates = []
    for cf_file in cf_paths:
        check_formating(cf_file)
        with open(cf_file) as cf:
            if len(cf.readlines()) >= min_elements:
                con_canidates.append(cf_file)
    for con in con_canidates:
        name = '{}_consensus.fa'.format(os.path.basename(con))
        con_path = os.path.join(output_dir, name)
        make_consensus(con, output_path=con_path,
                       min_elements=min_elements + 1)


def make_consensus_dirs(clstr_dir, output_dir, min_elements=2):
    #ALLOWED_FASTAS = set(['.fna', '.fasta', '.fa', '.fsa'])
    '''
    Extension on the make_consensus_clusters function. Instead of
    making consensus sequences for fasta_clstrs in a single dir, clstr_dir
    should be a directory containing multible sub dirs each with their own
    collection of fasta_clstr files. Output will mirror the file structure of
    the input with one top level dir containing many sub dirs each with their
    respective consensus sequences.
    '''
    if os.path.exists(output_dir) is False:
        os.mkdir(output_dir)

    for dir in os.listdir(clstr_dir):
        abs_dir = os.path.join(clstr_dir, dir)
        if os.path.isdir(abs_dir):
            sub_dir = os.path.join(output_dir, dir)
            if os.path.exists(sub_dir) is False:
                os.mkdir(sub_dir)
            make_consensus_clusters(
                abs_dir, sub_dir, min_elements=min_elements)
