import os

from Clustering.big_hit import run_cd_hit
from fasta_tools import check_formating


def get_intact_length(con_file_path):
    check_formating(con_file_path)
    length = []
    with open(con_file_path) as con:
        length = con.readlines()

    return len(length[1])

def make_clstr_name(path, clstr_dir):
    base = os.path.basename(path).split('.')[0]
    return os.path.join(clstr_dir, base)

def make_clstr(fasta, clstr_dir):
    clstr_file_name = make_clstr_name(fasta, clstr_dir)
    print('Clustering', fasta)
    run_cd_hit(clstr_file_name, fasta)
    return clstr_file_name + '.clstr'
