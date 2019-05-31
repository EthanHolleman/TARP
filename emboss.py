import os


def embosser(clustalized_file):
    # takes the file produced by clustalize to create conensus file
    file_parts = clustalized_file.split('_')
    consensus_name = file_parts[0] + file_parts[1] + '_consensus.fna'
    command = 'em_cons -sequence {} -outseq {}'.format(clustalized_file, consensus_name)
    os.system(command)


def formatConsensus(consensus_file):
    # reformats consensus from embosser to fasta file format
    lines = []
    with open(consensus_file, 'r') as con:
        lines = con.readlines()
    with open(consensus_file, 'w') as con:
        for line in lines:
            if line[0] == '>':
                con.write('> ' + consensus_file + '\n')
            else:
                con.write(line.strip())
