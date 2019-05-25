import subprocess

def embosser(clustalized_file):
    #takes the file produced by clustalize to create conensus file
    consensus_name = clustalized_file.split('_')[1] + '_consensus'
    command = 'em_cons -sequence {} -outseq {}'.format(clustalized_file, consensus_name)
    subprocess.call(command)

def formatConsensus(consensus_file):
    #reformats consensus from embosser to fasta file format
    lines = []
    with open(conensus_file, 'wr') as con:
        lines = con.readlines()

        for line in lines:
            if line[0] == '>': con.write('> ' + consensus_file + '\n')
            else:
                con.write(line.strip())
