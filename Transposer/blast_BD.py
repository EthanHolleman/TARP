# class for blast data bases
import subprocess


class Blast_DB():
    def __init__(self, path, accession_path):
        self.path = path
        self.acc = accession_path


    def search(self, start, end, entry, r_seq=True):
        '''
        search the blast db for a sequence in an entry
        '''

        seq_cmd = 'blastdbcmd -db {} -dbtype nucl -range {}-{} -entry {}'.format(
                   self.path, start, end, entry)
        try:
            output = subprocess.check_output(seq_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            return ''
        if r_seq:
            return ''.join(str(output).split('\\n')[1:])
            # returns string of just the sequence
        else:
            return output  # else return all output
