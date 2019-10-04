# class for blast data bases
import subprocess

def make_acc_dict(acc):
    try:
        chrs = {}
        with open(acc) as names:
            for line in names:
                line = line.strip()
                line = line.split("\t")
                number, CM = line
                chrs[CM] = number
        return chrs

    except FileNotFoundError as e:
        return e


class Blast_DB():
    def __init__(self, path, accession_path):
        self.path = path
        self.acc = accession_path
        self.acc2chr = make_acc_dict(self.acc)


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
