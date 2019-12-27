import csv
from backmap.feature import feature


def gene_reader(gene_txt, acc_to_chr_dict, delim='\t'):
    genes = []
    with open(gene_txt) as gene:
        reader = csv.reader(gene, delimiter=delim)
        reader.next()
        for row in reader:
            ID, chr, pos, type = row[2], acc_to_chr_dict[row[10]], row[11], 'GENE'
            genes.append(feature(ID, chr, pos, type))
    return genes



def vcf_reader(vcf_file, header='##'):
    '''
    Reads a vcf formated file and yeilds feature objects for each SNP. Each
    feature contains the ID, chromosome, position and a type which identifies
    the feaure as a SNP.
    '''
    with open(vcf_file) as snps:
        features = []
        reader = csv.reader(snps, delimiter='\t')
        head = True
        for i, row in enumerate(reader):
            if row[0][0:2] != header:
                if head:  # avoid including the header row
                    head = False
                else:
                    ID, chr, pos, type = row[2], row[0], row[1], 'SNP'
                    features.append(feature(ID, chr, pos, type))
        return features
