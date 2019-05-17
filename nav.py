from element import Element

def readElements(big_File):
    superfam_dict = {}
    seq = ''
    header = ''

    with open(big_File) as big:
        for line in big:

            if line[0] == '>':

                if seq != '':
                    name,start,end,length,status,superfam = headerParser(header)

                    if superfam in superfam_dict:
                        superfam_dict[superfam].append(Element(
                        name,'NONE',start,end,length,status,seq))
                    else:
                        superfam_dict[superfam] = [Element(
                        name,'NONE',start,end,length,status,seq)]
                    seq = ''
                header = line

            else:
                seq = seq + line.strip()

    return superfam_dict


def headerParser(header):
    header = header.split(' ')
    locs = header[-1].split(':')[-1].split('..')

    name = header[0].split('_')[-1]
    start = locs[0]
    end = locs[1]
    length = int(end) - int(start)
    status = header[-2].replace('Description=','')
    superfam = header[-4].replace('Super_Family=','')

    return (name,start,end,length,status,superfam)

#clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v
