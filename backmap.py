# temp file for backmapping methods
# stuff in here is likely to move
import csv
from Transposer.element import Element
from Transposer.search import get_seq_flanks

def search_BDB(start, end, entry, r_seq=True):
    '''
    search the blast db for a sequence in an entry
    '''

    seq_cmd = 'blastdbcmd -db {} -dbtype nucl -range {}-{} -entry {}'.format(
        self.BDB, start, end, entry)
    try:
        output = subprocess.check_output(seq_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        return ''
    if r_seq:
        return ''.join(str(output).split('\\n')[1:])
        # returns string of just the sequence
    else:
        return output  # else return all output

def get_flanks(element, n=20):
    '''
    Takes in a list of elements and a blast
    DB of the assembly and returns their flanking
    sequences out to n. Flanks returned as a tuple first index is left flank
    second is the right flank.
    '''
    lf = search_BDB(element.startLocation-(n+1), element.startLocation-1, element.accession)
    rf = search_BDB(element.endLocation + 1, element.endLocation + (n+1), element.accession)

    return tuple([lf, rf])

def flank_elements(element_list, n=20):
    '''
    Takes in a list of elements and calls get_flanks on each to find
    flanking sequences for all elements in the list. Returns a zipped
    list of elements and their flanks as a tuple.
    '''
    flanks = []
    for el in element_list:
        flanks.append(get_flanks(el, n=n))

    return zip(flanks, element_list)

def make_flank_dict(zipped_els):
    '''
    Returns a dictionary where keys are concat strings of left and flanking
    sequences and values are any elements that have those exact flanks. This
    is done in the case two elements do have same flanks since this is possible
    but highly unlikely. The flank_dict can then be used to find exact matches
    to flanking sequences from the outdated elements.
    '''
    flank_dict = {}
    for f, el in zipped_els:
        l, r = f  # unpack flank tuple
        f = l + r  # make string of both flanks to use as key for exact matches
        if f in flank_dict:
            flank_dict[f].append(el)
        else:
            flank_dict[f] = [el]

    return flank_dict


def test_exact_flanks(flank_dict, flank_tuple):
    # what to be able to create a list of only the elements that did not
    # get exact matches but keep them in order so can run the next backmap
    # with hamming distance on those
    l, r = flank_tuple
    if l+r in flank_dict:
        return flank_dict[l+r]
    else:
        return False


def test_ham(flank_a, flank_b):
    diffs = 0
    for a, b in zip(flank_a, flank_b):
        if a!= b:
            diffs += 1
    return diffs


def el_chr_dict(elements):
    '''
    Given a list of elements returns a dictionary key is chromosome and
    values are elements in that chromosome. For backmapping is used for organizing
    the old elements. NEEDS TO BE CHANGED VALUE SHOULD BE OLD ELEMENT FLANKS.
    '''
    chr = {}
    for el in elements:
        if el.chr not in chr:
            chr[el.chr] = [el]
        else:
            chr[el.chr].append(el)

    return chr

def rum_ham(elements, old_el_dict, dist=5):
    '''
    Need to compare the flanks of elements here maybe not worth doing the exact
    method and just do the hamming distance instead. Dictionary of old elements
    needs to be for the flanks but could calculate the flanks for new elements
    right away or include the flanks in the element type so dont have to so another
    search and just get the flanks right away.
    '''
    cur_chr = 0
    for el in elements:
        for old_el in old_el_dict[el.chr]:
            tf = el.left + el.right
            otf = old_el.left + old_el.right
            if test_ham(tf, otf) <= dist:
                yield tuple([el, old_el])
                break


def make_backmap(el_list, flank_dict):
    map = []
    for el in el_list:
        fs = get_flanks(el)
        m = test_exact_flanks(flank_dict, fs)
        if m != False:
            map.append(tuple([el, m]))

    return map


def make_chr_dict(acc_path):
    acc_dict = {}
    with open(acc_path) as names:
        for line in names:
            l = line.strip().split('\t')
            acc_dict[l[0]] = l[1]
    return acc_dict

def parse_soy_header(header):
    h = header.split(' ')
    name = h[0][1:0] # remove > from name
    chr = name.split('_')[2][2:].split('-')[0]
    acc = acc_dict(chr)
    start, end = tuple(h[-1].split(':')[-1].split('..'))
    length = end - start
    status = h[-2].split('=')[-1]

    return tuple([name, acc, chr, start, end, length, status])


def make_element(parser, header, seq, BDB):
    n, chr, acc, s, e, l, s = parser(header)
    seq, left, right = get_seq_flanks(s, e, acc, BDB=BDB)
    return Element(name=n, accession=acc, chr=chr, startLocation=s,
                   endLocation=e, length=l, status=s, seq=seq,
                   left=left, right=right)



def make_soy_element(header, seq, acc_dict):
    '''
    Parses the header and seq of a soybase element and returns and eleemnt
    object that is ready for backmap processing.
    '''
    h = header.split(' ')
    name = h[0][1:0] # remove > from name
    chr = name.split('_')[2][2:].split('-')[0]
    acc = acc_dict(chr)
    start, end = tuple(h[-1].split(':')[-1].split('..'))
    length = end - start
    status = h[-2].split('=')[-1]

    return Element(name, acc, chr, start, end, length, status, seq)

def make_general_element(delim, header, seq):
    h = header.split(delim)
    return Element(h[0], h[1], h[2], h[3], h[4], h[5], h[6], seq)
# probably should be working with the sorted elements
# that way they could be renamed with the info from fasta writer tool
# want to make that a seperate function and then fasta
# writer just runs through all that
def write_map(map, output):
    '''
    Writes a backmap file which is csv file where each row contains
    old and new elements that have been determined to be equivalent based
    on backmapping.
    '''
    with open(output, 'w') as out:
        writer = csv.writer(out)
        for a, b in map:
            row = [a.name, a.chr, a.startLocation, a.length, a.seq,
                     b.name, b.chr, b.startLocation, b.length, b.seq]
            writer.writerow(row)

    # return the element with the same flank
    # idea being new elements could be put into the flank
    # dict and then loop up old elements or the other way around
    # would work as well
