# temp file for backmapping methods
# stuff in here is likely to move
import csv
from Transposer.element import Element
from Transposer.search import get_seq_flanks
from Transposer.search import search_BDB


def test_unanchored(s):
    '''
    Test for reading in soy elements in soybase tab delim summary file.
    Checks if element is unanchored
    '''
    s = s[2:]
    try:
        s = int(s)
        return s  # remove GM header portion
    except ValueError:
        return 0
        # was not castable likely unanchored

def translate_chr(chr_dict, chr):
    try:
        return chr_dict[chr]
    except KeyError:
        return None


def make_chr_dict(acc_path):
    acc_dict = {}
    with open(acc_path) as names:
        next(names)
        for line in names:
            l = line.strip().split('\t')
            acc_dict[int(l[0])] = l[1]
    return acc_dict


def make_soy_elements(csv_path, acc_path, delim='\t'):
    '''
    Reads in a tab delim by defualt file that contains element descriptions
    in the soybase format and converts entries to element objects.
    '''
    elements = []
    chr_dict = make_chr_dict(acc_path)
    print(chr_dict[1])
    # TODO: add csv and backmap arguements
    with open(csv_path) as cp:
        reader = csv.reader(cp, delimiter=delim)
        next(reader)
        for row in reader:
            try:

                chr = test_unanchored(row[8])
                acc = translate_chr(chr_dict, chr)
                elements.append(Element(name=row[0], accession=acc,
                                        chr=chr, startLocation=row[10],
                                        endLocation=row[11],
                                        length=int(row[11]) - int(row[10]),
                                        status=row[7], seq=None))
            except IndexError as e:
                continue
    return elements


def add_flanks(els, BDB):
    '''
    Searches and retrieves flanks for elements created from csv file
    '''
    for e in els:
        e.seq, e.left, e.right = get_seq_flanks(e.startLocation,
                                                e.endLocation, e.acc, BDB)


def make_flank_dict(els):
    '''
    Returns a dictionary where keys are concat strings of left and flanking
    sequences and values are any elements that have those exact flanks. This
    is done in the case two elements do have same flanks since this is possible
    but highly unlikely. The flank_dict can then be used to find exact matches
    to flanking sequences from the outdated elements.
    '''
    return {el.left + el.right: el for el in els}


def test_exact(flank_dict, els):
    # m = matches, nmo = non matches old, nmn = non matches new
    m, nm = [], []
    for el in els:
        f = el.left + el.right
        if f in flank_dict:
            m.append((el, flank_dict.pop(flank_dict[f])))
        else:
            nm.append(el)
    return (m, nm)


def chr_indexer(flank_dict):
    '''
    Takes in list of old elements and returns a list of lists. Accessing
    index of list will equal to chromosome number -1. Base 0.
    '''
    max_chr = 0
    for f, el in flank_dict.items():
        if el.chr > max_chr:
            max_chr = el.chr

    index = [[] for i in range(0, max_chr)]
    for f, el in flank_dict.items():
        index[el.chr - 1].append(el)

    return index
    # since exact matches has already been run the index list will only contain
    # old elements that no exact matches were found to


def get_tf(el):
    return el.left + el.right


def rum_ham(index, nm, m, dist=5):
    for el in nm:
        f = get_tf(el)
        for i, oel in enumerate(index[el.chr - 1]):
            if test_ham(f, get_tf(oel)) <= dist:
                m.append((el, index[el.chr - 1].pop(i)))
                break
        # no matches found at this point
        nm.append(el)
    return (m, nm, index)
    # index can be used to get old elements with no matches

# need to also keep non matches old and new seperated
# could just have two different lists


def test_ham(flank_a, flank_b):
    diffs = 0
    for a, b in zip(flank_a, flank_b):
        if a != b:
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


def make_backmap(el_list, flank_dict):
    map = []
    for el in el_list:
        fs = get_flanks(el)
        m = test_exact_flanks(flank_dict, fs)
        if m != False:
            map.append(tuple([el, m]))

    return map


def parse_soy_header(header):
    h = header.split(' ')
    name = h[0][1:0]  # remove > from name
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


def make_soy_element(header, seq, acc_dict=None):

    return Element(name, acc, chr, start, end, length, status, seq)


#print(make_soy_element('>name=RLG_Gmr3_Scaffold328-1 Reference=Du et al. 2010 BMC Genomics 2010, 11:113 Class=I Sub_Class=I Order=LTR Super_Family=Gypsy Family=Gmr3 Description=SOLO', seq='AAA'))


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
