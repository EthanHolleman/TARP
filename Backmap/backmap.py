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
    # TODO: add csv and backmap arguements
    with open(csv_path) as cp:
        reader = csv.reader(cp, delimiter=delim)
        # next(reader)
        for i, row in enumerate(reader):
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
    for i, e in enumerate(els):
        e.seq, e.left, e.right = get_seq_flanks(e.startLocation,
                                                e.endLocation, e.accession, BDB, 20)

    return els


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
    return (m, nm, flank_dict)


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


def rum_ham(index, nm, m, dist=10):
    # old elements are stored in the index if they are matched they are poped
    # out
    nnm = []
    for el in nm:
        f = get_tf(el)  # both flanks as one string
        try:
            old_elements = index[el.chr - 1]
            if old_elements:
                for i, oel in enumerate(old_elements):
                    if test_ham(f, get_tf(oel)) <= dist:
                        m.append((el, old_elements.pop(i)))
            else:
                nnm.append(el)
        except IndexError:
            nnm.append(el)  # chromosoe of new element not in old el index
            continue

    return (m, nnm, index)
    # index can be used to get old elements with no matches

# need to also keep non matches old and new seperated


def test_ham(flank_a, flank_b):
    diffs = 0
    for a, b in zip(flank_a, flank_b):
        if a != b:
            diffs += 1
    return diffs


def el_chr_dict(elements):
    '''
    Given a list of elements returns a dictionary; key is chromosome and
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


def make_element(parser, header, seq, BDB):
    n, chr, acc, s, e, l, s = parser(header)
    seq, left, right = get_seq_flanks(s, e, acc, BDB=BDB)
    return Element(name=n, accession=acc, chr=chr, startLocation=s,
                   endLocation=e, length=l, status=s, seq=seq,
                   left=left, right=right)


def make_general_element(delim, header, seq):
    h = header.split(delim)
    return Element(h[0], h[1], h[2], h[3], h[4], h[5], h[6], seq)
