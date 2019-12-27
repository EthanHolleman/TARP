import math

def feature_sort(feature_iter, key=None, num_chr=20):
    '''
    Sorts features first into bins based on chromosome number and then by
    position within that chromosome. Returns a list of lists increasing indexes
    correspond to increasing chromosome number.
    '''
    chr_bins, cur_chr = [[] for i in range(0, num_chr+1)], None
    print(chr_bins)
    print(len(feature_iter), 'number of features')
    for i, chr in enumerate(chr_bins):
        print(len(chr), 'number of features in chr before add {}'.format(i))
        s = 0
    for f in feature_iter:
        if f.chr:
            if f.chr >= 0 and f.chr <= num_chr:
                chr_bins[f.chr-1].append(f)
                if f.chr == 1: s += 1
            else:
                print(type(f.chr))
    print(s, 'NUMBER CHR 1')
    for i, chr in enumerate(chr_bins):
        print(len(chr), 'number of features in chr before sort {}'.format(i))

    for bin, _ in enumerate(chr_bins):
        if key == None:
            chr_bins[bin] = sorted(chr_bins[bin], key=lambda x: x.position)
        else:
            chr_bins[bin] = sorted(chr_bins[bin], key=key)
    length_bins = sum([len(b) for b in chr_bins])
    for i, chr in enumerate(chr_bins):
        print(len(chr), 'number of features in chr {}'.format(i))
    return chr_bins


def feature_flank_search(l, start, end, chr):
    '''
    Search through an array of features using the start and end positions of
    an element. Returns the least distant flanking features to the element.
    '''
    chr -= 1  # convert to index
    low, high, midpoint = 0, len(l[chr])-1, None
    print(high, 'number elements in there')
    while low <= high:
            midpoint = int(math.floor((low + high) / 2))
            midpoint_value = l[chr][midpoint].position
            if start < midpoint_value:
                high = midpoint -1
            elif start > midpoint_value:
                low = midpoint+1
    # low and high are flipped at this point low in front of the target high
    # behind the target
    if low == len(l[chr]):
        low -= 1
        high -= 1
        print(low,high,'end of chromosome')
    else:
        print(l[chr][low].position, l[chr][high].position, 'positions in search')
        while l[chr][high].position < end:
            low+= 1
    return l[chr][high], l[chr][low]


def is_ambigous_feature(chr, left_feat_ind, right_feat_ind, index):
    '''
    If there is more than one element between the two features it is
    ambigous and further matching is required to make a definitive match
    between elements.
    '''
    c = 0
    print(left_feat_ind, right_feat_ind, 'feature locations')
    for el in index[chr -1]:
        if el.startLocation > left_feat_ind and el.startLocation < right_feat_ind:
            c+=1
            print(c, el.name, el.startLocation)
            if c > 1:
                return False
    return True




def compare_flanking_features(el_a, el_b, features, index):
    '''
    Searches a list of features for the flanking features of elements a and b.
    Then compares the flanking features for each element. If they are the same
    will return true otherwise returns false.
    '''
    print(el_a.startLocation, 'element a start location')
    left_a, right_a = feature_flank_search(features, el_a.startLocation,
                                           el_a.endLocation, el_a.chr)
    left_b, right_b = feature_flank_search(features, el_b.startLocation,
                                           el_b.endLocation, el_b.chr)
    if left_a == left_b and right_a == right_b:
        amb = is_ambigous_feature(el_a.chr, left_a.position, right_a.position, index)
        print(amb, 'AMB\n\n\n\n')
        return True
    else:
        return False

def anchored_backmap(index, nm, m, features):
    '''
    Index = chr index from backmap.py, nm= non matching elements, m = list of
    matching elements.
    '''
    nnm = []
    for el in nm:  # iterate through non-matched elements
        for i, oel in enumerate(index[el.chr-1]):
            if compare_flanking_features(el, oel, features, index):
                m.append((el, index[el.chr-1].pop(i)))
        nnm.append(el)
    return (m, nnm, index)
