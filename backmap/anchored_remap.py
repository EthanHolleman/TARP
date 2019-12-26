import math

def feature_sort(feature_iter, key=None, num_chr=20):
    '''
    Sorts features first into bins based on chromosome number and then by
    position within that chromosome. Returns a list of lists increasing indexes
    correspond to increasing chromosome number.
    '''
    chr_bins, cur_chr = [[]]*num_chr, None
    for f in feature_iter:
        if f.chr != None:
            if cur_chr != f.chr and f.chr <= num_chr and f.chr > 0:
                cur_chr = f.chr
                chr_bins[cur_chr-1].append(f)

    for bin, _ in enumerate(chr_bins):
        if key == None:
            chr_bins[bin] = sorted(chr_bins[bin], key=lambda x: x.position)
        else:
            chr_bins[bin] = sorted(chr_bins[bin], key=key)
    length_bins = sum([len(b) for b in chr_bins])
    print(length_bins, 'Total features after sort')
    return chr_bins


def feature_flank_search(l, start, end, chr):
    '''
    Search through an array of features using the start and end positions of
    an element. Returns the least distant flanking features to the element.
    '''
    chr -= 1  # convert to index
    low, high, midpoint = 0, len(l[chr])-1, None
    print(len(l), 'length of l')
    while low <= high:
            midpoint = int(math.floor((low + high) / 2))
            midpoint_value = l[chr][midpoint].position
            if start < midpoint_value:
                high = midpoint -1
            elif start > midpoint_value:
                low = midpoint+1
    print(l[chr], 'length chromosome')
    print(low, 'low value before while loop')
    if low == len(l[chr]):
        low -= 1
    else:
        while l[chr][low].position < end:
            low+= 1
    print(l[chr][high], start, end, l[chr][low], 'low, start, end, high')
    return l[chr][high], l[chr][low]


def compare_flanking_features(el_a, el_b, features):
    '''
    Searches a list of features for the flanking features of elements a and b.
    Then compares the flanking features for each element. If they are the same
    will return true otherwise returns false.
    '''
    left_a, right_a = feature_flank_search(features, el_a.startLocation,
                                           el_a.endLocation, el_a.chr)
    left_b, right_b = feature_flank_search(features, el_b.startLocation,
                                           el_b.endLocation, el_b.chr)
    if left_a == left_b and right_a == right_b:
        print('MATCH FOUND!')
        return True
    else:
        print('NO MATCH')
        return False

def anchored_backmap(index, nm, m, features):
    '''
    Index = chr index from backmap.py, nm= non matching elements, m = list of
    matching elements.
    '''
    nnm = []
    for el in nm:  # iterate through non-matched elements
        for i, oel in enumerate(index[el.chr-1]):
            if compare_flanking_features(el, oel, features):
                m.append((el, index[el.chr-1].pop(i)))
        nnm.append(el)
    return (m, nnm, index)
