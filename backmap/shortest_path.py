from backmap.anchored_remap import feature_sort
from Transposer.element import Element
from Transposer
new = []


def flatten(l):
    return [item for sublist in l for item in sublist]


def run_min_distance(new_elements, old_elements):
    matched, unmatched = [], []
    new_element_index = feature_sort(
        new_elements, key=lambda x: x.startLocation)
    old_element_index = feature_sort(
        old_elementsld, key=lambda x: x.startLocation)

    for i in range(len(old_elements_index)):  # len index dependent on num chr
        m, um = match_chromosome_elements(
            new_element_index[i], old_elements_index[i])
        matched.append(m)
        unmatched.append(um)

    return flatten(matched), flatten(unmatched)


def make_paired_index(new_element_list, index):
    '''
    Returns zipped iterator of lists representing the elements in each chromsome found
    by TARP (new) at index 0 of each tuple and the old elements at index 1.
    '''
    chr_bins = feature_sort(element_list, key=lambda x: x.startLocation)
    return zip(chr_bins, index)


def match_chromosome_elements(new_elements, old_elements):
    '''
    Given list of old elements in a chromsome, seperated them by type and
    then matches each type using shortest distance method. New and old elements
    must have the same types. (statuses)
    '''
    statuses, matches, unmatches = get_statuses(
        new_elements), [], []  # should be same as old
    for s in statuses:
        n, o = list(filter(lambda x: x.status, new_elements)), \
            list(filter(lambda x: x.status, old_elements))
        cs, short, long = minimize_distance(n, o)  # find min distnace allign
        matches += pair_matches(cs, short, long)
        unmatches += return_unmatched(cs, short, long)

    return matches, unmatches


def get_statuses(elements):
    statuses = set([])
    for el in elements:
        if el.status not in statuses:
            statuses.add(el.status)
    return statuses


def seperate_types(element_list):
    types_dict = {}
    for element in element_list:
        if element.status in types_dict:
            types_dict[element.status].append(element)
        else:
            types_dict[element.status] = [element]

    return types_dict


def pair_types(dict_one, dict_two):
    one_keys = list(dict_one.keys())
    return (dict_one[one_keys[0]], dict_two[one_keys[0]],
            dict_one[one_keys[1]], dict_two[one_keys[1]])


def minimize_distance(list_one, list_two):
    if len(list_one) >= len(list_two):
        longer_list = list_one  # if the same doesnt matter which is longer list
        shorter_list = list_two
    else:
        longer_list = list_two
        shorter_list = list_one

    shorter_list = sorted(shorter_list, key=lambda x: x.startLocation)
    longer_list = sorted(longer_list, key=lambda x: x.startLocation)

    d, cs = len(longer_list) - len(shorter_list), (float('inf'), 0)
    for i in range(0, d + 1):
        s = sum(abs(x.startLocation - y.startLocation)
                for x, y in zip(shorter_list, longer_list[i:len(shorter_list) + i]))
        if s < cs[0]:
            cs = (s, i)  # holds distance diff and i value of slice

    return (cs[1], shorter_list, longer_list)


def pair_matches(cs, shorter_list, longer_list):
    return [(a, b) for a, b in
            zip(shorter_list, longer_list[cs:len(shorter_list) + 1])]


def return_unmatched(cs, shorter_list, longer_list):
    return longer_list[:cs] + longer_list[len(shorter_list) + 1:]
