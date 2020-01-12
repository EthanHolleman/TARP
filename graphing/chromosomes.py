import csv

def read_matches(matches, delim=','):
    with open(matches) as m:
        reader = csv.reader(matches, delimiter=delim)
        return reader.readlines()

# need to get in chromosome lengths

def bin_sort_matches(matches):
    chr_bins, cur_chr = [[]]*num_chr, None
    for m in matches:
        if m[1] != None:
            if cur_chr != int(m[1]) and int(m[1]) <= num_chr and int(m[1]) > 0:
                cur_chr = int(m[1])
                chr_bins[cur_chr-1].append(m)

    for bin, _ in enumerate(chr_bins):
        chr_bins[bin] = sorted(chr_bins[bin], key=lambda x: x[3])

    return chr_bins

def key(item):
    if isinstance(item, feature):
        return item.position
    else:
        return item[3]


def merge_feats_matches(feats_bins, match_bins, key):
    merge_bins = [[]] * len(feats_bins)  # should be same as match bins
    for f, m in zip(feats_bins, match_bins):
        temp_list = f + m
        sorted_list = sorted(temp_list, key=key)
        merge_bins.append(sorted_list)

    return merge_bins

def ideo_line_writer(chr, start, end, stain):
    return 'chr{}\t{}\t{}\t{}\t{}'.format(chr, start, 'NA', end, stain)


def make_ideo(vcf_file, match_file, output, key):
    feats = vcf_reader(vcf_file)
    matches = bin_sort_matches(read_matches(match_file))
    big_list = merge_feats_matches(feats, matches, key)
    with open(output, 'w') as out:
        for chr in big_list:
            cur_el = chr.pop()
            while chr:
