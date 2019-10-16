# temp file for backmapping methods
# stuff in here is likely to move

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
    rf = search_BDB(element.endLocation + 1, element.endLocation + (n+1))

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
