import itertools
import csv

def index_teardown(index):
    '''
    Takes in the index used for ham searches after elements have been
    remapped. Rum ham will pop elements that are matched out of the index.
    Index teardown pulls the remaining elements out and reformats them into
    single level list.
    '''
    return (list(itertools.chain.from_iterable(index)))

def match_formater(matches):
    '''
    Takes in match tuples and returns genorator that is ready to be written by
    match writer function.
    '''
    for a, b in matches:
        x = [a.name, a.chr, a.accession, a.startLocation, a.length, a.left, a.right, a.seq]
        y = [b.name, b.chr, b.accession, b.startLocation, b.length, b.left, b.right, b.seq]
        yield x+y

def match_writer(gen_matches, output):
    HEAD = ['Name', 'Chromosome', 'Accession', 'Start', 'Length', 'Left Flank',
            'Right Flank', 'Sequence'] * 2
    print(HEAD)
    try:
        with open(output, 'w') as out:
            write = csv.writer(out, delimiter=',')
            write.writerow(HEAD)
            for m in gen_matches:
                write.writerow(m)
                print(m, 'writing now')
    except (FileNotFoundError, IsADirectoryError) as e:
        print('{} location not found or is a dir'.format(output))

def nomatch_writer(nm, output):
    pass
    '''
    Writes elements that had no match to a csv file. Should be able to add a
    column that tells if from the search or the old elements.
    '''
