import subprocess,random

def seperateTypes(element_list):
    #creates dictionary of elements, keys are element status
    INTACT = 'INTACT'
    SOLO = 'SOLO'
    pass
    status_dict = {INTACT : [], SOLO : []}
    for element in element_list:
        try:
            status_dict[element.status].append(element)
        except KeyError:
            print('Key error at seperateTypes: ' + element.status)

    return status_dict

def getRandomElements(elements, fam):
    #selects 10 random elements of the family to make fasta file of
    #list is passed to createConsensusFile
    #elements should be of one family
    rand_elements = []
    for i in range(0,11):
        rand = random.randint(0,len(elements)-1)
        rand_elements.append(elements[rand])

    return rand_elements

def clustalize(rep_elements):
    #takes in rep elements of a single type not whole family
    #makes the temp fasta input file and does the good good

    status = element_list[0].status #all elements should be of same status
    input_name = '{}_{}_repElements.fasta'.format(fam,status)
    output_name = '{}_{}_consensus.fasta'.format(fam,status)

    with open(input_name, 'w') as fasta:
        for element in rep_elements:
            fasta.write(element.toStringFasta())

    clustal_command = 'clustalo -i {} -o {} -v'.format(input_name,output_name)

    subprocess.call(clustal_command)
    #subprocess.call('rm input_name')


def makeFamConsensus(element_list):
    #runs other clustal methods to create family solo and intact consensuses
    INTACT = 'INTACT'
    SOLO = 'SOLO'

    types = seperateTypes(element_list)
    intacts = types[INTACT]
    solos = types[SOLO]

    rep_intacts = getRandomElements(intacts)
    rep_solos = getRandomElements(solos)

    clustalize(rep_intacts)
    clustalize(rep_solos)
