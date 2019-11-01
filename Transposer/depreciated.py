
def write_fasta(sels, output):
    with open(output, 'w') as out:
        for el in sels:
            print(el.status, 'type at write')
            out.write(el.get_header() + '\n' + el.seq + '\n')


def write_csv(sels, output):
    with open(output, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(['Name', 'Accession', 'Chr', 'Start', 'Length',
                         'Seq', 'Left Flank', 'Right Flank'])
        for el in sels:
            writer.writerow(el.get_row())



def make_jobs(self):
    '''
    generates a list of bowtie commands to be run based on the number of
    consensus clusters in the run object. Each consensus cluster will be
    bowtied and then needs to tripped of duplicate elements. Solo and intact
    will need to be seperated. Will be a list with first list being intact
    jobs and second list being solo jobs. Moved after created make jobs two.
    '''
    sam_dir = self.write_dirs[1]
    jobs = []
    ave_len = 0
    for c in self.cie_cons:
        sam_name = os.path.basename(c).split('.')[0] + '_intact.sam'
        sam_file = os.path.join(sam_dir, sam_name)
        l = get_intact_length(c)
        # redundancy here think about a way to reduceonly matches the parent
        # dont reall need to calculate for the intact elements
        jobs.append(Search(BTI=self.BTI, con_file=c,
                           out_file=sam_file, num_old_els=self.num_cie,
                           type='I', acc=self.cur_acc, BDB=self.cur_BDB,
                           intact_len=l))
        ave_len += l
    ave_len = round(ave_len / len(self.cie_cons))
    # TODO: edit the search objects so taking in all required parameters
    # calculate average length of all consensuses and use that for
    # calculating solo elements need to round to whole number
    if self.csi_cons is not None:
        for c in self.csi_cons:
            sam_name = os.path.basename(c).split('.')[0] + '_solo.sam'
            sam_file = os.path.join(sam_dir, sam_name)
            jobs.append(Search(BTI=self.BTI, con_file=c,
                               out_file=sam_file, num_old_els=self.num_csi,
                               type='S', acc=self.cur_acc, BDB=self.cur_BDB,
                               intact_len=ave_len))
        # need numbers of both the intact and solo files
    self.jobs = jobs  # jobs now stored in the run object



def prune_clstr_cons(self):
    cat_cie = os.path.join(self.write_dirs[0], 'intact_cons')
    cie_clstr = os.path.join(self.write_dirs[0], 'in_con_clstr')
    cat_csi = os.path.join(self.write_dirs[0], 'solo_cons')
    csi_clstr = os.path.join(self.write_dirs[0], 'solo_con_clstr')
    # make file names for concat consensus files

    c = ['cat']
    try:
        # make concat files of all consensuses
        cmd_cie = c + self.cie_cons + ['>', cat_cie]
        cmd_csi = c + self.cie_cons + ['>', cat_csi]
        os.system(' '.join(cmd_cie))
        os.system(' '.join(cmd_csi))
    except subprocess.CalledProcessError as e:
        return 1
        print('++++++++++ subprocess con cat error!!!')
    print(cat_cie)
    run_cd_hit(output=cie_clstr, input_file=cat_cie + '.clstr')
    run_cd_hit(output=csi_clstr, input_file=cat_csi + '.clstr')
    # get lists of single consensus clusters
    # need to actualy make the cluster then call
    # clstr file on the cluster file not just the fasta
    # last error was becasue never rand cd hit
    print(cie_clstr)
    cie_con_files = ClstrFile(cie_clstr + '.clstr').get_singles()
    csi_con_files = ClstrFile(csi_clstr + '.clstr').get_singles()
    # change con variables to reflect the clustering results
    # remove the first character as it will be a > due to writing
    # the clstr file
    self.cie_cons = [os.path.join(self.write_dirs[0], x[1:]) for x in cie_con_files]
    self.csi_cons = [os.path.join(self.write_dirs[0], x[1:]) for x in csi_con_files]


def search(con, index, outputFile, verbose=False):
    '''
    uses a consensus sequence and bowtie2 index to preform a global allignment
    to the given index, returns results as a sam file

    From tier.py and replaced with the search object bowtie funcitonality
    '''

    try:
        commandSearch = "bowtie2 -x {} -f {} -a --non-deterministic -S {}".format(index, con, outputFile)
        print(commandSearch)
        samCommand = "Samtools view -S -b {} ".format(outputFile)
        os.system(commandSearch)
        outputFile = toSortedBam(outputFile, verbose)

        return outputFile
        # will run the command and return a samfile
    except FileNotFoundError:
        print("FileNotFoundError at search in search.py")

import csv
import os
import subprocess
from Transposer.element import Element
# want this to work with eleemnt objects


def cigarParser(cigar):
    '''
    reads the cigar info from an Bowtie allignment and translates
    into the length on the alligned reference
    '''
    length = 0
    temp = ""
    for char in cigar:
        temp = temp + "" + char
        if char == "M" or char == "D":
            temp = temp[:-1]
            length += int(temp)
            temp = ""
        elif char == "I" or char == "H":
            temp = ""

    return length


def process_rname(r_name):
    return r_name.split('|')[3]


def get_chr_num(acc_path, acc):
    try:
        lines = []
        with open(acc_path) as names:
            for line in names:
                l = line.strip().split('\t')
            if l[1] == acc:
                return int(l[0])
    except FileNotFoundError as e:
     return -1

def sort_elements(elements):
    chr_dict = {}
    sorted_list = []
    for e in elements:
        if e.chr in chr_dict:
            chr_dict[e.chr].append(e)
        else:
            chr_dict[e.chr] = [e]

    for chr in sorted(chr_dict):
        sorted_list += sorted(chr_dict[chr], key= lambda e: (e.startLocation))

    return sorted_list


class Sam():
    def __init__(self, path, bdb, accession_path, type, element_set=set([])):
        self.path = path
        self.acc_path = accession_path
        self.name = os.path.basename(path)
        self.element_set = element_set
        self.dbd = Blast_DB(bdb, accession_path)
        self.type = type

    def make_element_set(self):
        try:
            element_set = set([])
            with open(self.path) as path:
                reader = csv.reader(path, delimiter='\t')
                for row in reader:
                    if row[0][0] != '@':
                        name = row[0]
                        acc = process_rname(row[2])
                        chr = get_chr_num(self.acc_path, acc)
                        start = int(row[3])  # 1 based start
                        length = cigarParser(row[5])
                        end = start + length
                        sequence = self.dbd.search(start, end, acc)
                        element_set.add(
                            Element(name, acc, chr, start, end, length, type, sequence))
            self.element_set = element_set
            return 0
        except FileNotFoundError as e:
            return set([])

    def remove_dups(self):
        '''
        Removes elements from the element set that are duplicates within an individual
        sam file. Additional processing required at the run level to seperate
        elements that are duplicate between sam files.
        '''
        id_set = set([])
        nr_list = set([])
        for e in self.element_set:
            id = str(e.chr) + str(e.startLocation)
            if id not in id_set:
                nr_list.add(e)
                id_set.add(id)

        self.element_set = nr_list


    def type_elements(self, len_intact, allowance=25):
        '''
        If Sam is solo type then sort solo elements and remove those that are
        actually the ends of an LTR. If sam is an intact type then change the
        status of all elements in the set to intact.
        '''
        if self.type == 'S':
            solo_set = set([])
            sort_e = sort_elements(self.element_set)
            # sort solos by start location
            i = 0
            n = len(sort_e)-1
            print(n, 'number of element')
            while i < n:
                print(i)
                current = sort_e[i]
                next = sort_e[i+1]  # watch out of bounds error
                if current.chr != next.chr:
                    solo_set.add(current)
                    i+=1
                elif next.startLocation - len_intact - allowance <= current.endLocation:
                    # elements are close enough to be intact
                    i += 2
                else:
                    current.status = 'S'
                    solo_set.add(current)
                    i += 1
                if i == n:  # last element if only one must be a solo
                    sort_e[n].status = 'S'
                    solo_set.add(sort_e[n])
                    break
            self.element_set = solo_set
        else:
            for e in self.element_set:
                e.status = 'I'


class Blast_DB():
    def __init__(self, path, accession_path):
        self.path = path
        self.acc = accession_path


    def search(self, start, end, entry, r_seq=True):
        '''
        search the blast db for a sequence in an entry
        '''

        seq_cmd = 'blastdbcmd -db {} -dbtype nucl -range {}-{} -entry {}'.format(
                   self.path, start, end, entry)
        try:
            output = subprocess.check_output(seq_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            return ''
        if r_seq:
            return ''.join(str(output).split('\\n')[1:])
            # returns string of just the sequence
        else:
            return output  # else return all output


# from interp.py

def toSortedBam(samFile, verbose):
    # # TODO: figure out if actually need this one
    # might be good to use for verify consensus methods
    '''
    reads sam file and converts to a sorted bam
    then removes redundant files
    '''
    try:
        outputFile = samFile
        if 'sam' not in outputFile:
            outputFile = outputFile + ".sam"

        bam = samFile.replace(".sam", ".bam")
        sortedBam = samFile.replace(".sam", ".sorted.bam")

        # commands to be executed
        commandBam = "samtools view -S -b {} > {}".format(samFile, bam)
        commandSort = "samtools sort {} -o {}".format(bam, sortedBam)

        if verbose:
            print("Running command: " + commandBam)
        os.system(commandBam)
        if verbose:
            print("Running command: " + commandSort)
        os.system(commandSort)

        try:
            # removes redundant files but leaves the sorted bam for user
            os.system("rm {}".format(samFile))
            os.system("rm {}".format(bam))

        except FileNotFoundError:
            print("FileNotFoundError at toSortedBam when trying to remove files")

        return sortedBam   # returns a sorted bam file name corresponding to one created in method

    except FileNotFoundError:
        print("FileNotFoundError at toSortedBam in search.py")


def toTxt(samFile, verbose):
    '''
    takes sam file and converts to a txt file to easier parsing
    '''
    outputFile = samFile + ".txt"

    try:
        command = "samtools view {} > {}".format(samFile, outputFile)
        if verbose:
            print("Running Command: " + command)
        os.system(command)
        return outputFile   # returns txt file name corresponding to txt file created

    except FileNotFoundError:
        print("FileNotFoundError at toTxt in search.py")


def createAllignmentList(bamToTxtFile, dict, verbose, curBlastDB):
    '''
    takes the converted txt (sam -> txt) file and returns a list of
    Element objects
    '''
    try:
        elementList = []
        with open(bamToTxtFile, "r") as txt:

            for line in txt:
                line = line.split("\t")  # splits each line into a list
                print(line)
                name = line[2]
                start = line[3]
                length = line[5]  # gives cigar to the length as a temp holder

                elementList.append(
                    Element(name, name, start, 0, length, "NONE", "ATGC"))

        for element in elementList:
            print(element.name)
            print(element.accession)
            print('--------------------------')
            element.length = cigarParser(element.length)
            element.endLocation = element.startLocation + element.length - \
                1  # changed calculation of length using CIGAR
            element.status = "INTACT"  # defualt status is INTACT
            element.seq = getElementSeq(curBlastDB, element)

            # to get the element sequence need a blast DB of the new assembly and assencion number

            try:
                # uses provided dictionary to set name to chr number
                element.name = dict[element.name]

            except KeyError:
                if verbose:
                    print("keyError at createAllignmentList")
                    print("key used was " + element.name)

        # removes txt version as no longer used
        os.system("rm {}".format(bamToTxtFile))

        return elementList

    except FileNotFoundError:
        print("FileNotFoundError at makeAllignDict in samParser.py")


def translateName(assenstionNums):
    '''
    reads GenBank accenstion number file and returns a dictionary
    which can be used to translate chr numbers
    '''
    try:
        chrs = {}
        with open(assenstionNums) as names:
            for line in names:
                line = line.strip()
                line = line.split("\t")
                number, CM = line
                chrs[CM] = number
        return chrs

    except FileNotFoundError:
        print("FileNotFoundError" + " at translateName")


def mergeLists(soloList, conList, familyName):
    '''
    merges the solo and intact element lists
    renames elements based on new relative positions
    and returns as a list of Elements
    '''
    # possible check for None type for either list
    merge = soloList + conList
    mergeDict = {}
    finalList = []

    for element in merge:
        if element.name not in mergeDict:
            mergeDict[element.name] = [element]
        else:
            mergeDict[element.name].append(element)

    for key in sorted(mergeDict.keys()):
        mergeDict[key] = sorted(mergeDict[key], key=lambda e: e.startLocation)
        for i, element in enumerate(mergeDict[key], 1):
            element.name = str(familyName) + " " + str(key) + "-" + str(i)
        finalList += mergeDict[key]

    return finalList


def findSolos(LTRList, completeCon, allowance):
    '''
    takes list of LTR elements and determines if each element is a solo element
    based on a given allowance in base pairs and the complete element length
    '''

    soloList = []
    i = 0
    completeConLength = 0

    try:
        with open(completeCon) as con:
            completeConLength = len(con.readline())
    except FileNotFoundError:
        return "FileNotFoundError at findSolos opening completeCon file"

    while i < len(LTRList):

        reach = LTRList[i].endLocation + completeConLength
        diff = (LTRList[i + 1].startLocation - reach)

        if(diff > -(completeConLength) and diff <= 0):  # test for truncated element
            LTRList[i].status = "INTACT"
            LTRList[i + 1].status = "INTACT"
            i += 1
            if i == len(LTRList) - 1:
                break

        elif (abs(diff) > allowance):
            LTRList[i].status = "SOLO"
            soloList.append(LTRList[i])
            if i == len(LTRList) - 2:
                soloList.append(LTRList[i + 1])
                break
        else:
            LTRList[i].status = "INTACT"
            LTRList[i + 1].status = "INTACT"
            i += 1
            if i == len(LTRList) - 1:
                break
        i += 1

    return soloList


def cigarParser(cigar):
    '''
    reads the cigar info from an Bowtie allignment and translates
    into the length on the alligned reference
    '''
    length = 0
    temp = ""
    for char in cigar:
        temp = temp + "" + char
        if char == "M" or char == "D":
            temp = temp[:-1]
            length += int(temp)
            temp = ""
        elif char == "I" or char == "H":
            temp = ""

    return length


# from tier.py
def score_consensus(con, old_index, old_elements, verbose=False):
    TEMP_FILE = os.path.join(os.getcwd(), 'temp.fasta')
    '''
    This function will attempt to validate that a consensus sequence will return
    good results when looking for elements in the new assembly by searching
    the consensus against the old assembly and comparing the elements returned
    against the old elements. The assumption being that if it can yeild the
    old elements it should be able to find the new elements.
    '''
    bowtie_results = toSortedBam(search(con, old_index, TEMP_FILE), verbose)
    # need to change the toTxt becuase output will not be in a fasta format so
    # it will be better to read number of elements using a subprocess command
    allign_cmd = ['ssamtools', 'view', bowtie_results, '|', 'wc', '-l']
    num_alligns = subprocess.check_output(allign_cmd)
    num_old_elements = len(ft.read_as_tuples(old_elements))

    return len(elements / old_elements)

# from remap.py

def transpose(BT_index, solo_con, intact_con, cur_acc, output, name,
              cur_BDB, old_BDB, old_acc, old_elements, allowance,
              verbose):

    sortedTxtCon, ElementListLTR = [], []
    ElementListLTR, ElementListCon = [], []
    finalList = []

    temp_dir = os.getcwd()
    temp_solo = os.path.join(temp_dir, 'temp_solo.sam')
    temp_intact = os.path.join(temp_dir, 'temp_intact.sam')

    chrDict = translateName(cur_acc)
    # need to check if either of the consensus are none and avoid
    # running if that is the case and deal with consequences later in
    # the run
    if solo_con is not None:
        # print("Searching for solo elements")
        sortedTxtLTR = toTxt(search(solo_con, BT_index,
                                    temp_solo, verbose), verbose)
        ElementListLTR = createAllignmentList(
            sortedTxtLTR, chrDict, verbose, cur_BDB)

        soloList = findSolos(ElementListLTR, intact_con, allowance)

    if intact_con is not None:
        #print("Searching for intact elements")
        sortedTxtCon = toTxt(search(intact_con, BT_index,
                                    temp_intact, verbose), verbose)
        ElementListCon = createAllignmentList(
            sortedTxtCon, chrDict, verbose, cur_BDB)
    # point where None needs consideration

    # ElementListLTR = createAllignmentList(
        # sortedTxtLTR, chrDict, verbose, cur_BDB)
    # ElementListCon = createAllignmentList(
    #    sortedTxtCon, chrDict, verbose, cur_BDB)

    #soloList = findSolos(ElementListLTR, intact_con, allowance)

    # if only one list cannot merge list so may need another sorting method

    # need to check if merged list will work with an empty list is True then
    # will not need to do any more if statements for testing

    finalList = mergeLists(soloList, ElementListCon, name)
    # final list is given to the backmapping part so final list needs to come
    # out regardless of what is included either LTR or solo only
    if old_BDB and old_elements:
        # read the previous elements into a list
        backDict = backMapElements(readPreviousElements(old_elements, cur_accPrev),
                                   finalList, old_BDB, cur_BDB, cur_accPrev, cur_acc)
        matchsToTxt(backDict)

    if output == True:
        for element in finalList:
            print(element.toStringFasta())
    else:
        with open(output, "w") as out:
            for element in finalList:
                out.write(element.toStringFasta)
