
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
