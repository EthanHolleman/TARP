#main file for integrating the transposer into the other software

def runTransposer:
    '''
    runs the transposer program for all families in the dictionary given
    requires: dictionary, the accession numbers/files, blastDB and bowtie index
    for the whole genome (only need 1 or each), both family consensus files
    need a way to get the file names for each family 
    '''
    pass

def createCommand:
    #creates the command to be run by transposer
    pass

def createR:
    #takes output files from transposer run and makes data file for easy graphing in r
    pass

def makeLocationPlot:
    #uses createR data and calls to R to make plots
    #need to figure out if can save the plots
    #could also generate R code to create each of the plots in markdown style
    pass
