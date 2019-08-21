import sys, argparse, os

def argsRemap():
    #control structure for remap functions
    #add ability to give the assenstion file which will translate the chromosome number
    ### TODO: add parameter and control structure to do create index
    parser = argparse.ArgumentParser()
    parser.add_argument("-index", "-i", help = "Set Bowtie index to be used")
    parser.add_argument("-LTRcon","-l", defualt=None, help = "LTR consensus sequence")
    parser.add_argument("-seqCon", "-s",defualt=None, help = "Complete element consensus")
    parser.add_argument("-allowance", "-a", type = int,
    help = "Number base pair variance in LTR map", nargs='?',const=150)
    parser.add_argument("-name", "-n", help = "Name of element family")
    parser.add_argument("-newIndex", "-b", help= "Tells program to create new bowtie index")
    parser.add_argument("-chrKeys", "-k", help = "Acts as a key to translate NCBI assention file chromosome names to numbers")
    parser.add_argument("-outputFile", "-o", help = "Specify name of output file to write results to, will be in fasta style")
    parser.add_argument("-verbose", "-v", help = "Transposer is quiet by defualt, typing True will set program to verbose")
    parser.add_argument("-prevBlastDB", "-p", help = "Nuc BLAST database created from the same reference as the elements to be remapped")
    parser.add_argument("-curBlastDB", "-c", help = "Nuc BLAST database created from the most recent reference")
    parser.add_argument("-prevElements", "-e", help = "Nuc BLAST database created from the most recent reference")
    parser.add_argument("-chrKeysPrev", "-m", help = "Nuc BLAST database created from the most recent reference")

    args = parser.parse_args()

    if not args.index:
        print("Need index")
        sys.exit()
    if not args.chrKeys:
        print("Please provide a file of keys for chromosomes")
        sys.exit()
    if not args.outputFile:
        args.outputFile = True
    if not args.curBlastDB:
        print("Please provide a blastdb nucl type of the most current assembly")
        sys.exit()
    if not args.prevBlastDB or not args.prevElements:
        print("No previous assembly and or previous element list provided so new elements will not be backmapped")



    return args
