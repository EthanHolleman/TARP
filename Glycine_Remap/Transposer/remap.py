#!/usr/bin/python
import sys, argparse
from interp import *
from argsRemap import *
from element import Element
from tier import *
from flanker import *

def main():

    args = argsRemap()
    chrDict = translateName(args.chrKeys)
    print("Searching for solo elements")
    sortedTxtLTR = toTxt(search(args.LTRcon, args.index, "TestRunLTR.sam", args.verbose ),args.verbose)
    print("Searching for intact elements")
    sortedTxtCon = toTxt(search(args.seqCon, args.index, "TestRunCon.sam",args.verbose ),args.verbose)


    ElementListLTR = createAllignmentList(sortedTxtLTR, chrDict, args.verbose, args.curBlastDB)
    ElementListCon = createAllignmentList(sortedTxtCon, chrDict, args.verbose, args.curBlastDB)


    soloList = findSolos(ElementListLTR, args.seqCon, args.allowance)
    finalList = mergeLists(soloList, ElementListCon, args.name)

    if args.prevBlastDB and args.prevElements:
        #read the previous elements into a list
        backDict = backMapElements(readPreviousElements(args.prevElements, args.chrKeysPrev),finalList,args.prevBlastDB, args.curBlastDB,args.chrKeysPrev,args.chrKeys)
        matchsToTxt(backDict)

    if args.outputFile == True:
        for element in finalList:
            print(element.toStringFasta())
    else:
        with open(args.outputFile, "w") as out:
            for element in finalList:
                out.write(element.toStringFasta)


if __name__ == "__main__":
    main()
