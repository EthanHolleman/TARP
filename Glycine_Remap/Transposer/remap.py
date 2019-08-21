#!/usr/bin/python
import sys
import argparse
from interp import *
from argsRemap import *
from element import Element
from tier import *
from flanker import *


def main():
    sortedTxtCon, ElementListLTR = [], []
    ElementListLTR, ElementListCon = [], []
    finalList = []

    args = argsRemap()
    chrDict = translateName(args.chrKeys)
    # need to check if either of the consensus are none and avoid
    # running if that is the case and deal with consequences later in
    # the run
    if args.LTRcon is not None:
        #print("Searching for solo elements")
        sortedTxtLTR = toTxt(search(args.LTRcon, args.index,
                                    "TestRunLTR.sam", args.verbose), args.verbose)
        ElementListLTR = createAllignmentList(
            sortedTxtLTR, chrDict, args.verbose, args.curBlastDB)

        soloList = findSolos(ElementListLTR, args.seqCon, args.allowance)

    if args.seqCon is not None:
        #print("Searching for intact elements")
        sortedTxtCon = toTxt(search(args.seqCon, args.index,
                                    "TestRunCon.sam", args.verbose), args.verbose)
        ElementListCon = createAllignmentList(
            sortedTxtCon, chrDict, args.verbose, args.curBlastDB)
    # point where None needs consideration

    #ElementListLTR = createAllignmentList(
        #sortedTxtLTR, chrDict, args.verbose, args.curBlastDB)
    #ElementListCon = createAllignmentList(
    #    sortedTxtCon, chrDict, args.verbose, args.curBlastDB)

    #soloList = findSolos(ElementListLTR, args.seqCon, args.allowance)

    # if only one list cannot merge list so may need another sorting method

    # need to check if merged list will work with an empty list is True then
    # will not need to do any more if statements for testing
    
    finalList = mergeLists(soloList, ElementListCon, args.name)
    #final list is given to the backmapping part so final list needs to come
    # out regardless of what is included either LTR or solo only
    if args.prevBlastDB and args.prevElements:
        # read the previous elements into a list
        backDict = backMapElements(readPreviousElements(args.prevElements, args.chrKeysPrev),
                                   finalList, args.prevBlastDB, args.curBlastDB, args.chrKeysPrev, args.chrKeys)
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
