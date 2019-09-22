#!/usr/bin/python3

import sys
import argparse

from Transposer.interp import *
from Transposer.element import Element
from Transposer.tier import *
from Transposer.flanker import *

# need to get rid of the args and change to parameters so can run those
# from the Glycine Remap module
# TODO: Refactor other files in case they use anything from arg parser
# remap must be runnable as a function with no args only parameters
# or could create list of parser objects

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
