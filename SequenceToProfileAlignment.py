#!/usr/bin/env python3.7

import argparse
import numpy as np

def getSeqList():

    seqList = list()

    with open(alignedSequencesFileName, "r") as inFile:
        for line in inFile.readlines():
            seqName, seq = line.strip('\n').split()
            seqList.append(seq)

        inFile.close()

    return seqList

def getProfileScore(seqList):

    noSeq = len(seqList)
    lenSeq = len(seqList[0])

    profileScore = np.zeros((5, lenSeq), dtype='float')

    for profileIdx in range(lenSeq):

        noGap, noA, noC, noG, noT = 0.0, 0.0, 0.0, 0.0, 0.0

        for seqIdx in range(noSeq):

            if seqList[seqIdx][profileIdx] == '-':
                noGap += 1.0
            elif seqList[seqIdx][profileIdx] == 'A':
                noA += 1.0
            elif seqList[seqIdx][profileIdx] == 'C':
                noC += 1.0
            elif seqList[seqIdx][profileIdx] == 'G':
                noG += 1.0
            elif seqList[seqIdx][profileIdx] == 'T':
                noT += 1.0

        profileScore[0][profileIdx] = noGap / noSeq
        profileScore[1][profileIdx] = noA / noSeq
        profileScore[2][profileIdx] = noC / noSeq
        profileScore[3][profileIdx] = noG / noSeq
        profileScore[4][profileIdx] = noT / noSeq

    return profileScore

def getProfile(profileScore):

    profile = ""
    for profileIdx in range(len(profileScore[0])):

        currProfile = ""
        maxScore = 0.0

        if profileScore[0][profileIdx] > maxScore:
            maxScore = profileScore[0][profileIdx]
            currProfile = "-"
        if profileScore[1][profileIdx] > maxScore:
            maxScore = profileScore[1][profileIdx]
            currProfile = "A"
        if profileScore[2][profileIdx] > maxScore:
            maxScore = profileScore[2][profileIdx]
            currProfile = "C"
        if profileScore[3][profileIdx] > maxScore:
            maxScore = profileScore[3][profileIdx]
            currProfile = "G"
        if profileScore[4][profileIdx] > maxScore:
            currProfile = "T"

        profile += currProfile

    return profile

def getSeqFasta():
    return open(seqFastaFileName,"r").readlines()[1]

def getScore(currChar, idx):

    score = 0.0
    charList = ["-", "A", "C", "G", "T"]

    for charIdx in range(5):
        if currChar == "-" or charIdx == 0:
            alignmentScore = gapPenalty
        elif currChar == charList[charIdx]:
            alignmentScore = matchScore
        else:
            alignmentScore = mismatchPenaltyScore
        score = score + (alignmentScore * profileScore[charIdx][idx-1])

    return score

def NaiveNeedlemanWunsch():

    lenSeq = len(seqFasta) + 1
    lenProfile = len(profileFasta) + 1

    alignmentMatrix = np.zeros((lenSeq, lenProfile), dtype = 'float')
    tracebackMatrix = np.zeros((lenSeq, lenProfile), dtype='str')

    for profileIdx in range(lenProfile):
        alignmentMatrix[0][profileIdx] = profileIdx * gapPenalty
        tracebackMatrix[0][profileIdx] = '1' # Left

    for seqIdx in range(lenSeq):
        alignmentMatrix[seqIdx][0] = seqIdx * gapPenalty
        tracebackMatrix[seqIdx][0] = '2' # Up

    alignmentMatrix[0][0] = 0.0

    for seqIdx in range(1, lenSeq):
        for profileIdx in range(1, lenProfile):

            diagonalValue = alignmentMatrix[seqIdx - 1][profileIdx - 1] + getScore(seqFasta[seqIdx - 1], profileIdx)
            upValue = alignmentMatrix[seqIdx - 1][profileIdx] + gapPenalty
            leftValue = alignmentMatrix[seqIdx][profileIdx - 1] + getScore('-', profileIdx)

            maxValue = max(diagonalValue, upValue, leftValue)

            alignmentMatrix[seqIdx][profileIdx] = maxValue

            if maxValue == diagonalValue:
                tracebackMatrix[seqIdx][profileIdx] = '0'

            elif maxValue == leftValue:
                tracebackMatrix[seqIdx][profileIdx] = '1'

            elif maxValue == upValue:
                tracebackMatrix[seqIdx][profileIdx] = '2'

    return alignmentMatrix, tracebackMatrix

def TraceBack():

    alignedSequence = ""
    alignedProfile = ""

    rowIdx = alignmentMatrix.shape[0] - 1
    colIdx = alignmentMatrix.shape[1] - 1

    while rowIdx > 0 and colIdx > 0:

            if tracebackMatrix[rowIdx][colIdx] == '0':
                alignedSequence += seqFasta[rowIdx - 1]
                alignedProfile += profileFasta[colIdx - 1]
                rowIdx -= 1
                colIdx -= 1

            elif tracebackMatrix[rowIdx][colIdx] == '1':
                alignedSequence += "-"
                alignedProfile += profileFasta[colIdx - 1]
                colIdx -= 1

            elif tracebackMatrix[rowIdx][colIdx] == '2':
                alignedSequence += seqFasta[rowIdx - 1]
                alignedProfile += "-"
                rowIdx -= 1

    while rowIdx > 0:
        alignedSequence += seqFasta[rowIdx - 1]
        rowIdx -= 1

    while colIdx > 0:
        alignedProfile += profileFasta[colIdx - 1]
        colIdx -= 1

    return alignedSequence[::-1], alignedProfile[::-1]

def createNewAlignedSequenceFile():

    with open(allAlignedSequencesFileName, "w") as outFile:

        with open(alignedSequencesFileName, "r") as inFile:

            for line in inFile.readlines():
                outFile.write(line)

            inFile.close()

        outFile.write("\nsequence\t" + alignedSequence)

        outFile.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fasta", help = "fasta file to be aligned to the given profile")
    parser.add_argument("--aln", help = "aln-formatted file containing all given alignments")
    parser.add_argument("--out", help = "aln-formatted file containing all given alignments including the newly aligned sequence")
    parser.add_argument("--gap", help = "gap penalty score")
    parser.add_argument("--mismatch", help = "mismatch penalty score")
    parser.add_argument("--match", help = "matching score")

    args = parser.parse_args()

    seqFastaFileName = str(args.fasta)
    alignedSequencesFileName = str(args.aln)
    allAlignedSequencesFileName = str(args.out)
    gapPenalty = float(args.gap)
    mismatchPenaltyScore = float(args.mismatch)
    matchScore = float(args.match)

    profileScore = getProfileScore(getSeqList())
    profileFasta = getProfile(profileScore)
    seqFasta = getSeqFasta()

    alignmentMatrix, tracebackMatrix = NaiveNeedlemanWunsch()
    alignedSequence, alignedProfile = TraceBack()

    createNewAlignedSequenceFile()