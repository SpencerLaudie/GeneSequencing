#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean
    # that tells you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        score, P = self.getEditDistance(seq1, seq2)
        alignment1, alignment2 = self.aligned(seq1, seq2, P)

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

    # takes two sequences as arguments
    # returns banded edit distance if self.banded is true; otherwise, returns unrestricted edit distance
    def getEditDistance(self, seq1, seq2):
        if self.banded:
            return self.getEditDistanceBanded(seq1, seq2)
        else:
            return self.getEditDistanceUnrestricted(seq1, seq2)

    #  Edit Distance Algorithm (Banded Implementation)
    #  takes two sequences as arguments
    #  returns banded edit distance and 2D array of back-pointers
    #  t: O(kn)  outer loop of double for-loop iterates n times. inner for-loop iterates k times.
    #  s: O(kn)  two 2-dimensional arrays (E & P) of size n x k are created.
    def getEditDistanceBanded(self, seq1, seq2):
        seq1 = seq1[:self.MaxCharactersToAlign]
        seq2 = seq2[:self.MaxCharactersToAlign]
        len1 = len(seq1) + 1  # n = len1 = len(seq1) (or align length)
        len2 = 2 * MAXINDELS + 1  # k = len2 = 2d + 1
        if abs(len(seq1) - len(seq2)) > MAXINDELS:
            return math.inf, []

        E = []
        P = []
        for i in range(len1):  # t: O(kn)  s: O(kn)
            E.append([])
            P.append([])
            for j in range(len2):
                if len(seq2) >= i + j - MAXINDELS >= 0:
                    if i == 0:
                        E[i].append(5 * (j - MAXINDELS))
                        P[i].append(1)
                    else:
                        if j == len2 - 1:
                            e = min(self.diff(seq1[i - 1], seq2[i + j - MAXINDELS - 1]) + E[i - 1][j], INDEL + E[i][j - 1])
                            if e == INDEL + E[i][j - 1]:
                                P[i].append(1)
                            else:
                                P[i].append(0)
                        elif j == 0:
                            e = min(self.diff(seq1[i - 1], seq2[i + j - MAXINDELS - 1]) + E[i - 1][j],
                                    INDEL + E[i - 1][j + 1])
                            if e == INDEL + E[i - 1][j + 1]:
                                P[i].append(2)
                            else:
                                P[i].append(0)
                        else:
                            e = min(self.diff(seq1[i - 1], seq2[i + j - MAXINDELS - 1]) + E[i - 1][j], INDEL + E[i][j - 1],
                                    INDEL + E[i - 1][j + 1])
                            if e == INDEL + E[i][j - 1]:
                                P[i].append(1)
                            elif e == INDEL + E[i - 1][j + 1]:
                                P[i].append(2)
                            else:
                                P[i].append(0)
                        E[i].append(e)
                else:
                    E[i].append(math.inf)
                    P[i].append(math.inf)
        ret_j = len2 - 1
        ret = E[len1 - 1][ret_j]
        while ret == math.inf:
            ret_j -= 1
            ret = E[len1 - 1][ret_j]
        return ret, P

    #  Edit Distance Algorithm (Unrestricted Implementation)
    #  takes two sequences as arguments
    #  returns unrestricted edit distance and 2D array of back-pointers
    #  t: O(nm)  outer loop of double for-loop iterates n times. inner for-loop iterates m times.
    #  s: O(nm)  two 2-dimensional arrays (E & P) of size n x m are created.
    def getEditDistanceUnrestricted(self, seq1, seq2):
        len1 = min(len(seq1), self.MaxCharactersToAlign) + 1  # n = len1 = len(seq1) (or align length)
        len2 = min(len(seq2), self.MaxCharactersToAlign) + 1  # m = len2 = len(seq2) (or align length)
        E = []
        P = []
        for i in range(len1):  # t: O(n)  s: O(2n)
            E.append([5 * i])
            P.append([0])
        for j in range(1, len2):  # t: O(m)  s: O(2m)
            E[0].append(5 * j)
            P[0].append(1)
        for i in range(1, len1):  # t: O(nm)  s: O(2nm)
            for j in range(1, len2):
                e = min(self.diff(seq1[i - 1], seq2[j - 1]) + E[i - 1][j - 1], INDEL + E[i][j - 1], INDEL + E[i - 1][j])
                E[i].append(e)
                if e == INDEL + E[i][j - 1]:
                    P[i].append(1)
                elif e == INDEL + E[i - 1][j]:
                    P[i].append(2)
                else:
                    P[i].append(0)
        return E[len1 - 1][len2 - 1], P

    #  takes two characters as arguments
    #  returns MATCH if characters match; otherwise, returns SUB
    def diff(self, i, j):  # t: O(1)
        if i == j:
            return MATCH
        else:
            return SUB

    #  takes two sequences and a 2D array of back-pointers as arguments
    #  returns first 100 characters of each sequence aligned using the back-pointers
    def aligned(self, seq1, seq2, P):
        if not P:
            return "No Alignment Possible", "No Alignment Possible"
        if self.banded:
            return self.alignedBanded(seq1, seq2, P)
        else:
            return self.alignedUnrestricted(seq1, seq2, P)

    #  Alignment Extraction Algorithm (Banded Implementation)
    #  takes two sequences and a 2D array of back-pointers as arguments
    #  returns first 100 characters of each sequence aligned using the back-pointers
    def alignedBanded(self, seq1, seq2, P):
        i = len(P) - 1
        j = len(P[i]) - 1
        while P[i][j] == math.inf:
            j -= 1

        seq1_aligned = seq2_aligned = ''

        while not (i == 0 and j == MAXINDELS):
            if P[i][j] == 2:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = "-" + seq2_aligned
                i -= 1
                j += 1
            elif P[i][j] == 0:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = seq2[i + j - MAXINDELS - 1] + seq2_aligned
                i -= 1
            else:
                seq1_aligned = "-" + seq1_aligned
                seq2_aligned = seq2[i + j - MAXINDELS - 1] + seq2_aligned
                j -= 1

        return seq1_aligned[:100], seq2_aligned[:100]

    #  Alignment Extraction Algorithm (Unrestricted Implementation)
    #  takes two sequences and a 2D array of back-pointers as arguments
    #  returns first 100 characters of each sequence aligned using the back-pointers
    def alignedUnrestricted(self, seq1, seq2, P):
        i = len(P) - 1
        j = len(P[i]) - 1

        seq1_aligned = seq2_aligned = ''

        while i > 0 or j > 0:
            if P[i][j] == 2:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = "-" + seq2_aligned
                i -= 1
            elif P[i][j] == 0:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned
                i -= 1
                j -= 1
            else:
                seq1_aligned = "-" + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned
                j -= 1

        return seq1_aligned[:100], seq2_aligned[:100]
