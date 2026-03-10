
import numpy as np


def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    '''Calculate score for position (i,j) in scoring matrix, also record move to trace back

    Args:
        matrix (numpy array): scoring matrix
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        i (int): current row number
        j (int): current column number

    Returns:
        score in position (i,j)
        move to trace back: 0-END, 1-DIAG, 2-UP, 3-LEFT

    Pseudocode:
        Calculate scores based on upper-left, up, and left neighbors:
            diag_score = upper-left + (match or mismatch)
            up_score = up + gap
            left_score = left + gap
        score = max(0, diag_score, up_score, left_score)
        traceback = maximum direction or end

    '''
    end = 0
    diag = 1
    up = 2
    left = 3



    if i == 0:
        score = 0
        traceback = end
        return score, traceback
    elif j == 0:
        score = 0
        traceback = end
        return score, traceback
    else:
        if seq1[i - 1] == seq2[j - 1]:
            diag_score = matrix[i - 1, j - 1] + match
        else:
            if matrix[i - 1][j - 1] == 0:
                diag_score = 0
            else:
                diag_score = matrix[i - 1][j - 1] + mismatch


        if matrix[i - 1][j] == 0:
            up_score = 0
        else:
            up_score = matrix[i - 1][j] + gap

        if matrix[i][j - 1] == 0:
            left_score = 0
        else:
            left_score = matrix[i][j - 1] + gap

    scores = {"diag_score": diag_score, "up_score": up_score, "left_score": left_score}

    score = max(scores, key=scores.get)

    max_score = scores[score]

    if max_score == 0:
        traceback = end
        return max_score, traceback
    elif score == "diag_score" and max_score != 0:
        traceback = diag
        return max_score, traceback
    elif score == "up_score" and max_score != 0 :
        traceback = up
        return max_score, traceback
    elif score == "left_score" and max_score != 0:
        traceback = left
        return max_score, traceback


def max_scores(matrix):
    ''' Finds max score locations in matrix

    Args:
        matrix (numpy array): scoring matrix

    Returns:
        max_locs (list of tuples): list of tuples containing locations of max scores in scoring matrix
    '''

    score_dict = {}
    for i, row in enumerate(matrix):
        for j, col in enumerate(row):
            score_dict[(i, j)] = int(matrix[i, j])


    max_score = max(score_dict.values())
    scores = []
    for key in score_dict.keys():
        if score_dict[key] == max_score:
            scores.append(key)

    return scores


def traceback(seq1, seq2, traceback_matrix, maximum_position):
    '''Find the optimal path through scoring marix

        diagonal: match/mismatch
        up: gap in seq1
        left: gap in seq2

    Args:
        seq1 (str) : First sequence being aligned
        seq2 (str) : Second sequence being aligned
        traceback_matrix (numpy array): traceback matrix
        maximum_position (tuple): starting position to trace back from - highest score node

    Returns:
        aligned_seq1 (str): e.g. GTTGAC
        aligned_seq2 (str): e.g. GTT-AC

    Pseudocode:
        current_row, current_col = maximum_position
        while current_move != END:
            current_move = traceback_matrix[current_row][current_col]
            if current_move == DIAG:
                aligned_seq1 = seq1[current_col] + aligned_seq1
                aligned_seq2 = seq1[current_row] + aligned_seq2
                current_row -= 1
                current_col -= 1
            elif current_move == UP:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[current_row] + aligned_seq2
                current_col -= 1
            elif current_move == LEFT:
                aligned_seq1 = seq1[current_row] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                current_col -= 1

    '''
    END = 0
    DIAG = 1
    UP = 2
    LEFT = 3

    print(traceback_matrix)
    row = maximum_position[0]
    column = maximum_position[1]
    aligned_seq1 = ''
    aligned_seq2 = ''
    current_move = traceback_matrix[row][column]
    while current_move != END:
        if current_move == DIAG:
            p_seq1 = seq1[row - 1]
            aligned_seq1 += p_seq1
            p_seq2 = seq2[column - 1]
            aligned_seq2 += p_seq2

            row = row - 1
            column = column - 1
            current_move = traceback_matrix[row][column]

        elif current_move == UP:
            p_seq1 = seq1[row - 1]
            aligned_seq1 += p_seq1
            p_seq2 = '-'
            aligned_seq2 += p_seq2

            row = row - 1
            current_move = traceback_matrix[row][column]
        elif current_move == LEFT:
            p_seq1 = '-'
            aligned_seq1 += p_seq1
            p_seq2 = seq2[column - 1]
            aligned_seq2 += p_seq2

            column = column - 1
            current_move = traceback_matrix[row][column]

    return aligned_seq1[::-1], aligned_seq2[::-1]

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    '''Smith-Waterman algorithm for local alignment

    Args:
        seq1 (str): input seq 1
        seq2 (str): input seq 2
        match: default = +1
        mismatch: default = -1
        gap: default = -1

    Returns:
        alignments (list of tuples): [(aligned_seq1 (str), aligned_seq2 (str))]
        score_matrix (numpy array): scoring matrix
        max_score (int): max score in score_matrix

    Pseudocode:
        Define matrices based on seq lengths
        Call cal_score() to calculate scoring matrix and record traceback matrix
        all max_score() to determine max score locations in matrix
        Over each max location, to generate multiple best alignments:
        Call traceback() to generate local alignment
        Return values
        '''

    matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    traceback_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))


    for i, row in enumerate(matrix):
        for j, col in enumerate(row):
            max_score, traceback_val = cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap)
            matrix[i, j] = int(max_score)
            traceback_matrix[i, j] = int(traceback_val)

    print(matrix)
    scores = max_scores(matrix)

    aligned_seqs = []
    for score in scores:
        aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, score)
        aligned_seqs.append((aligned_seq1, aligned_seq2))


    return aligned_seqs


seq1 = 'TACTTAG'
seq2 = 'CACATTAA'
aligned_seqs = smith_waterman(seq1, seq2, match = 1, mismatch=-1, gap=-1)
print(aligned_seqs)
