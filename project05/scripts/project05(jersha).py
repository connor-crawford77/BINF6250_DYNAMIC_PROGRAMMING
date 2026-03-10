import numpy as np
def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    END = 0
    DIAG = 1
    UP = 2
    LEFT = 3

    # to compare the current character from both the sequences
    if seq1[i - 1] == seq2[j - 1]:
        diag_score = matrix[i - 1][j - 1] + match
    else:
        diag_score = matrix[i - 1][j - 1] + mismatch

    # scoring for top to down
    up_score = matrix[i - 1][j] + gap

    # scoring for left to right
    left_score = matrix[i][j - 1] + gap

    # to find the highest scoring value
    score = max(0, diag_score, up_score, left_score)

    # if 0 is the best score, stop the traceback
    if score == 0:
        move = END
    elif score == diag_score:
        move = DIAG
    elif score == up_score:
        move = UP
    else:
        move = LEFT

    return score, move

def max_scores(matrix):

    # find the highest value in the matrix
    max_score = np.max(matrix)

    # to find all the positions where that particular value occurs
    # np.where to return 2 arrays - array of rows and array of columns
    position = np.where(matrix == max_score)

    max_position = list(zip(position[0], position[1]))

    return max_position


def traceback(seq1, seq2, traceback_matrix, max_position):
    END = 0
    DIAGONAL = 1
    UP = 2
    LEFT = 3

    aligned_seq1 = ""
    aligned_seq2 = ""

    current_row, current_col = max_position

    current_move = traceback_matrix[current_row][current_col]

    # tracebacking until END is reached
    while current_move != END:

        if current_move == DIAGONAL:
            aligned_seq1 = seq1[current_row - 1] + aligned_seq1
            aligned_seq2 = seq2[current_col - 1] + aligned_seq2
            current_row -= 1
            current_col -= 1

        elif current_move == UP:
            aligned_seq1 = seq1[current_row - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            current_row -= 1

        elif current_move == LEFT:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[current_col - 1] + aligned_seq2
            current_col -= 1

        current_move = traceback_matrix[current_row][current_col]

    return aligned_seq1, aligned_seq2


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    num_rows = len(seq1) + 1
    num_cols = len(seq2) + 1

    score_matrix = np.zeros((num_rows, num_cols), dtype=int)
    traceback_matrix = np.zeros((num_rows, num_cols), dtype=int)

    for i in range(1, num_rows):
        for j in range(1, num_cols):
            score, move = cal_score(score_matrix, seq1, seq2,
                                    i, j, match, mismatch, gap)
            score_matrix[i][j] = score
            traceback_matrix[i][j] = move

    max_position = max_scores(score_matrix)
    max_score = score_matrix[max_position[0][0]][max_position[0][1]]

    alignments = []
    for position in max_position:
        aligned_seq1, aligned_seq2 = traceback(seq1, seq2,
                                               traceback_matrix, position)
        alignments.append((aligned_seq1, aligned_seq2))

    return alignments, score_matrix

seq1 = 'TACTTAG'
seq2 = 'CACATTAA'

alignments, score_matrix = smith_waterman(seq1, seq2)

for aligned_seq1, aligned_seq2 in alignments:
    print(aligned_seq1)
    print(aligned_seq2)

print(score_matrix)
