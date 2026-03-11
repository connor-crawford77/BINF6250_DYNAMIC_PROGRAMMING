# Introduction
The aim of this project was to create an algorithm to implement pairwise alignment. This program uses the Smith-Waterman method for local alignment. The program takes in two sequences and finds the best local alignments between the two sequences. The process of sequence alignment becomes too complex if you try to naively make every comparison between two sequences. Dynamic programming implementation breaks alignment in to smaller sub-problems and reduces the time-complexity of alignment.

# Pseudocode

```
Make two matrices, read in 2 seqs, scores 


Starting at i=1, j=1, iterate over matrix to len(seq1), len(seq2): 

    cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap): 
        Calculate scores based on upper-left, up, and left neighbors: 
            diag_score = upper-left + (match or mismatch) 
            up_score = up + gap 
            left_score = left + gap 
        score = max(0, diag_score, up_score, left_score) 
        traceback = maximum direction or end (0-END, 1-DIAG, 2-UP, 3-LEFT) 
        return score, move 
    
    Update score matrix and traceback matrix at i, j 
 

Find maxmimum scores and make a list of them:

Def max_scores(): read in score matrix and return max score as well as list of tuples of max score locations

for list of max scores: 
Aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, maximum_position) 
    #Note: max position = (current_row, current_column), first time is max score 

Alignseq1, 2 = [] 

current_row, current_col = maximum_position 
    while current_move != END:   
        current_move = traceback_matrix[current_row][current_col]
        if current_move == DIAG: 
            aligned_seq1 = seq1[current_col] + aligned_seq1  
            aligned_seq2 = seq1[current_row] + aligned_seq2 
            current_row = current_row -1 
            current_col = current_col -1 

        elif current_move == UP: 
            aligned_seq1 = "-” + aligned_seq1  
            aligned_seq2 = seq1[current_row] + aligned_seq2 
            current_col -= 1 

        elif current_move == LEFT: 
            as above except opposite 
            ... 

 	  Returns: 
     	   aligned_seq1 (str): e.g. GTTGAC 
        	aligned_seq2 (str): e.g. GTT-AC 


Get each alignment for each maximum, separately. Could do a list of tuples 

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1): 
    Max Score ()
    Define the matrices etc 
    Call score 
    Call find max 
    Call traceback for all maxes 
    Return list of tuples and scoring matrix 

 

Outside functions:  
Print all alignments 
Print scoring matrix 
Print max score if we want because we already know it we just need to pass it back from smithwaterman 
```

# Successes
Our group had great success from a productive psuedocode session. Due to this, we all felt comfortable with the implementation step. Being able to visualize the algorithm flow and having only two sequences at a time helped wrapping our brains around things.

# Struggles
We did have a few stumbles around the indexing, where we needed extra examples to test various correct alignment possibilities to debug. Making sure we had enough positions in the matrices was crucial. Initially we were constructing a matrix with the length of the rows and columns being equal to the length of the sequences. We had to adjust this and add one extra row and column to represent our initial diagonal state. We also had to be careful when accessing neighboring cells (diagonal, up, and left) to ensure we were referencing the correct positions in the matrix.

# Personal Reflections
## Group Leader
## Aaronie Jersha Jenyfred
A lot of our progress came from the pseudocode discussions we had together as a group. Having previously studied some of these concepts in genomics helped me frame the logic behind the core parts of the alignment algorithm. Compared to some of the earlier projects, this one felt less intense, which gave me more time to experiment with the code and understand how each function interacted with the others.

## Group members
## Connor Crawford
The scale of the data definitely made the project smoother this week. Especially when it came to debugging, when I ran into stumbling blocks it was much easier to figure out what was going wrong in my code. Spending a couple hours with my other group members making the pseudocode made implementation very smooth.

## Victoria Van Berlo
This week's project was made a lot easier by having great pseudocode. I feel like I'm getting better at the planning phase. This is one of the most tangibly-understandable projects for me, so I felt pretty confident about our implementation.

# Generative AI Appendix
Our group did not use generative AI for this project. 
