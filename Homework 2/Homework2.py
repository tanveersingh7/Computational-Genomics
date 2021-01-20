####################################
# Created by : Tanveer Singh Virdi #
####################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import random


def needleman_wunsch(s1, s2) :
    nw_score_matrix, match_mismatch_matrix = fill_matrix(s1, s2)
    # Alignment score is the value in the bottom right corner of the Needleman-Wunsch scoring matrix
    alignment_score = nw_score_matrix[nw_score_matrix.shape[0]-1][nw_score_matrix.shape[1]-1]
    # Alignment sequences obtained
    al1, al2 = get_alignments(s1, s2, nw_score_matrix, match_mismatch_matrix)

    return alignment_score, al1, al2

def fill_matrix(s1, s2) :
    rows = len(s1) + 1
    cols = len(s2) + 1

    score_matrix = np.zeros((rows, cols))
    match_mismatch_matrix = np.zeros((len(s1), len(s2)))

    # Fill up the match_mismatch_matrix based on the matches and mismatches between the two sequences
    for i in range(len(s1)) :
        for j in range(len(s2)) :
            # Fill the current cell with MATCH_SCORE if the character in the ith position of sequence1 matches with the
            # jth position of sequence 2. Otherwise fill the current cell with MISMATCH_SCORE
            match_mismatch_matrix[i][j] = MATCH_SCORE if s1[i] == s2[j] else MISMATCH_SCORE

    # Filling the score_matrix

    # Filling first column
    for i in range(rows) :
        score_matrix[i][0] = i * GAP_PENALTY
    # Filling first row
    for j in range(cols) :
        score_matrix[0][j] = j * GAP_PENALTY

    for i in range(1, rows) :
        for j in range(1, cols) :
                # Filling the score in the current cell
                score_matrix[i][j] = max(score_matrix[i-1][j-1] + match_mismatch_matrix[i-1][j-1], score_matrix[i-1][j]
                                         + GAP_PENALTY, score_matrix[i][j-1] + GAP_PENALTY)

    return score_matrix, match_mismatch_matrix

def get_alignments(s1, s2, score_matrix, match_mismatch_matrix) :
    # Variables for storing alignments
    a1 = ""
    a2 = ""
    # We will start tracing back from the bottom right corner of the scores matrix
    i = len(s1)
    j = len(s2)

    # Iterate till we reach the top left corner of the scores matrix
    while i > 0 and j > 0 :
        # Update the alignments if the score in the current cell was obtained from the prev diagonal cell
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + match_mismatch_matrix[i - 1][j - 1] :
            a1 = s1[i-1] + a1
            a2 = s2[j-1] + a2
            i -= 1
            j -= 1
        # Update the alignments if the score in the current cell was obtained from the prev top cell
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + GAP_PENALTY :
            a1 = s1[i-1] + a1
            a2 = GAP_CHAR + a2
            i -= 1
        # Update the alignments if the score in the current cell was obtained from the prev left cell
        elif j > 0 and score_matrix[i][j] == score_matrix[i][j - 1] + GAP_PENALTY :
            a1 = GAP_CHAR + a1
            a2 = s2[j-1] + a2
            j -= 1

    # Continue tracing to the top left corner
    while i > 0:
        a1 = s1[i - 1] + a1
        a2 = GAP_CHAR + a2
        i -= 1
    while j > 0:
        a1 = GAP_CHAR + a1
        a2 = s2[j - 1] + a2
        j -= 1

    return a1, a2

def anchoring_needleman_wunsch(s1, s2, anchor_positions) :
    # Variables for storing the final alignments
    al1 = ""
    al2 = ""
    # Total Alignment score
    alignment_score = 0
    # Start index for iterating over sequence s1
    start_index1 = 0
    # Start index for iterating over sequence s2
    start_index2 = 0

    # Iterating through all the anchor positions specified in different rows of the Matches file
    for position in anchor_positions :
        # Start position of the anchor segment of sequence s1
        ap_start1 = position[0]
        # End position of the anchor segment of sequence s1
        ap_end1 = position[1]
        # Start position of the anchor segment of sequence s2
        ap_start2 = position[2]
        # End position of the anchor segment of sequence s2
        ap_end2 = position[3]

        # Perform Needleman Wunsch on the non-anchor segments of both sequences
        sub_score, sub_al1, sub_al2 = needleman_wunsch(s1[start_index1:ap_start1], s2[start_index2:ap_start2])
        # Adding the alignment score for the non-anchor segments and alignment score for the anchor segments of both
        # sequences to the total alignment score
        alignment_score += sub_score + anchor_seq_score(s1[ap_start1:ap_end1+1], s2[ap_start2:ap_end2+1])
        # Adding the alignments and anchor segments of the sequences to the final alignments
        al1 += sub_al1 + s1[ap_start1:ap_end1+1]
        al2 += sub_al2 + s2[ap_start2:ap_end2+1]

        # Updating the start positions to the start of next non-anchor segments of both sequences
        start_index1 = ap_end1 + 1
        start_index2 = ap_end2 + 1

    # Perform Needleman Wunsch on the last non-anchored segment of both sequences
    sub_score, sub_al1, sub_al2 = needleman_wunsch(s1[start_index1:], s2[start_index2:])
    # Adding alignment score of the final non-anchored segments to the total alignment score
    alignment_score += sub_score
    # Adding the alignments of the last non-anchor segments to the final alignments
    al1 += sub_al1
    al2 += sub_al2

    return alignment_score, al1, al2

def anchor_seq_score(s1, s2) :
    # Variable storing the total alignment score for the anchor segments of the two sequences
    score = 0
    for i in range(len(s1)) :
        # MATCH_SCORE added to score if the corresponding positions in the anchor segments match
        if s1[i] == s2[i] :
            score += MATCH_SCORE
        else :
            # MISMATCH_SCORE added to score if the corresponding positions in the anchor segments don't match
            score += MISMATCH_SCORE
    return score

def get_sequence(file_name) :
    with open(file_name, 'r') as f:
        # iterating through all the lines in the input file
        for curline in f:
            # ignoring the current line if it starts with ">"(fasta header)
            if curline.startswith(">") or curline.startswith("\n"):
                continue
            else:
                #returning the sequence in the file
                return curline.strip('\n').strip()

def get_anchor_positions(file_name) :
    # Array of all anchor positions in the matches file
    anchor_positions = []
    with open(file_name, 'r') as f :
        # Skipping the first line of the file
        next(f)
        # Iterating through all lines of the file having anchor positions for both sequences
        for curline in f :
            positions = curline.strip('\n').strip().split()
            positions = [int(i) - 1 for i in positions]
            # Storing the anchor positions in the array
            anchor_positions.append(positions)

    return anchor_positions

def get_permutation(seq) :
    l = list(seq)
    random.shuffle(l)
    # Returning a random permutation of the sequence
    return ''.join(l)

def permutations(s1, s2) :
    #Array for storing all the 100 alignment scores
    alignment_scores = []
    for i in range(100) :
        #Permuting both sequences
        s1 = get_permutation(s1)
        s2 = get_permutation(s2)
        #Random Alignment score for this particular permutation
        score = needleman_wunsch(s1, s2)[0]
        #Adding the random alignment score to the array
        alignment_scores.append(score)
    return alignment_scores

def write_results(file_name, alignment_score, a1, a2, flag) :
    with open(file_name, 'w') as f :
        if flag :
            f.write("Results for Anchor version of Needleman Wunsch :\n ")
        else :
            f.write("Results for Standard version of Needleman Wunsch :\n ")
        f.write("\nAlignment score : "+ str(alignment_score))
        f.write("\nAlignment 1 :\n "+a1)
        f.write("\nAlignment 2 :\n "+a2)
        f.write("\n")


#Main function
if __name__ == '__main__':
    MATCH_SCORE = 1
    MISMATCH_SCORE = -3
    GAP_PENALTY = -2
    GAP_CHAR = "_"

    seq1_file_name = sys.argv[1]
    seq2_file_name = sys.argv[2]

    #Getting the sequences from the FASTA sequence files
    seq1 = get_sequence(seq1_file_name)
    seq2 = get_sequence(seq2_file_name)

    #Output file where the results will be written to
    output_file_name = "Output.txt"
    if (len(sys.argv) > 3):
        print("Running Anchor version of Needleman Wunsch: ")
        matches_file = sys.argv[3]
        # Getting the anchor positions from the Matches file
        anchor_positions = get_anchor_positions(matches_file)
        # Running the Anchor version of Needleman Wunsch on the sequences and getting the required results
        alignment_score, alignment_seq1, alignment_seq2 = anchoring_needleman_wunsch(seq1, seq2, anchor_positions)
        # Writing the results to the output file
        write_results(output_file_name, alignment_score, alignment_seq1, alignment_seq2, True)
        print("Results written to Output.txt")
    elif (len(sys.argv) == 3):
        print("Running Standard version of Needleman Wunsch : ")
        # Running the standard version of Needleman Wunsch on the sequences and getting the required results
        alignment_score, alignment_seq1, alignment_seq2 = needleman_wunsch(seq1, seq2)
        # Writing the results to the output file
        write_results(output_file_name, alignment_score, alignment_seq1, alignment_seq2, False)
        print("Results written to Output.txt")

        print("Permuting the nucleotides in the sequences and getting the random alignment scores for plotting the "
              "histogram. ")
        alignment_scores = permutations(seq1, seq2)
        # Plotting the distribution of random alignment scores
        plt.title("Histogram showing the distribution of random alignment scores")
        plt.xlabel("Random Alignment Scores")
        plt.hist(alignment_scores, color= 'r')
        plt.axvline(x=alignment_score, color="black", label='Optimal Score')
        plt.show()
