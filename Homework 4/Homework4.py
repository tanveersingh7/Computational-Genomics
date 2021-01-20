##############################################
#  Created by : Tanveer Singh Virdi          #
#  Discussion with : Mohana Krishna Vutukuru #
##############################################

import sys
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks, peak_widths

'''
Function to get all Sequences from the input FASTA file
'''
def get_sequences(file_name) :
    #list which stores all sequences in the file
    sequences = []
    with open(file_name, 'r') as f:
        # iterating through all the lines in the input file
        for curline in f:
            curline = curline.rstrip()
            # ignoring the current line if it starts with ">"(fasta header)
            if curline.startswith(">") or curline.startswith("\n"):
                continue
            else:
                #adding the sequence in the current line to the list
                sequences.append(curline)
    return sequences

'''
Function to get "counts_map" : a list of dictionaries for each position maintaining the count of each conserved(A,C,G,T)
and non conserved character in each position of the sequence
'''
def get_count_map(seqs) :
    #list of dictionaries maintaining the count of each conserved(A,C,G,T) and non conserved(Non) characters
    #in each position of the sequence
    counts_map = [{} for sub in range(len(seqs[0]))]

    #Iterating through all the sequences
    for i in range(len(seqs)):
        #Iterating through all the positions in the sequence
        for j in range(len(seqs[i])):
            #Incrementing the count of the character in the current position of the current sequence
            if seqs[i][j] == 'A':
                if not 'A' in counts_map[j]:
                    counts_map[j]['A'] = 1
                else:
                    counts_map[j]['A'] += 1
            elif seqs[i][j] == 'C':
                if not 'C' in counts_map[j]:
                    counts_map[j]['C'] = 1
                else:
                    counts_map[j]['C'] += 1
            elif seqs[i][j] == 'G':
                if not 'G' in counts_map[j]:
                    counts_map[j]['G'] = 1
                else:
                    counts_map[j]['G'] += 1
            elif seqs[i][j] == 'T':
                if not 'T' in counts_map[j]:
                    counts_map[j]['T'] = 1
                else:
                    counts_map[j]['T'] += 1
            else:
                if not 'Non' in counts_map[j]:
                    counts_map[j]['Non'] = 1
                else:
                    counts_map[j]['Non'] += 1

    for count_map in counts_map :
        #Getting the total count of all characters in the current position
        total = sum(count_map.values())

        for k, v in count_map.items():
            #Percentage presence of each character in the current position
            count_map[k] = v / total

    return counts_map

'''
Returns a list containing tuples having the conserved character and the conservation rate for each position
'''
def get_conserved_rates(counts_map) :
    conserved = []
    conserved_sorted = []
    #Sort the count map of each position in descending order
    for count_map in counts_map:
        conserved_sorted.append(sorted(count_map, key=lambda x: count_map[x], reverse=True))
    #For each position in the sequence get the most common base and its conservation rate and add it to the tuple
    for i in range(len(conserved_sorted)) :
        if conserved_sorted[i][0] == 'Non' :
            conserved.append((conserved_sorted[i][1], counts_map[i][conserved_sorted[i][1]]))
        else :
            conserved.append((conserved_sorted[i][0], counts_map[i][conserved_sorted[i][0]]))

    return conserved

'''
Using the 'savgol_filter' function in scipy for smoothing
'''
def smooth_function(conservation_rates) :
    conservation_rates_smoothed = savgol_filter(conservation_rates, window_length=95, polyorder=2)
    return conservation_rates_smoothed

'''
Plotting the variability against the position in the gapped alignment. Smoothing function has been used for the 
conservation rates of each position 
 '''
def variability_plot(conservation_rates_smoothed) :
    plt.plot(conservation_rates_smoothed)
    plt.xlabel("Position in Gapped Alignment")
    plt.ylabel("Conservation Rate %")
    plt.title("Question 2 : Variability plot")
    plt.show()

'''
Plotting the variability against the position in the gapped alignment along with the variable regions
'''
def variability_plot_with_variable_regions(conservation_rates_smoothed, var_region_left, var_region_right) :
    var_region_y = [50] * len(var_region_left)
    plt.plot(conservation_rates_smoothed)
    plt.hlines(var_region_y, var_region_left, var_region_right, linestyles='solid')
    plt.xlabel("Position in Gapped Alignment")
    plt.ylabel("Conservation Rate %")
    plt.title("Question 4 : Variability plot with variable regions")
    plt.show()

'''
Returns the lists of left and right coordinates of the variable regions
'''
def find_variable_regions(conservation_rates_smoothed) :
    # Coordinates of the peaks('find_peaks' function from the scipy package is used to identify the local
    # minimas(variable regions)
    peak_regions, _ = find_peaks([(-1 * x) for x in conservation_rates_smoothed], distance = 130)
    # Left and Right coordinates of the variable regions
    _ , _ , left, right = peak_widths([(-1 * x) for x in conservation_rates_smoothed], peak_regions, rel_height=0.7)

    #Storing the left and right coordinates of the variable regions
    var_region_left = []
    var_region_right = []

    for i in range(len(left)):
        var_region_left.append(int(left[i]))
        var_region_right.append(int(right[i]))

    return var_region_left, var_region_right


#Main function
if __name__ == '__main__':
    conserved_file_name = 'solution-problem-1.txt'
    variable_region_coordinates_file = 'solution-problem-3.txt'
    seq_file_name = sys.argv[1]
    #Getting the sequences from the FASTA sequence files
    seqs = get_sequences(seq_file_name)
    # getting list of dictionaries maintaining the count of each conserved(A,C,G,T) and non conserved(Non) characters
    # in each position of the sequence
    counts_map = get_count_map(seqs)
    # getting a list containing tuples having the conserved character and the conservation rate for each position
    conserved = get_conserved_rates(counts_map)

    # Question 1 : Write the conservation rate and the conserved character of each position to the output file
    with open(conserved_file_name, 'w') as file :
        for value in conserved:
            file.write(str(value[0]) + '\t' + str(value[1]*100) + '\n')

    conservation_rates = [x[1] * 100 for x in conserved]
    #Performing smoothing on the conservation rates
    conservation_rates_smoothed = smooth_function(conservation_rates)

    #Question 2 : Plot the variability against the position in the gapped alignment.
    variability_plot(conservation_rates_smoothed)

    #Question 3 : Find the left and right coordinates of the variable regions
    var_region_left, var_region_right = find_variable_regions(conservation_rates_smoothed)
    #Write the left and right coordinates of each variable region to the output file
    with open(variable_region_coordinates_file, 'w') as file :
        for i in range(len(var_region_left)):
            file.write(str(var_region_left[i]) + '\t' + str(var_region_right[i]) + '\n')

    #Question 4 : Plot the variability against the position in the gapped alignment along with the variable regions.
    variability_plot_with_variable_regions(conservation_rates_smoothed, var_region_left, var_region_right)
