import pandas as pd
import sys

#Storing the input file name and the output file name entered as parameters to the script
ip_file = sys.argv[1]
op_file = sys.argv[2]

#Dictionary that maintains the count of all the codons encountered in all the lines containing sequence
overall_count_map = {}

#Function that returns a dictionary containing the count of all the codons in the current sequence
def codon_count(str) :
    n = 3
    # Dictionary that maintains the count of all the codons encountered in the current sequence
    count_map = {}
    for i in range(0, len(str)-n-1, n):
        # ignore the extra characters when the size of the sequence is not a multiple of 3
        if len(str[i:i + n]) < 3:
            continue

        # current codon
        sub = str[i:i+n]
        # counter for the frequency of occurence of the current codon in the sequence
        count = 0
        if not sub in count_map :
            # iterate through the entire sequence
            for j in range(0, len(str)-n, n) :
                if str[j:j+n] == sub :
                    count+=1
            # update the dictionary with the counter for the current codon
            count_map[sub] = count
    return count_map

#Function that converts a dictionary to a csv file
def dict_to_csv(dict, filename) :
    df = pd.DataFrame.from_dict(dict, orient='index')
    df.to_csv(filename, header=False)


with open(ip_file, 'r') as f :
    # iterating through all the lines in the input file
    for curline in f:
        # ignoring the current line if it
        # starts with ">"
        if curline.startswith(">") or curline.startswith("\n"):
            continue
        else:
            # map containing the count of all the codons in the current sequence
            count_map = codon_count(curline)
            # updating the dictionary containing the overall count of all the codons with the values in the dictionary for the current line
            overall_count_map = {key: overall_count_map.get(key, 0) + count_map.get(key, 0)
                                                    for key in set(overall_count_map) | set(count_map)}

#Converting the dictionary containg the overall count for all codons to the output csv file
dict_to_csv(overall_count_map, op_file)
