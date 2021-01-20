####################################
# Created by : Tanveer Singh Virdi #
####################################

import sys
import numpy as np
import pandas as pd
import random

#Class defining a tree structure storing the nodes and their children and their corresponding distances
class TreeNode :

    def __init__(self, id):
        self.id = id
        self.children = None

    def add_child(self, child_id, distance):
        if not self.children :
            self.children = {}
        self.children[child_id] = distance


#Reads the fasta file and stores the sequences and their labels
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

#Calaculates the genetic distance score between two sequences
def genetic_distance_score(s1, s2) :
    # Variable storing the genetic distance score for the two sequences
    score = 0
    for i in range(len(s1)) :
        if s1[i] != s2[i] :
            # MISMATCH_SCORE added to score if the corresponding positions in the segments don't match
            score += MISMATCH_SCORE
    return score/len(s1)

#Generates the genetic distance matrix for all the input sequences
def get_dist_matrix(seqs) :
    gen_dist_matrix = []
    for i in range(len(seqs)):
        matrix = []
        for j in range(len(seqs)):
            matrix.append(genetic_distance_score(seqs[i], seqs[j]))
        gen_dist_matrix.append(matrix)
    return gen_dist_matrix

# Calculates the Q matrix from the distance matrix and returns the two closest nodes and their distance
def get_qmatrix(gen_dist_matrix):
    N = gen_dist_matrix.shape[0]
    #Initializing the Q matrix
    q_matrix = [[0] * N for _ in range(N)]
    #Variables for storing the indexes of the two closest nodes and the min distance
    min_idx1 = 0
    min_idx2 = 0
    min_dist = float('inf')
    for i in range(N) :
        for j in range(N) :
            if i != j :
                q_matrix[i][j] = (N-2) * gen_dist_matrix.iloc[i,j] - gen_dist_matrix.iloc[i].sum() - gen_dist_matrix.iloc[j].sum()
                #Updating the min distance and the indexes of the two closest nodes
                if q_matrix[i][j] < min_dist :
                    min_dist = q_matrix[i][j]
                    min_idx1,min_idx2 = i,j
    min_i, min_j = gen_dist_matrix.index.values[min_idx1], gen_dist_matrix.index.values[min_idx2]
    return q_matrix, min_dist, min_i, min_j

#Calculates the distance of each of the nodes being joined to the new node being created
def get_dist_to_new_node(gen_dist_matrix, min_idx_i, min_idx_j) :
    N = gen_dist_matrix.shape[0]
    dist1 = 1 / 2 * gen_dist_matrix.loc[min_idx_i, min_idx_j] + 1 / (2 * (N - 2)) * (gen_dist_matrix.loc[min_idx_i].sum() - gen_dist_matrix.loc[min_idx_j].sum())
    dist2 = gen_dist_matrix.loc[min_idx_i, min_idx_j] - dist1
    return dist1, dist2

#Updating the genetic distance matrix at each step
def get_next_gd_matrix(df_gd_matrix, min_idx_i, min_idx_j, seq_count) :
    #Mainting a copy of the current genetic distance pandas dataframe
    df_gd_matrix_copy = df_gd_matrix.copy()
    #Dropping the rows and columns from the genetic distance matrix corresponding to the nodes having the min distance
    df_gd_matrix = df_gd_matrix.drop(index=[min_idx_i, min_idx_j], columns=[min_idx_i, min_idx_j])
    new_size = df_gd_matrix.shape[0]
    #Adding a row of zeros and a column of zeros corresponding to the new node being created
    zero_row = pd.DataFrame(np.zeros((1, new_size)), index=[seq_count], columns=df_gd_matrix.index.values)
    zero_col = pd.DataFrame(np.zeros((new_size, 1)), index=df_gd_matrix.index.values, columns=[seq_count])
    #Adding this zero row and column to the distance matrix
    df_gd_matrix = zero_row.append(df_gd_matrix)
    df_gd_matrix = pd.concat([zero_col, df_gd_matrix], axis=1)

    indexes = list(df_gd_matrix.index.values)
    indexes.remove(seq_count)
    #Reindexing the pandas dataframe so that the first row and first column correspond to the new node
    df_gd_matrix = df_gd_matrix.reindex(index = [seq_count] + indexes, columns = [seq_count] + indexes)
    #For every other node we calculate their distance to this new node
    for index in indexes:
        dist = 1 / 2 * (df_gd_matrix_copy.loc[min_idx_i, index] + df_gd_matrix_copy.loc[min_idx_j, index] -
                        df_gd_matrix_copy.loc[min_idx_i, min_idx_j])
        df_gd_matrix.loc[index, seq_count] = dist
        df_gd_matrix.loc[seq_count, index] = dist

    df_gd_matrix.loc[seq_count, seq_count] = 0

    return df_gd_matrix

#Implementing the Nei-Saitou neighbor joining algorithm
def neighbor_joining(df_gd_matrix) :
    N = df_gd_matrix.shape[0]
    #Numbering of the internal nodes starting from 120
    seq_count = 120
    #Edges map storing the mapping of the child nodes to their parent nodes and their distances
    edges_map = {}
    #Variable storing the final root id
    root_id = 0
    #Iterate until the size of the genetic distance matrix shrinks to 3
    while N > 2 :
        if N == 3 :
            #The root node will have three children. Edges map will store the mapping of the new node created and its
            # distance to the three children
            dist1 = 1/2 * df_gd_matrix.iloc[0,1] + 1/(2 * (N - 2)) * (df_gd_matrix.iloc[0].sum() - df_gd_matrix.iloc[1].sum())
            dist2 = df_gd_matrix.iloc[0, 1] - dist1
            dist3 = df_gd_matrix.iloc[0,2] - dist1
            edges_map[df_gd_matrix.index.values[0]] = (seq_count,dist1)
            edges_map[df_gd_matrix.index.values[1]] = (seq_count,dist2)
            edges_map[df_gd_matrix.index.values[2]] = (seq_count, dist3)
            #Root id is the sequence number of the final node added
            root_id = seq_count
            break
        else :
            #Variables storing the q matrix and the two closest nodes
            q_matrix, min_dist, min_idx1, min_idx2 = get_qmatrix(df_gd_matrix)
            #Variables storing the distance of the two nodes to the new node created
            dist1, dist2 = get_dist_to_new_node(df_gd_matrix, min_idx1, min_idx2)
            #Mapping the new node created and its distance to the two nodes in the edges map
            edges_map[min_idx1] = (seq_count, dist1)
            edges_map[min_idx2] = (seq_count, dist2)

            #Updating the genetic distance matrix
            df_gd_matrix = get_next_gd_matrix(df_gd_matrix, min_idx1, min_idx2, seq_count)
            #Id of the internal node decremented by 1
            seq_count -= 1
            N = df_gd_matrix.shape[0]

    return edges_map, root_id

#Constructing a tree using the edges map using Breadth First Search
def create_tree(edges_map, root_id) :
    root = TreeNode(root_id)
    #Adding the root node to the queue first
    queue = [root]
    #Looping till the queue becomes empty
    while len(queue) > 0:
        node_list = []
        for node in queue:
            for child, (parent, distance) in edges_map.items():
                #If the current node id matches the parent
                if parent == node.id:
                    #Create new node corresponding to the child
                    childNode = TreeNode(child)
                    #Add the new child node to the parent node
                    node.add_child(childNode, distance)
                    #add the child node to the list
                    node_list.append(childNode)
        for node in node_list:
            #Remove the mappings from the edges map whose nodes have been created
            del edges_map[node.id]
        queue = node_list
    return root

#Preorder traversal of the tree
def pre_order(root, pre_order_edges_list) :
    if root is None or root.children is None :
        return
    for child, distance in root.children.items() :
        #Add the node to the list first
        pre_order_edges_list.append((root.id, child.id, distance))
        #Traverse all the children recursively
        pre_order(child, pre_order_edges_list)
    return pre_order_edges_list

#Potorder traversal of the tree
def postorder_traversal(root):
    if root is None :
        return
    if root.children is None:
        return id_label_map[root.id]
    visited_nodes = []
    # Traverse the children first recursively
    for child, distance in root.children.items():
        visited_nodes.append(postorder_traversal(child) + ':' + str(distance))
    newick_string = '(' + ','.join(visited_nodes) + ')'
    return newick_string

#Generating bootstrap sequences from the input sequences
def bootstrap(seqs) :
    #Variable storing the bootstrap sequences
    bootstrap_seqs = []
    #Variable storing the length of sequence
    len_seq = len(seqs[0])
    #Indices picked for repeated sampling with replacement
    indices_to_pick = []
    for i in range(len_seq) :
        indices_to_pick.append(random.randint(0,len_seq-1))

    #generating the bootstrap sequences
    for i in range(len(seqs)) :
        bootstrap_seq = ''
        for index in indices_to_pick :
            bootstrap_seq += seqs[i][index]
        bootstrap_seqs.append(bootstrap_seq)
    return bootstrap_seqs

#Generating a list of partition nodes
def partitions_of_node(node, nodes_list) :
    if node is None or node.children is None :
        return
    if nodes_list is None :
        nodes_list = []
    nodes_list.append(node.id)
    for node_child in node.children.keys() :
        partitions_of_node(node_child, nodes_list)
    return nodes_list

#Generating a map of partitions from the root node
def map_of_partitions(root) :
    q = []
    q.append(root)
    partitions_map = {}
    #Breadth First Search used for generating the partitons map
    while q :
        node = q.pop(0)
        if node.children is None :
            continue
        partitions_list = partitions_of_node(node, None)
        partitions_map[node.id] = partitions_list
        for node_child in node.children.keys():
            q.append(node_child)
    return partitions_map

#Main function
if __name__ == '__main__':
    MISMATCH_SCORE = 1
    # Output files where the results will be written to
    edges_file_name = "edges.txt"
    tree_file_name = "tree.txt"
    bootstrap_file_name = "bootstrap.txt"
    #Input file
    seq_file_name = sys.argv[1]
    # Lists for storing the sequences and their labels
    labels = []
    seqs = []

    #Getting the sequences from the FASTA sequence files
    with open(seq_file_name) as fp:
        for name, seq in read_fasta(fp):
            labels.append(name[1:])
            seqs.append(seq)

    #Sequence labels are mapped to tip ids 1,2,...,61
    ids = [i for i in range(1, len(labels)+1)]
    #Map containing mapping of tip ids to sequence labels
    id_label_map = dict(zip(ids, labels))
    #Generating the genetic distance matrix
    gen_dist_matrix = get_dist_matrix(seqs)
    #Converting the genetic distance matrix to a pandas dataframe and writing this to a tab-delimited file with rows and
    # columns labeled by the sequence identifiers.
    gen_dist_matrix_df = pd.DataFrame(data = np.array(gen_dist_matrix), index = labels, columns = labels)
    gen_dist_matrix_df.to_csv("genetic_distance.txt", sep='\t', mode='w')

    # Converting the genetic distance matrix to a pandas dataframe
    df_gd_matrix = pd.DataFrame(data=np.array(gen_dist_matrix),index=ids, columns=ids)

    #Running the Nei-Saitou neighbor joining algorithm and getting the edges map having the mappings from the
    # child nodes to their parent nodes and their corresponding distances
    edges_map, root_id = neighbor_joining(df_gd_matrix)
    #Creating a tree from the edges_map having the parent child mapping and their distances and getting the root node
    # from this tree
    root_node = create_tree(edges_map, root_id)

    #Preorder traversal for the edges matrix
    pre_order_edges_list = []
    pre_order_edges_list = pre_order(root_node, pre_order_edges_list)
    #Writing the preorder traversal to the edges file
    with open(edges_file_name, 'w') as file :
        for edge in pre_order_edges_list:
            file.write(str(edge[0]) + '\t' + str(edge[1]) + '\t' + str(edge[2]) + '\n')
    #Variables storing the three children of the root node and their distances to the root node
    (child1, dist1), (child2, dist2), (child3, dist3) = root_node.children.items()
    #Variable storing the postorder traversal in Newick format
    post_order_newick_output = '((' + postorder_traversal(child2) + ':' + str(dist2) + ',' + postorder_traversal(child3) + ':' + str(dist3) + '):' + str(dist1) + ');'
    # Writing the postorder traversal to the newick file
    with open(tree_file_name, 'w') as file:
        file.write(post_order_newick_output)

    #List of partitions of the root node and the partitions under them
    root_partitions = partitions_of_node(root_node, None)
    #Map storing the mapping of all the nodes to their partitions
    partitons_map = map_of_partitions(root_node)

    matching_partitions = [0] * 59 #The tree has 59 internal nodes
    #100 inferences
    for i in range(100) :
        print("Iteration : "+ str(i))
        #Generating the bootstrap sequences
        bootstrap_seqs = bootstrap(seqs)
        #Genetic distance matrix for the bootstrap sequences
        bootstrap_gen_dist_matrix = get_dist_matrix(bootstrap_seqs)
        # Converting the genetic distance matrix to a pandas dataframe
        df_bootstap_gd_matrix = pd.DataFrame(data=np.array(bootstrap_gen_dist_matrix),index=ids, columns=ids)

        # Running the Nei-Saitou neighbor joining algorithm and getting the edges map having the mappings from the
        # child nodes to their parent nodes and their corresponding distances
        btstr_edges_map, btstr_root_id = neighbor_joining(df_bootstap_gd_matrix)
        # Creating a tree from the edges_map having the parent child mapping and their distances and getting the root node
        # from this tree
        btstr_root_node = create_tree(btstr_edges_map, btstr_root_id)
        #Generating the partitions for the tree nodes corresponding to the bootstrap sequences
        btstr_partitions_map = map_of_partitions(btstr_root_node)
        #Comparing the original and bootstrap partitions
        for idx, node_id in enumerate(root_partitions) :
            if partitons_map[node_id] == btstr_partitions_map[node_id] :
                matching_partitions[idx] += 1

    #Bootstrap confidences
    btstr_percentages = [matches / 100.0 for matches in matching_partitions]
    #Writing the bootstrap confidences to the bootstrap file
    with open(bootstrap_file_name, 'w') as file:
        for percent in btstr_percentages:
            file.write(str(percent) + '\n')
