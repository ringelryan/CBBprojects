# Ryan P. Ringel - Affine Scoring

OPEN_PENALTY = -4
EXTEND_PENALTY = -1
NEG_INFINITY = -1000 # big enough ?

import re

SCORE = [[5, -1, -4, -4, -2], [-1, 5, -4, -4, -2], [-4, -4, 5, -1, -2], [-4, -4, -1, 5, -2], [-2, -2, -2, -2, -1]]

#read seq from the fasta files
def read_fasta(filename):
    sequences = []
    headers = []
    with open(filename, "r") as file:
        sequence = ""
        header = None
        for line in file:
            if line.startswith(">"):  #new seq header
                if sequence:  
                    #if completed seq when new header; add it and header
                    sequences.append(sequence)
                    headers.append(header)
                header = line.strip()[1:]  #header without ">"
                sequence = ""  
            else:
                sequence += line.strip()  
        if sequence:  
            sequences.append(sequence) #final seq: add it
            headers.append(header)
    return headers, sequences

def alignTwoSequences(seq1, seq2):

    # rows & cols +1 to account for boundary conditions
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # create Match matrix: Score of best alignment ending in match
    M = [[0] * cols for _ in range(rows)]
    traceback_M = [[""] * cols for _ in range(rows)]

    # create Insertion matrix: Score of best alignment ending in Insertion
    I = [[0] * cols for _ in range(rows)]
    traceback_I = [[""] * cols for _ in range(rows)]

    # create Deletion matrix: Score of best alignment ending in Deletion
    D = [[0] * cols for _ in range(rows)]
    traceback_D = [[""] * cols for _ in range(rows)]

    # Base case visualization, < is negative infinity

    # M Matrix
    # * * T A G
    # * < < < <
    # C <-4
    # A <
    # T <
    #______________
    # I Matrix:  - - -
    #            C A T

    # * * T A G
    # * < < < <
    # C-4
    # A-5
    # T-6
    #______________
    # D Matrix:  T A G
    #            - - -

    # * * T A G
    # * <-4-5-6
    # C <
    # A <
    # T <


    # Initialize base cases for matrices 
    # __________________________________

    # first column
    j = 0
    for i in range(rows):
        M[i][j] = NEG_INFINITY
        I[i][j] = NEG_INFINITY
        D[i][j] = OPEN_PENALTY + (j-1) * EXTEND_PENALTY

        traceback_M[i][j] = "MDU"
        traceback_D[i][j] = "DDU"


    # first row
    i = 0
    for j in range(cols):
        M[i][j] = NEG_INFINITY
        I[i][j] = OPEN_PENALTY + (i-1) * EXTEND_PENALTY
        D[i][j] = NEG_INFINITY

        traceback_M[i][j] = "MIL"
        traceback_I[i][j] = "IIL"

    # corners
    M[0][0] = 0  # allows for M[1][1] to be score of first match/mismatch
    I[0][0] = NEG_INFINITY
    D[0][0] = NEG_INFINITY

    traceback_M[0][0] = ""
    traceback_I[0][0] = ""
    traceback_D[0][0] = ""


    # Base Cases for traceback matrices: 
    # Top of D should trace to starting with deletion
    # Side of I should trace to starting with insertion
    # Top of M should trace to starting with deletion
    # Side of M should trace to starting with insertion



    #________________________________________________


    # Iterate through and fill out matrices
    for i in range(1, rows):
        for j in range(1, cols):
            si = getLetterIndex(seq1[i-1])
            sj = getLetterIndex(seq2[j-1])
            s = SCORE[si][sj]

            # Calculate values in matrices based on formulas
            M[i][j] = max((M[i-1][j-1] + s), (I[i-1][j-1] + s), (D[i-1][j-1] + s))
            I[i][j] = max((M[i][j-1] + OPEN_PENALTY), (I[i][j-1] + EXTEND_PENALTY))
            D[i][j] = max((M[i-1][j] + OPEN_PENALTY), (D[i-1][j] + EXTEND_PENALTY))

            # store values in traceback matrices based on max value: First letter is which matrix, second letter is diag, up, or left

            # traceback for M matrix
            if (M[i-1][j-1] + s) >= (I[i-1][j-1] + s) and (M[i-1][j-1] + s) >= (D[i-1][j-1] + s):
                traceback_M[i][j] = "MMD" # traceback points to M matrix on diagonal
            elif (I[i-1][j-1] + s) >= (M[i-1][j-1] + s) and (I[i-1][j-1] + s) >= (D[i-1][j-1] + s):
                traceback_M[i][j] = "MID" # traceback came from  I matrix on diagonal
            elif (D[i-1][j-1] + s) >= (M[i-1][j-1] + s) and (D[i-1][j-1] + s) >= (I[i-1][j-1] + s):
                traceback_M[i][j] = "MDD" # traceback came from D matrix on diagonal

            # traceback for I matrix
            if (M[i][j-1] + OPEN_PENALTY) >= (I[i][j-1] + EXTEND_PENALTY):
                traceback_I[i][j] = "IML" # traceback came from M matrix from left
            elif (I[i][j-1] + EXTEND_PENALTY) >= (M[i][j-1] + OPEN_PENALTY):
                traceback_I[i][j] = "IIL" # traceback came from I matrix from left

            # traceback for D matrix
            if (M[i-1][j] + OPEN_PENALTY) >= (D[i-1][j] + EXTEND_PENALTY):
                traceback_D[i][j] = "DMU" # traceback came from M matrix from Up
            elif (D[i-1][j] + EXTEND_PENALTY) >= (M[i-1][j] + OPEN_PENALTY):
                traceback_D[i][j] = "DDU" # traceback came from D matrix from Up


    #After Done, print out best scores
    print("Best M: " + str(M[rows-1][cols-1]))
    print("Best I: " + str(I[rows-1][cols-1]))
    print("Best D: " + str(D[rows-1][cols-1]))

    # print("M: " + str(M))
    # print("I: " + str(I))
    # print("D: " + str(D))

    # print("trace_M" + str(traceback_M))
    # print("trace_I" + str(traceback_I))
    # print("trace_D" + str(traceback_D))

    # call traceback on correct matrix
    if (M[rows-1][cols-1]) >= (I[rows-1][cols-1]) and (M[rows-1][cols-1]) >= (D[rows-1][cols-1]):
        print("M highest")
        printTraceback(seq1, seq2, traceback_M, traceback_I, traceback_D, traceback_M, M[rows-1][cols-1])
    elif (I[rows-1][cols-1]) >= (M[rows-1][cols-1]) and (I[rows-1][cols-1]) >= (D[rows-1][cols-1]):
        print("I highest")
        printTraceback(seq1, seq2, traceback_M, traceback_I, traceback_D, traceback_I, I[rows-1][cols-1])
    elif (D[rows-1][cols-1]) >= (M[rows-1][cols-1]) and (D[rows-1][cols-1]) >= (I[rows-1][cols-1]):
        print("D highest")
        printTraceback(seq1, seq2, traceback_M, traceback_I, traceback_D, traceback_D, D[rows-1][cols-1])

# Print alignment of the two sequences, startMatrix should be traceback matrix of the matrix with highest alignment score in bottom right corner
def printTraceback(seq1, seq2, traceback_M, traceback_I, traceback_D, curMatrix, score):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # lines to print out; store as lists for appending, reverse at end
    line1 = []
    line2 = []
    line3 = []

    i = rows - 1
    j = cols - 1

    # get code for first traceback
    trace_code = curMatrix[i][j]

    # code structure: ex. "MMD"
    # 1st letter indicates M, I, or D at current i, j:
    # 2nd letter indicates which matrix you came from, either M, D, or I
    # 3rd letter indicates what direction you came from, either D (diagonal), U (up), or L (left)

    while trace_code != "":
        letter_one = trace_code[0]
        letter_two = trace_code[1]
        letter_three = trace_code[2]

        # sequences are 0 indexed, retrieve using i-1/ j-1

        if letter_one == 'M':
            #print("match: " + trace_code)
            # i, j is a Match/ Mismatch: 
            line1.append(seq1[i-1])

            # Determine if a match or mismatch
            if seq1[i-1] == seq2[j-1]:
                line2.append('|')
            else:
                line2.append('*')

            line3.append(seq2[j-1])
        elif letter_one == 'I':
            #print("Insertion: " + trace_code)
            # i, j is insertion: line3 will have -
            line1.append('-')
            line2.append(' ')
            line3.append(seq2[j-1])
        elif letter_one == 'D':
            #print("deletion: " + trace_code)
            # i, j is deletion: line1 will have -
            line1.append(seq1[i-1])
            line2.append(' ')
            line3.append('-')

        # Update i and j based on direction
        if letter_three == 'D':
            # diagonal
            i = i-1
            j = j-1
        elif letter_three == 'U':
            # up
            i = i-1
        elif letter_three == 'L':
            # left
            j = j-1

        # retrieve new trace_code from curMatrix and direction
        if letter_two == 'M':
            
            trace_code = traceback_M[i][j]
            #print("letter 2 is M, new code is: "+trace_code)
        elif letter_two == 'I':
            trace_code = traceback_I[i][j]
            # print("traceback_I: " + str(traceback_I))
            #print("letter 2 is I, new code is: "+trace_code)
        elif letter_two == 'D':
            trace_code = traceback_D[i][j]
            #print("letter 2 is D, new code is: "+trace_code)

        


    line1_string = "".join(reversed(line1))
    line2_string = "".join(reversed(line2))
    line3_string = "".join(reversed(line3))

    printStats(seq1, seq2, line2_string, score)

    # divide into chunks of 60 and print

    chunk_size = 60
    for i in range(0, len(line1_string), chunk_size):
        print(line1_string[i:i + chunk_size])
        print(line2_string[i:i + chunk_size])
        print(line3_string[i:i + chunk_size])
        print()

    # print(line1_string)
    # print(line2_string)
    # print(line3_string)




def printStats(seq1, seq2, line2, score):
    # count number of |'s in line2
    matches = line2.count('|')
    
    space_segments = re.findall(r' +', line2)  # Matches one or more spaces

    # Count the number of segments
    num_segments = len(space_segments)

    # Calculate the mean length of space segments
    if num_segments > 0:
        mean_length = sum(len(segment) for segment in space_segments) / num_segments
    else:
        mean_length = 0  # Handle the case where there are no spaces

    print("Alignment #1:")
    print()
    print("Sequence #1: ")
    print("Sequence #2: ")
    print("Matches: "+str(matches))
    print("Percent Identity: " + str(100*(matches/len(line2))) + "%")
    print("Indels: number=" + str(num_segments) + ", mean length =" + str(mean_length))
    print("Alignment Length: " + str(len(line2)))
    print("Score=" + str(score))
    print()


# Take in one letter, return 0, 1, 2, 3, or 4
def getLetterIndex(letter):
    if letter == 'A':
        return 0
    elif letter == 'C':
        return 1
    elif letter == 'G':
        return 2
    elif letter == 'T':
        return 3
    elif letter == 'N':
        return 4
        


def main():
    # seq1 = "ATCCGATATGCGCGATATGGGGTACCCATAATTTAACCGAGAGCAGATAAGACACCCAGTATA"
    # seq2 = "ATGCAATATTCAGAGGGGCAAATACATAGACCAGCATTACAGGACATAATACCCCATTTAGAGACCTA"

    # Read Sequences from class folder

    close_headers1, close_seqs1 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-first.fasta")
    close_headers2, close_seqs2 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-second.fasta")
    distant_headers1, distant_seqs1 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-first.fasta")
    distant_headers2, distant_seqs2 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-second.fasta")

    # seq1 = close_seqs1[9]
    # seq2 = close_seqs2[9]


    seq1 = distant_seqs1[9]
    seq2 = distant_seqs2[9]

    alignTwoSequences(seq1, seq2)



if __name__ == "__main__":
    main()