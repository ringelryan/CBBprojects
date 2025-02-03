
OPEN_PENALTY = -4
EXTEND_PENALTY = -1
NEG_INFINITY = -1000 # big enough ?

SCORE = [[5, -1, -4, -4, -2], [-1, 5, -4, -4, -2], [-4, -4, 5, -1, -2], [-4, -4, -1, 5, -2], [-2, -2, -2, -2, -1]]



def alignTwoSequences(seq1, seq2):

    # rows & cols +1 to account for boundary conditions
    rows = len(string1) + 1
    cols = len(string2) + 1

    # create Match matrix: Score of best alignment ending in match
    M = [[0] * cols for _ in range(rows)]

    # create Insertion matrix: Score of best alignment ending in Insertion
    I = [[0] * cols for _ in range(rows)]

    # create Deletion matrix: Score of best alignment ending in Deletion
    D = [[0] * cols for _ in range(rows)]

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
    for i in rows:
        M[i][j] = NEG_INFINITY
        I[i][j] = OPEN_PENALTY + (i-1) * EXTEND_PENALTY
        D[i][j] = NEG_INFINITY

    # first row
    i = 0
    for j in cols:
        M[i][j] = NEG_INFINITY
        I[i][j] = NEG_INFINITY
        D[i][j] = OPEN_PENALTY + (j-1) * EXTEND_PENALTY

    # corners
    M[0][0] = NEG_INFINITY
    I[0][0] = NEG_INFINITY
    D[0][0] = NEG_INFINITY

    M[1][1] = # Score of match/mismatch


        




    

