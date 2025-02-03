
OPEN_PENALTY = 4
EXTEND_PENALTY = 1

SCORE = [[5, -1, -4, -4, -2], [-1, 5, -4, -4, -2], [-4, -4, 5, -1, -2], [-4, -4, -1, 5, -2], [-2, -2, -2, -2, -1]]



def alignTwoSequences(seq1, seq2):

    # rows & cols +1 to account for boundary conditions
    rows = len(string1) + 1
    cols = len(string2) + 1

    # create Match matrix
    M = [[0] * cols for _ in range(rows)]

    # create Insertion matrix
    I = [[0] * cols for _ in range(rows)]

    # create Deletion matrix
    D = [[0] * cols for _ in range(rows)]

    # Base case visualization, < is negative infinity

    # M Matrix
    # * * T A G
    # * 0 < < <
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
     for i in rows:
        for j in cols:




    

