# read close sequences
def read_fasta_close(filename):
    sequences = []
    with open(filename, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):  
                if sequence:  
                    sequences.append(sequence)
                sequence = ""  
            else:
                sequence += line.strip()  
        if sequence:  
            sequences.append(sequence)
    return sequences

# parse close
close_first_sequences = read_fasta_close("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-first.fasta")
close_second_sequences = read_fasta_close("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-second.fasta")
print(close_first_sequences)
print(close_second_sequences)


print("hello")


# read distant sequences 
def read_fasta_distant(filename):
    sequences = []
    with open(filename, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):  
                if sequence:  
                    sequences.append(sequence)
                sequence = ""  
            else:
                sequence += line.strip()  
        if sequence:  
            sequences.append(sequence)
    return sequences

# parse distant
distant_first_sequences = read_fasta_distant("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-first.fasta")
distant_second_sequences = read_fasta_distant("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-second.fasta")
print(distant_first_sequences)
print(distant_second_sequences)


print("hello")

#read substitution matrix
def read_substitution_matrix(filename):
    matrix = {}
    with open(filename, "r") as file:
        lines = file.readlines()
        headers = lines[0].split()[1:] #need to to skip the top left corner 
        
        for line in lines[1:]:
            row = line.split()
            if not row:
                continue  # Skip empty lines
            
            row_label = row[0]  # First column is the row label
            scores = list(map(int, row[1:]))  # Remaining columns are scores

            for col_label, score in zip(headers, scores):
                matrix[(row_label, col_label)] = score  
        
    return matrix

substitution_matrix = read_substitution_matrix("/hpc/group/coursess25/CS561-CS260/DATA/project1/matrix.txt")
print(substitution_matrix)


# needleman-wunsch algorithm : global alignment aligning two seqs 
def needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty):
    n = len(seq1) + 1 #rows seq1 + gap 
    m = len(seq2) + 1 #cols seq2 + gap 

    scoring_matrix = [[0] * m for _ in range(n)]
    traceback_matrix = [[""] * m for _ in range(n)]

    #filling in the first row and column with gap
    for i in range(n):
        scoring_matrix[i][0] = -i * gap_penalty
        traceback_matrix[i][0] = 'U'  #up
    for j in range(m):
        scoring_matrix[0][j] = -j * gap_penalty
        traceback_matrix[0][j] = 'L'  #left

    #filling it in
    for i in range(1, n):
        for j in range(1, m):
            match = scoring_matrix[i-1][j-1] + substitution_matrix.get((seq1[i-1], seq2[j-1]), -gap_penalty)
            delete = scoring_matrix[i-1][j] - gap_penalty
            insert = scoring_matrix[i][j-1] - gap_penalty

            #get max score + update traceback
            scoring_matrix[i][j] = max(match, delete, insert)

            if scoring_matrix[i][j] == match:
                traceback_matrix[i][j] = 'D'  #(mis/match)
            elif scoring_matrix[i][j] == delete:
                traceback_matrix[i][j] = 'U'  
            else:
                traceback_matrix[i][j] = 'L'  
