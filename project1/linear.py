
# Ashley E. Deleon - Linear Scoring 

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

#read substitution matrix; rets a dict 

def read_substitution_matrix(filename):
    #Reads the substitution matrix file and returns a dict scores
    matrix = {}
    
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip()]  #remove empty lines
        headers = lines[0].split()  #column headers

        for line in lines[1:]:  #rows
            row = line.split()
            row_label = row[0]  #first entry is the row label
            scores = list(map(int, row[1:]))  #convert score to ints
            
            for col_label, score in zip(headers, scores):
                matrix[(row_label, col_label)] = score  #score in dict

    #(A, B) - (B, A)
    for (a, b) in list(matrix.keys()):
        if (b, a) not in matrix:  #if reverse missing, add it
            matrix[(b, a)] = matrix[(a, b)]

    print("Substitution Matrix:")
    for key, value in matrix.items():
        print(f"{key}: {value}") 

    return matrix


substitution_matrix = read_substitution_matrix("/hpc/group/coursess25/CS561-CS260/DATA/project1/matrix.txt")




# needleman-wunsch (linear scoring) global alignment; rets aligned seqs and optimal alignment scores 
def needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty):
    n = len(seq1)  
    m = len(seq2)  

    #start of scoreing and traceback (size (n+1) x (m+1))
    scoring_matrix = [[0] * (m+1) for _ in range(n+1)]
    traceback_matrix = [[""] * (m+1) for _ in range(n+1)]

    #up and left incur gap penalty 
    for i in range(1, n+1):
        scoring_matrix[i][0] = -i * gap_penalty
        traceback_matrix[i][0] = 'U'  #up, gap in seq2

    for j in range(1, m+1):
        scoring_matrix[0][j] = -j * gap_penalty
        traceback_matrix[0][j] = 'L'  #left, gap in seq1

    #filling in the matrices; iterate each cell and fill optimal score  
    for i in range(1, n+1):
        for j in range(1, m+1):
            #calcs scores for mis/match, deletion, and insertion
            match = scoring_matrix[i-1][j-1] + substitution_matrix[(seq1[i-1], seq2[j-1])] 
            delete = scoring_matrix[i-1][j] - gap_penalty #score if we delete from seq1 (gap in seq2)
            insert = scoring_matrix[i][j-1] - gap_penalty #gap in seq1

            #get the best score from the three
            max_score = max(match, delete, insert)
            scoring_matrix[i][j] = max_score
          

            #track the traceback direction
            if max_score == match:
                traceback_matrix[i][j] = 'D'  #diagonal (mis/match)
                
            elif max_score == delete:
                traceback_matrix[i][j] = 'U'  #up
                
            else:  
                traceback_matrix[i][j] = 'L'  #left 
                

    #traceback to construct alignment 
    align1, align2 = "", "" #to store aligned seqs 
    i, j = n, m  #bottom right

   # while i >= 0 or j >= 0:

       # if traceback_matrix[i][j] == 'D': #aligning in both; add chars to the front
        #    align1 = seq1[i-1] + align1
         #   align2 = seq2[j-1] + align2
         #   i -= 1
         #   j -= 1
        #elif traceback_matrix[i][j] == 'U': #aligning char from seq1 to gap in seq2, add char in seq1, gap to seq2
        #    align1 = seq1[i] + align1
        #    align2 = "-" + align2
        #    i -= 1
       # else:  #left
        #    align1 = "-" + align1 #aligning gap in seq1 with char in seq2; add gap to seq1, char from deq2
        #    align2 = seq2[j] + align2
        #    j -= 1
    #optimal in the bottom right 
    #return align1, align2, scoring_matrix[-1][-1]


    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 'D': #mis/match
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i][j] == 'U': #gap in se1
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1 #gap in seq1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, scoring_matrix[n][m] #aligned seq and score 


#either pass og ines or take out the gaps 
#calcs alignment stats (matches, percent, indels, mean indel, alignment length, score); for the output 
def calc_sum_stats(align1, align2, score, ): #aligned seq and score 
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-' and b != '-')  #counts num of pos where chars are matched

    indels = sum(1 for a, b in zip(align1, align2) if a == "-" or b == "-") #counts num of gaps in either
    
    indel_lengths = []
    current_indel = 0
    in_indel = False

    #indel lengths 
    for a, b in zip(align1, align2):
        if a == b and a != '-' and b != '-':
            if in_indel:
                indel_lengths.append(current_indel) 
                current_indel = 0
                in_indel = False
        elif a == "-" or b == "-":
            current_indel += 1
            in_indel = True

    if in_indel:
        indel_lengths.append(current_indel)

    #calcs avg length of indels; no indels ret 0 
    mean_indel_length = sum(indel_lengths) / len(indel_lengths) if indel_lengths else 0 

    #remove the gaps
    seq1_length = len(align1.replace("-", ""))
    seq2_length = len(align2.replace("-", ""))

    average_length = (seq1_length + seq2_length) / 2

    #calcs the percent identity
    percent_identity = (matches / average_length) * 100

 
    #rets dict
    return {
        "matches": matches,
        "percent_identity": percent_identity,
        "indels": indels,
        "mean_indel_length": mean_indel_length,
        "alignment_length": len(align1),
        "score": score
    }

#formatting the alignment output 
def format_alignment_output(alignment_num, seq1_name, seq2_name, align1, align2, stats):

    #create match line | for match; * mismatch, " " gap; visual rep 

    match_line = "".join("|" if a == b and a != '-' and b != '-' else 
                     " " if a == '-' or b == '-' else 
                     "*" for a, b in zip(align1, align2))   
     
    #split alignemnts into 60 char
    width = 60
    align1_chunks = [align1[i:i+width] for i in range(0, len(align1), width)]
    match_chunks = [match_line[i:i+width] for i in range(0, len(match_line), width)]
    align2_chunks = [align2[i:i+width] for i in range(0, len(align2), width)]


    
    #directly formatting the output string; alignment num, seq name, stats 
    formatted = f"Alignment #{alignment_num}:\n"
    formatted += f"Sequence #1: {seq1_name}\n"
    formatted += f"Sequence #2: {seq2_name}\n"
    formatted += f"Matches: {stats['matches']}\n"
    formatted += f"Percent identity: {stats['percent_identity']:.2f}%\n"
    formatted += f"Indels: number={stats['indels']}, mean length={stats['mean_indel_length']:.1f}\n"
    formatted += f"Alignment Length: {stats['alignment_length']}\n"
    formatted += f"Score={stats['score']}\n\n"
    
    #add aligned seqs to string, iterates over chuncks in both and match line
    for a, m, b in zip(align1_chunks, match_chunks, align2_chunks):
        #chunck of first aligned seq, corresponding chunck of match line, corresponding second aligned seq
        formatted += f"{a}\n{m}\n{b}\n\n"

    #print(align1)
    #print(align2)
    return formatted

#sequences
close_headers1, close_seqs1 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-first.fasta")
close_headers2, close_seqs2 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/close-second.fasta")
distant_headers1, distant_seqs1 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-first.fasta")
distant_headers2, distant_seqs2 = read_fasta("/hpc/group/coursess25/CS561-CS260/DATA/project1/distant-second.fasta")

#substitution matrix
substitution_matrix = read_substitution_matrix("/hpc/group/coursess25/CS561-CS260/DATA/project1/matrix.txt")

#func aligns seq and writes the formatted output into a file; handles mismatched sequence counts
def align_and_write_output(output_file, headers1, seqs1, headers2, seqs2):
    with open(output_file, "w") as out: #to close it after
        for i in range(len(seqs1)):
            align1, align2, score = needleman_wunsch(seqs1[i], seqs2[i], substitution_matrix, 1) #aligned seqs and the score 
            stats = calc_sum_stats(align1, align2, score)
            out.write(format_alignment_output(i + 1, headers1[i], headers2[i], align1, align2, stats))


#output files
align_and_write_output("linear_close_output.txt", close_headers1, close_seqs1, close_headers2, close_seqs2)
align_and_write_output("linear_distant_output.txt", distant_headers1, distant_seqs1, distant_headers2, distant_seqs2)



#testing matches 
seq1 = "AGTCNAC"
seq2 = "AGTCNTA"
gap_penalty = 1

align1, align2, score = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
stats = calc_sum_stats(align1, align2, score)
output = format_alignment_output(1, "test1", "test2", align1, align2, stats)
print(output)

#testing matches 
seq1 = "TAGA"
seq2 = "AGAN"
gap_penalty = 1

align1, align2, score = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
stats = calc_sum_stats(align1, align2, score)
output = format_alignment_output(1, "test3", "test4", align1, align2, stats)
print(output)


#testing matxches 
seq1 = "AAAGNAGN"
seq2 = "AAATAATG"
gap_penalty = 1

align1, align2, score = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)
stats = calc_sum_stats(align1, align2, score)
output = format_alignment_output(1, "test3", "test4", align1, align2, stats)
print(output)


    

  