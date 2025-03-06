import itertools
import numpy as np

#reads a file and rets a lst of seqs
def read_fasta(filename):
    sequences = []
    with open(filename, "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):  #header ignore
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)  #add last seq
    return sequences


#generates all possible k-mers of length k using ACGT
def generate_kmers(k):
    bases = ["A", "C", "G", "T"] #DNA seqs made from ACGT 
    return ["".join(p) for p in itertools.product(bases, repeat=k)]


#computes the feature matrix where each seq is represented by k-mer counts
def compute_feature_matrix(sequences, k):
    kmers = generate_kmers(k)  #all possible k-mers
    feature_matrix = np.zeros((len(sequences), len(kmers)), dtype=int)  #initilize feature matrix

    for i, seq in enumerate(sequences):
        kmer_counts = {kmer: 0 for kmer in kmers}  #initilize k-mer counts 
        for j in range(len(seq) - k + 1):  
            kmer = seq[j:j+k]  #extract k-mer
            if kmer in kmer_counts:  #increment count if valid
                kmer_counts[kmer] += 1
        
        feature_matrix[i] = [kmer_counts[kmer] for kmer in kmers]  #convert dict to vector

    return feature_matrix

#computes the spectrum kernel --<> dot product similarity matrix from the feature matrix
def compute_spectrum_kernel(feature_matrix):
    return np.dot(feature_matrix, feature_matrix.T)  #kernel matrix with dot product

if __name__ == "__main__":
    fasta_file = "/hpc/group/coursess25/CS561-CS260/DATA/project2/kmeans.fasta" 
    k = 4  # k-mer size

    sequences = read_fasta(fasta_file)  
    feature_matrix = compute_feature_matrix(sequences, k)  
    kernel_matrix = compute_spectrum_kernel(feature_matrix)  

    print("Feature Matrix Shape:", feature_matrix.shape)
    print("Kernel Matrix Shape:", kernel_matrix.shape)
    print(kernel_matrix)

print("test")
#TESTING 
test_sequences = [
    "ATTGGCAA",  # k=2' AT, TT, TG, GG, GC, CA, AA
    "GTTACAGT", #k=2 GT, TT, TA, AC, CA, AG, GT
]

k = 2  # k-mer size; 4^2

#feature matrix
feature_matrix = compute_feature_matrix(test_sequences, k)
print("Feature Matrix:\n", feature_matrix)

#spectrum kernel
kernel_matrix = compute_spectrum_kernel(feature_matrix)
print("Kernel Matrix:\n", kernel_matrix)

for seq in test_sequences:
    print([seq[i:i+k] for i in range(len(seq) - k + 1)])