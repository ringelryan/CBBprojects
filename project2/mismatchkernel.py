import itertools
import numpy as np

# Reads a FASTA file and returns a list of sequences
def read_fasta(filename):
    sequences = []
    with open(filename, "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):  # Header line, ignore
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)  # Add last sequence
    return sequences

# Generates all possible k-mers of length k using ACGT
def generate_kmers(k):
    bases = ["A", "C", "G", "T"]
    return ["".join(p) for p in itertools.product(bases, repeat=k)]

# Generates all k-mers that are within d mismatches from a given k-mer
def generate_mismatch_kmers(kmer, d):
    bases = ["A", "C", "G", "T"]
    mismatch_kmers = set()

    # Generate all possible k-mers with at most d mismatches
    def mutate_kmer(kmer, mismatches_left, index, current_kmer):
        if mismatches_left == 0 or index == len(kmer):
            mismatch_kmers.add("".join(current_kmer))
            return
        # Keep the original letter
        mutate_kmer(kmer, mismatches_left, index + 1, current_kmer)
        # Substitute with other bases
        for base in bases:
            if base != kmer[index]:  # Only change if different
                new_kmer = list(current_kmer)
                new_kmer[index] = base
                mutate_kmer(kmer, mismatches_left - 1, index + 1, new_kmer)

    mutate_kmer(kmer, d, 0, list(kmer))
    return mismatch_kmers

# Computes the feature matrix considering mismatches
def compute_feature_matrix_mismatch(sequences, k, d):
    kmers = generate_kmers(k)  # Generate all possible k-mers
    feature_matrix = np.zeros((len(sequences), len(kmers)), dtype=int)  # Initialize matrix

    kmer_to_index = {kmer: i for i, kmer in enumerate(kmers)}  # Map k-mers to indices

    for i, seq in enumerate(sequences):
        kmer_counts = {kmer: 0 for kmer in kmers}  # Initialize k-mer counts
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j + k]  # Extract k-mer
            if kmer in kmer_counts:  # Exact match
                kmer_counts[kmer] += 1
            # Consider mismatched k-mers
            for mismatch_kmer in generate_mismatch_kmers(kmer, d):
                if mismatch_kmer in kmer_counts:
                    kmer_counts[mismatch_kmer] += 1

        feature_matrix[i] = [kmer_counts[kmer] for kmer in kmers]  # Convert dictionary to vector

    return feature_matrix

# Computes the mismatch kernel (similarity matrix)
def compute_mismatch_kernel(feature_matrix):
    return np.dot(feature_matrix, feature_matrix.T)  # Kernel matrix using dot product

if __name__ == "__main__":
    fasta_file = "/hpc/group/coursess25/CS561-CS260/DATA/project2/kmeans.fasta"  
    k = 4  # k-mer size
    d = 1  # Maximum allowed mismatches

    sequences = read_fasta(fasta_file)  
    feature_matrix = compute_feature_matrix_mismatch(sequences, k, d)  
    kernel_matrix = compute_mismatch_kernel(feature_matrix)  

    print("Feature Matrix Shape:", feature_matrix.shape)
    print("Kernel Matrix Shape:", kernel_matrix.shape)
    print(kernel_matrix)

    # TESTING
    test_sequences = [
        "ATTGGCAA",  # k=2
        "GTTACAGT",  # k=2
    ]
    
    k = 2  # Small k-mer size for testing
    d = 1  # Allow 1 mismatch

    feature_matrix = compute_feature_matrix_mismatch(test_sequences, k, d)
    print("Feature Matrix:\n", feature_matrix)

    kernel_matrix = compute_mismatch_kernel(feature_matrix)
    print("Kernel Matrix:\n", kernel_matrix)

    for seq in test_sequences:
        print([seq[i:i+k] for i in range(len(seq) - k + 1)])

    print("Similarity between seq1 and seq2:", kernel_matrix[0, 1])