# Ashley E. Deleon 

import itertools
import numpy as np
from spectrumkernel import compute_spectrum_kernel


#reads the FASTA file and assign labels (exon = 1, intron = 0)
def read_fasta_with_labels(fasta_files, label=None):
    sequences = []
    labels = []
    for filename in fasta_files:
        with open(filename, "r") as f:
            seq = ""
            seq_label = label
            for line in f:
                if line.startswith(">"):  #header line, extract label if available
                    if seq:
                        sequences.append(seq)
                        labels.append(seq_label)
                        seq = ""
                    if label is None:  #get label from defline if not provided
                        parts = line.strip().split("/")
                        for part in parts:
                            if part.startswith("class="):
                                seq_label = 1 if part.split("=")[1] == "exon" else 0 if part.split("=")[1] == "intron" else -1  # -1 for intergenic or unknown
                else:
                    seq += line.strip()
            if seq:
                sequences.append(seq)
                labels.append(seq_label)
    return sequences, labels

#generate all k-mers of length k using ACGT
def generate_kmers(k):
    bases = ["A", "C", "G", "T"]
    return ["".join(p) for p in itertools.product(bases, repeat=k)]

# compute the feature matrix -- spectrum kernel feature vectors -- each sequence into num vector based on  k-mer count
def compute_feature_matrix(sequences, k):
    kmers = generate_kmers(k)  #all possible k-mers of length k
    kmer_indices = {kmer: i for i, kmer in enumerate(kmers)}  #map each k-mer to its index in the feature vector
    feature_matrix = np.zeros((len(sequences), len(kmers)), dtype=int)

    for i, seq in enumerate(sequences):
        kmer_counts = {kmer: 0 for kmer in kmers}  #initialize k-mer counts for the current sequence

        for j in range(len(seq) - k + 1):
            kmer = seq[j:j+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1

        #assing k-mer counts to the feature matrix
        for kmer, count in kmer_counts.items():
            feature_matrix[i, kmer_indices[kmer]] = count

    return feature_matrix

#this computes the Spectrum Kernel -- similiarity matrx between seqs using dot product of their feature vectors
def compute_spectrum_kernel(feature_matrix):
    return np.dot(feature_matrix, feature_matrix.T)


#KNN classification using spectrum kernel similarity -- classifies test sequences based on their similarity to the knn in the training set
#Uses similarity scores as weights so closer neighbors have more influence

def knn_classify(train_features, train_labels, test_features, k):
    num_test = test_features.shape[0]
    predictions = []

    for i in range(num_test):
        similarities = np.dot(train_features, test_features[i])  #similarity scores
        nearest_neighbors = np.argsort(similarities)[-k:]  # k closest neighbors (highest similarities)
        
        neighbor_labels = np.array([train_labels[j] for j in nearest_neighbors]) 
        neighbor_similarities = similarities[nearest_neighbors]  #similarity scores for neighbors
        
        # weighted scores (higher similarity → higher weight)
        exon_score = np.sum(neighbor_similarities[neighbor_labels == 1])
        intron_score = np.sum(neighbor_similarities[neighbor_labels == 0])
        
        #class with the highest weighted similarity sum
        if exon_score > intron_score:
            prediction = 1  # Exon
        elif intron_score > exon_score:
            prediction = 0  # Intron
        else:
            prediction = np.random.choice([0, 1])  #resolve ties randomly
        
        predictions.append(prediction)

    return np.array(predictions)

# KNN and accuracy  -- accuracy of the KNN classifier for given parameters 
def evaluate_knn(train_seqs, train_labels, test_seqs, test_labels, k, k_mer_size):
    train_features = compute_feature_matrix(train_seqs, k_mer_size)
    test_features = compute_feature_matrix(test_seqs, k_mer_size)
    
    predictions = knn_classify(train_features, train_labels, test_features, k)
    accuracy = np.mean(predictions == test_labels)
    return accuracy

if __name__ == "__main__":
    #file paths for each sample size; paths are dictionaries mapping sample sizes to file paths for exon and intron training data
    training_files_exons = {
        10: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-exons10.fasta",
        30: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-exons30.fasta",
        100: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-exons100.fasta"
    }
    training_files_introns = {
        10: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-introns10.fasta",
        30: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-introns30.fasta",
        100: "/hpc/group/coursess25/CS561-CS260/DATA/project2/train-introns100.fasta"
    }

    #the path to the test data
    test_file = "/hpc/group/coursess25/CS561-CS260/DATA/project2/test.fasta"

    #reads test data with labels from the test file
    test_sequences, test_labels = read_fasta_with_labels([test_file])  # Extract labels from defline

    #parameter values 
    k_mer_sizes = [2, 4, 6, 8]  #k-mer sizes
    knn_values = [1, 3, 5, 7]  #k values -- lists the num of nearest neighbors to consider for classification

    print("Sample Size | K-mer Size | K | Accuracy")
    print("----------------------------------------")

    #iterate over the three sample sizes (10, 30, 100)
    for sample_size in [10, 30, 100]:
        #read the corresponding training data for the current sample size
        train_exons, exon_labels = read_fasta_with_labels([training_files_exons[sample_size]], label=1) #exon - 1
        train_introns, intron_labels = read_fasta_with_labels([training_files_introns[sample_size]], label=0) #intron - 0 

        #combine exon and intron training sets
        train_sequences = train_exons + train_introns
        train_labels = exon_labels + intron_labels

        #run KNN for all combinations of k-mer size and k values
        for k_mer in k_mer_sizes:
            for k in knn_values:
                accuracy = evaluate_knn(train_sequences, train_labels, test_sequences, test_labels, k, k_mer)
                print(f"{sample_size:^12} | {k_mer:^10} | {k:^1} | {accuracy:.4f}")


        
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