import numpy as np
import itertools
from collections import Counter

import mismatchkernel as mk
import spectrumkernel as sk

def k_means_kernel(kernel_matrix, k, max_iters=100):
    num_sequences = kernel_matrix.shape[0]
    centers = np.random.choice(num_sequences, size=k, replace=False)
    labels = np.zeros(num_sequences)
    
    for iteration in range(max_iters):
        for i in range(num_sequences):
            similarities = kernel_matrix[i, centers]
            labels[i] = np.argmin(similarities)
        
        new_centers = []
        for cluster_id in range(k):
            cluster_indices = np.where(labels == cluster_id)[0]
            if len(cluster_indices) > 0:
                cluster_matrix = kernel_matrix[cluster_indices][:, cluster_indices]
                center = np.mean(cluster_matrix, axis=0)
                new_center_index = np.argmin(np.sum(center))
                new_centers.append(cluster_indices[new_center_index])
            else:
                new_centers.append(centers[cluster_id])
        centers = np.array(new_centers)
    
    return labels

def compute_class_proportions(sequences, labels, class_map):
    num_clusters = len(np.unique(labels))
    proportions = []
    
    for cluster_id in range(num_clusters):
        cluster_indices = np.where(labels == cluster_id)[0]
        cluster_classes = [class_map[sequences[i]] for i in cluster_indices]
        class_counts = Counter(cluster_classes)
        total = len(cluster_classes)
        
        cluster_proportions = {
            "exon": class_counts.get("exon", 0) / total,
            "intron": class_counts.get("intron", 0) / total,
            "intergenic": class_counts.get("intergenic", 0) / total
        }
        proportions.append(cluster_proportions)
    
    return proportions

# Assuming you have your sequence data and class labels
if __name__ == "__main__":
    fasta_file = "/path/to/your/kmeans.fasta"  # Specify your FASTA file
    k = 4  # K-mer size
    d = 1  # Maximum mismatches allowed
    
    sequences = sk.read_fasta(fasta_file)
    class_map = {seq: "exon" for seq in sequences}  # Dummy class map for example
    
    feature_matrix_spectrum = sk.compute_feature_matrix(sequences, k)
    kernel_matrix_spectrum = sk.compute_spectrum_kernel(feature_matrix_spectrum)
    
    # Apply k-means with the spectrum kernel
    k_clusters = 3
    labels_spectrum = k_means_kernel(kernel_matrix_spectrum, k_clusters)
    proportions_spectrum = compute_class_proportions(sequences, labels_spectrum, class_map)
    print("Proportions for Spectrum Kernel:", proportions_spectrum)
    
    feature_matrix_mismatch = mk.compute_feature_matrix_mismatch(sequences, k, d)
    kernel_matrix_mismatch = mk.compute_mismatch_kernel(feature_matrix_mismatch)
    
    # Apply k-means with the mismatch kernel
    labels_mismatch = k_means_kernel(kernel_matrix_mismatch, k_clusters)
    proportions_mismatch = compute_class_proportions(sequences, labels_mismatch, class_map)
    print("Proportions for Mismatch Kernel:", proportions_mismatch)
