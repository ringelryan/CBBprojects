# Ryan Ringel

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
    class_counts = []  # Store the total number of sequences per cluster

    for cluster_id in range(num_clusters):
        cluster_indices = np.where(labels == cluster_id)[0]
        cluster_classes = [class_map[sequences[i]] for i in cluster_indices]
        class_count = len(cluster_classes)
        class_counts.append(class_count)  # Store the number of sequences in this cluster

        class_freq = Counter(cluster_classes)  # Count occurrences of each class
        
        # Compute proportions
        cluster_proportions = {
            "exon": class_freq.get("/class=exon", 0) / class_count if class_count > 0 else 0,
            "intron": class_freq.get("/class=intron", 0) / class_count if class_count > 0 else 0,
            "intergenic": class_freq.get("/class=intergenic", 0) / class_count if class_count > 0 else 0
        }
        proportions.append(cluster_proportions)

    return proportions, class_counts  # Return both proportions and cluster sizes

def create_class_map(fasta_file):
    class_map = {}  # Dictionary to store sequence -> class mapping
    with open(fasta_file, "r") as f:
        current_class = None
        current_sequence = ""
        
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if current_sequence:  
                    class_map[current_sequence] = current_class  # Store previous sequence and its class
                
                parts = line[1:].split()  # Remove '>' and split header line
                current_class = parts[-1]  # Assume the last part of the header is the class label
                current_sequence = ""  # Reset sequence
            
            else:
                current_sequence += line  # Append sequence data
        
        if current_sequence:  # Store the last sequence
            class_map[current_sequence] = current_class
    
    return class_map

def print_cluster_proportions(kernel_name, k, k_clusters, proportions, class_counts):
    print(f"{kernel_name.upper()} (KMER={k}, CLUSTERS={k_clusters}):")

    for cluster_id, cluster_data in enumerate(proportions):  
        print(f"CLUSTER {cluster_id + 1}:")  

        for class_name, proportion in cluster_data.items():
            count = round(proportion * class_counts[cluster_id])  # Compute count from proportion
            print(f"{class_name} = {proportion:.2f} ({count})")

        # print()  # Blank line between clusters

    print()




if __name__ == "__main__":
    fasta_file = "/hpc/group/coursess25/CS561-CS260/DATA/project2/kmeans.fasta"  # Specify FASTA file
    k = 4  # K-mer size
    d = 1  # Maximum mismatches allowed
    
    sequences = sk.read_fasta(fasta_file)
    class_map = create_class_map(fasta_file)

    feature_matrix_spectrum = sk.compute_feature_matrix(sequences, k)
    kernel_matrix_spectrum = sk.compute_spectrum_kernel(feature_matrix_spectrum)
    
    # Apply k-means with the spectrum kernel
    k_clusters = 3
    labels_spectrum = k_means_kernel(kernel_matrix_spectrum, k_clusters)
    proportions_spectrum, class_counts_spectrum = compute_class_proportions(sequences, labels_spectrum, class_map)
    print_cluster_proportions("SPECTRUM KERNEL", k, k_clusters, proportions_spectrum, class_counts_spectrum)
    
    feature_matrix_mismatch = mk.compute_feature_matrix_mismatch(sequences, k, d)
    kernel_matrix_mismatch = mk.compute_mismatch_kernel(feature_matrix_mismatch)
    
    # Apply k-means with the mismatch kernel
    labels_mismatch = k_means_kernel(kernel_matrix_mismatch, k_clusters)
    proportions_mismatch, class_counts_mismatch = compute_class_proportions(sequences, labels_mismatch, class_map)
    print_cluster_proportions("MISMATCH KERNEL", k, k_clusters, proportions_mismatch, class_counts_mismatch)
