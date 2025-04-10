import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models, callbacks
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import pandas as pd
import os

#to check if im actully in gpu
print("Num GPUs Available:", len(tf.config.list_physical_devices('GPU')))
print(tf.config.list_physical_devices('GPU'))

# DNA base encoding
base_dict = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'N': [0, 0, 0, 0]
}

results = []  #store evaluation results

def one_hot_encode(seq):
    return np.array([base_dict.get(base.upper(), [0, 0, 0, 0]) for base in seq.strip()])

#loads sequences and their labels 
def load_fasta_and_labels(fasta_path):
    sequences = []
    labels = []
    current_seq = []
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(one_hot_encode(''.join(current_seq))) #encode previous seq
                    current_seq = []
                label = int(line.strip().split('=')[-1])  # Extract class label from defline
                labels.append(label)
            else:
                current_seq.append(line.strip()) #append sequence data
        
        if current_seq:
            sequences.append(one_hot_encode(''.join(current_seq)))  # last sequence

    #pad all seqs to the length of the longest sequence

    max_len = max(len(seq) for seq in sequences)
    padded_seqs = np.zeros((len(sequences), max_len, 4))
    for i, seq in enumerate(sequences):
        padded_seqs[i, :len(seq), :] = seq

    return padded_seqs, np.array(labels)

#Builds CNN MODEL 
def get_cnn_model(input_shape):
    model = models.Sequential([
        layers.Input(shape=input_shape),  #input takes the one hot encoded DNA sequence
        layers.Conv1D(64, 7, activation='relu', padding='same'),  # 1D
        layers.BatchNormalization(),  
        layers.MaxPooling1D(2),  #reduce dimensionality by a factor of 2
        layers.Conv1D(128, 5, activation='relu', padding='same'),  # 2D
        layers.BatchNormalization(),  
        layers.GlobalMaxPooling1D(),  # reduce the entire sequence to a single vector
        layers.Dense(128, activation='relu'),  #connected layer with ReLU activation
        layers.Dropout(0.5),  #dropoutfor regularization to prevent overfitting
        layers.Dense(1, activation='sigmoid')  #output layer with sigmoid activation -> binary classification
    ])
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])  
    return model

#Builds Transformer MODEL
def get_transformer_model(input_shape):
    inputs = layers.Input(shape=input_shape)  #input layer
    x = layers.Dense(64)(inputs)  #dense layer to transform the input
    x = layers.MultiHeadAttention(num_heads=4, key_dim=64)(x, x)  #to learn contextual relationships
    x = layers.LayerNormalization()(x)  #to stabilize training
    x = layers.GlobalAveragePooling1D()(x)  
    x = layers.Dense(128, activation='relu')(x)  # Fully connected layer with ReLU
    x = layers.Dropout(0.4)(x)  # Dropout for regularization
    outputs = layers.Dense(1, activation='sigmoid')(x)  # Output layer with sigmoid activation
    model = models.Model(inputs, outputs)  # Create the model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])  #compile Adam optimizer and binary crossentropy loss
    return model

# Function to train and evaluate the model
def train_and_evaluate(model, X_train, Y_train, X_valid, Y_valid, X_test, Y_test, name="Model"):
    #early stopping callback to prevent overfitting
    es = callbacks.EarlyStopping(patience=10, monitor="val_loss", restore_best_weights=True)
    model.fit(X_train, Y_train, 
              validation_data=(X_valid, Y_valid),
              batch_size=128, epochs=100, 
              callbacks=[es], verbose=1)
    
    #predicts on test data and calculate accuracy and AUC
    preds = model.predict(X_test, batch_size=128).flatten()
    acc = np.mean((preds > 0.5) == Y_test)  # Accuracy: compare predictions to actual labels
    auc = roc_auc_score(Y_test, preds) # AUC score

    # ROC Plot
    fpr, tpr, _ = roc_curve(Y_test, preds)
    plt.figure()
    plt.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve - {name} ({SIM})')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"roc_{SIM}_{name}.png", bbox_inches='tight')
    plt.close()

    print(f"{name} Test Accuracy: {acc:.4f}, AUC: {auc:.4f}")

    #store results in the results list
    results.append({
        "Model": name,
        "Dataset": SIM,
        "Accuracy": acc,
        "AUC": auc
    })

# === Main ===
for SIM in ["sim1", "sim2", "sim6", "sim7"]:
    print(f"\n=== Running on dataset: {SIM} ===")
    path = f"/hpc/group/coursess25/CS561-CS260/DATA/project3/{SIM}"

    X_train, Y_train = load_fasta_and_labels(os.path.join(path, "train.fasta"))
    X_valid, Y_valid = load_fasta_and_labels(os.path.join(path, "validation.fasta"))
    X_test, Y_test = load_fasta_and_labels(os.path.join(path, "test.fasta"))

    print(f"Data loaded: {X_train.shape}, Labels: {np.unique(Y_train, return_counts=True)}")

    #CNN and Transformer models and evaluate them
    cnn_model = get_cnn_model(X_train.shape[1:])
    train_and_evaluate(cnn_model, X_train, Y_train, X_valid, Y_valid, X_test, Y_test, name="CNN")

    transformer_model = get_transformer_model(X_train.shape[1:])
    train_and_evaluate(transformer_model, X_train, Y_train, X_valid, Y_valid, X_test, Y_test, name="Transformer")

#save all results
results_df = pd.DataFrame(results)
results_df.to_csv("results_all.csv", index=False)
print("All results saved to results_all.csv")
