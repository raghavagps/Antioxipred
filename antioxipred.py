import sys
import os
import numpy as np
import pandas as pd
import argparse
import re
import pickle
import subprocess
import datetime
import joblib
from catboost import CatBoostClassifier

parser = argparse.ArgumentParser()
parser.add_argument("--file", "-f", type = str, help='Input of peptide sequences in FASTA format or single sequence per line in single letter code', required = True)
parser.add_argument("--output", "-o", type =str, help = 'Path to output folder', required = True)
parser.add_argument("--model", "-m", type =int, help = 'Model_selection: (1) for OHE profile + PAAC based catboost model, (2) For OHE profile + PAAC based catboost model integrated with BLAST', default = 1)
parser.add_argument("--threshold", "-t", type = float, help = 'Threshold for classification. Value ranges from 0 to 1', default =0.51)

args = parser.parse_args()
fasta_loc = args.file
output_loc = args.output
model_chosen = args.model
threshold = args.threshold
directory_path = os.path.dirname(os.path.abspath(__file__))


# Read Fasta file or plain sequences
def read_fasta(file_path):
    sequences = []
    headers = []
    with open(file_path, 'r') as file:
        lines = file.readlines()

    if any(line.startswith('>') for line in lines):
        sequence =''
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
                headers.append(line[1:])    
            else:
                sequence += line
        if  sequence:
            sequences.append(sequence)
    else:
        for i, line in enumerate (lines):
            seq =line.strip()
            if seq:
                sequences.append(seq)
                headers.append(f"Seq_{i+1}")
    return headers, sequences

# padding the fasta sequences
def padding_sequence(sequences, padded_length=20):
    padded_sequences = []
    for seq in sequences:
        if len(seq)> padded_length:
            padded_seq = seq[:padded_length]
        else:
            padded_seq = seq.ljust(padded_length, 'X')
        padded_sequences.append(padded_seq)
    return padded_sequences

# Funtion to one hot encode the sequences.
def one_hot_encoding_profile(seq, padded_length=20):
    aa_list = list("ACDEFGHIKLMNPQRSTVWYX")
    aa_index = {aa: i for i, aa in enumerate(aa_list)}
    profile = np.zeros((padded_length, len(aa_index)), dtype=int)

    for pos, aa in enumerate(seq):
        profile[pos, aa_index[aa]] = 1

    return profile.flatten()
# Generate feature column names
def feature_column_names(padded_length =20):
    aa_list = list("ACDEFGHIKLMNPQRSTVWYX")
    column_names = []
    for i in range(1, padded_length + 1):
        for amino in aa_list:
            column_names.append(f"{amino}_{i}")
    return column_names

# Main Function to process the Fasta file/ Plain seqeunces and save the DataFrame
def ohef_fasta(fasta_file):
    headers, sequences = read_fasta(fasta_file)
    padded_sequences = padding_sequence(sequences)
    
    bpf_sequences = [one_hot_encoding_profile(seq) for seq in padded_sequences]
    column_names = feature_column_names()
    ohe_df = pd.DataFrame(bpf_sequences, index=headers, columns=column_names)
    ohe_df.to_csv(os.path.join("feature_used","One_hot_encoding.csv"),index_label="Peptide")
    return ohe_df

#fasta_file = "AOPP.test.2023.fasta"
#output_csv = "binary_independent.csv"
ohe_df = ohef_fasta(fasta_loc)

#====================#====================#=================#

# Extract feature PAAC

def run_pfeature_paac(fasta_file, output_file):
    cmd = ["python", "pfeature_comp.py", 
           "-i", fasta_file,
           "-o", output_file,
           "-j", "PAAC"]
    subprocess.run(cmd, check =False)

def paac_pfeature(fasta_file):
    temp_output = os.path.join("feature_used","PAAC.csv")
    
    run_pfeature_paac(fasta_file, temp_output)
    
    if not os.path.exists(temp_output):
        raise FileNotFoundError("PAAC output not generated!")
    
    paac_df = pd.read_csv(temp_output)
    
    os.remove(temp_output)
    return paac_df

#fasta_file = "AOPP.test.2023.fasta"
paac_df = paac_pfeature(fasta_loc)
paac_df.index = ohe_df.index

#Combine both one hot encoded sequences + PAAC features dataframes

combined_df = pd.concat([ohe_df, paac_df], axis = 1)
combined_df.to_csv(os.path.join("feature_used","Ohe+paac.csv"), index_label="Peptide")
#reading the combined dataframe csv file
combine_df = pd.read_csv(os.path.join("feature_used", "Ohe+paac.csv"))
#276 features are selected from the combined dataframe.
Selected = pd.read_csv(os.path.join(directory_path, "Feature", "selected_features.csv"))
selected_features = Selected['Feature'].tolist()

combined_selected_df = combine_df[selected_features]
combined_selected_df.to_csv("selected_ohe+paac.csv", index = False)

#Prediction

#Function to load the model
def model(file_path):
    model_c = joblib.load(file_path)
    return model_c

#Loading the model
model_path = f"{directory_path}/model/catboost_model_server.pkl"
load_model = model(model_path)

co_df = pd.read_csv("selected_ohe+paac.csv")

cc_df = co_df.drop(columns=['Peptide'])
feature_array = np.array(cc_df)
prob = load_model.predict_proba(feature_array)

y_pred = prob[:,1]

prob_df = pd.DataFrame({'Sequence_ID': co_df['Peptide'],
                        'Prediction_score': y_pred})

if model_chosen == 1:
    #prob_df.to_csv(output_loc + "/prob_df.csv", index =False)
    classification_df = prob_df.copy()
    classification_df['Class']=classification_df['Prediction_score'].apply(lambda x: 'Antioxidant' if x> threshold else 'Non-Antioxidant')

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{output_loc}/classification_ml_{timestamp}.csv"
    classification_df.to_csv(output_file, index = False)
    print("Classification results saved to", output_file)

else:
   #===========================ML+BLAST================================================#
   db_path = f"{directory_path}/blast_db/antioxidant_Training_db"
   output_path = f"{directory_path}/blast_db/100.txt"
   """ 
   blastp_path = os.path.join(directory_path, "ncbi_blast", "bin", "blastp")

   # Add .exe only on Windows
   if os.name == "nt":
       blastp_path += ".exe"
   """     
   #BLast command
    
   blast_command = [
    "blastp",
    "-query", fasta_loc,
    "-task", "blastp-short",
    "-db", db_path,
    "-out", output_path,
    "-outfmt", "6",
    "-evalue", "100",
    "-max_target_seqs", "1",
    "-max_hsps", "1"]

   result = subprocess.run(blast_command, capture_output=True, text=True)
   print("BLAST STDOUT:", result.stdout)
   print("BLAST STDERR:", result.stderr)
   print("Exit code:", result.returncode)

# ALWAYS define blast_df, even if empty
   if result.returncode != 0 or not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
      print("BLAST failed or no hits - using ML only")
      blast_df = pd.DataFrame(columns=[0,1,2])  # Empty but with correct columns
   else:
       try:
            blast_df = pd.read_csv(output_path, sep="\t", header=None)
       except pd.errors.EmptyDataError:
        print("BLAST output empty")
        blast_df = pd.DataFrame(columns=[0,1,2])

# NOW safe to process

   blast_df = blast_df[[0,1,2]]
   blast_df.columns = ["Query_ID", "Database_ID", "Percent_identity"]
    
   blast_df["adjustment"] = 0.0
   
   blast_df.loc[
                blast_df["Database_ID"].str.contains("pos", na =False),
                "adjustment"
                ] = 0.5
   
   blast_df.loc[
                blast_df["Database_ID"].str.contains("neg", na =False),
                "adjustment"
                ] = -0.5
   
   prob_blast_df = prob_df.copy()
   

   prob_blast_df = prob_df.merge(
    blast_df[['Query_ID','adjustment']],
    left_on='Sequence_ID',
    right_on='Query_ID',
    how='left')
   prob_blast_df["adjustment"] = prob_blast_df["adjustment"].fillna(0)
   prob_blast_df["Hybrid_score"] = (prob_blast_df["Prediction_score"] + prob_blast_df["adjustment"]).clip(0, 1)

   classification_hybrid_df = prob_blast_df.copy()
   classification_hybrid_df["Class"] = classification_hybrid_df["Hybrid_score"].apply(lambda x: "Antioxidant" if
                                                                                       x > threshold else "Non-Antioxidant")

   classification_hybrid_df = classification_hybrid_df.iloc[:, [0,1,3,4,5]]
   timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
   output_file = f"{output_loc}/classification_hybrid{timestamp}.csv"
   classification_hybrid_df.to_csv(output_file, index = False)
   print("Classification results saved to", output_file)
