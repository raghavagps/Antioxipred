import sys
import os
import numpy as np
import pandas as pd
import argparse
from itertools import repeat
import re
import gc
import time
import math
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
os.makedirs(output_loc, exist_ok=True)
FEATURE_DIR = os.path.join(os.getcwd(), "feature_used")
os.makedirs(FEATURE_DIR, exist_ok=True)

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
    feature_dir = FEATURE_DIR
    ohe_df.to_csv(os.path.join(feature_dir, "One_hot_encoding.csv"),index_label="Peptide")
    return ohe_df

#fasta_file = "AOPP.test.2023.fasta"
#output_csv = "binary_independent.csv"
ohe_df = ohef_fasta(fasta_loc)

#====================#====================#=================#

# Extract feature PAAC
std = list("ACDEFGHIKLMNPQRSTVWY")
def split(file,v):
    filename, file_extension = os.path.splitext(os.path.basename(file))
    headers, sequences = read_fasta(file)
    df2 = pd.DataFrame([seq.upper() for seq in sequences])
    for i in range(0,len(df2)):
        ss = len(df2[0][i])
        if ss == 0:
            print('\nSequence number',i+1,'has length of',ss,'. Hence, the number of splits should be between 2 to',ss,'. Kindly provide number of splits in the suggested range.')
            sys.exit()
    k1 = []
    for e in range(0,len(df2)):
        s = 0
        k2 = []
        r = 0
        if len(df2[0][e])%v == 0:
            k2.extend(repeat(int(len(df2[0][e])/v),v))
        else:
            r = int(len(df2[0][e])%v)
            k2.extend(repeat(int(len(df2[0][e])/v),v-1))
            k2.append((int(len(df2[0][e])/v))+r)
        for j in k2:
            df3 = df2[0][e][s:j+s]
            k1.append(df3)
            s = j+s
    df4 = pd.DataFrame(k1)
    #df4.to_csv(filename+".split", index = None, header = False, encoding = 'utf-8')
    return df4
def aac_comp(file,out):
    filename, file_extension = os.path.splitext(os.path.basename(file))
    f = open(out, 'w')
    sys.stdout = f
    headers, sequences = read_fasta(file)
    zz = [seq.upper() for seq in sequences]
    print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()
	
def aac_split(file,v,out):
    file1 = split(file,v)
    file2 = file1.iloc[:,0]
    file2.to_csv(os.path.join(FEATURE_DIR, 'sam_input.csv'), index=None, header=None)
    aac_comp(fasta_loc, os.path.join(FEATURE_DIR, 'tempfile_out'))
    df4_1 = pd.read_csv(os.path.join(FEATURE_DIR, 'tempfile_out'))
    df4 = df4_1.iloc[:,:-1]
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),(v*20)):
        pp.append(ss[i:i+(v*20)])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove(os.path.join(FEATURE_DIR, 'sam_input.csv'))                             # FIXED
    os.remove(os.path.join(FEATURE_DIR, 'tempfile_out'))



def val(AA_1, AA_2, aa, mat):
    return sum([(mat[i][aa[AA_1]] - mat[i][aa[AA_2]]) ** 2 for i in range(len(mat))]) / len(mat)
def paac_1(file,lambdaval,w=0.05):
    data1 = pd.read_csv(os.path.join(directory_path, "Data", "data"), sep = "\t") 
    filename, file_extension = os.path.splitext(os.path.basename(file)) 
    headers, sequences = read_fasta(file)  
    df1 = pd.DataFrame([seq.upper() for seq in sequences])
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std)):
        aa[std[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        head.append('_lam' + str(n))
    head = ['PAAC'+str(lambdaval)+sam for sam in head]
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df1)):
        cc = []
        pseudo1 = [] 
        for n in range(1,lambdaval+1):
            cc.append(sum([val(df1[0][k][p], df1[0][k][p + n], aa, dd) for p in range(len(df1[0][k]) - n)]) / (len(df1[0][k]) - n))
            qq = pd.DataFrame(cc)
        pseudo = pseudo1 + [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    ii.to_csv(os.path.join(FEATURE_DIR, f"{filename}.lam"), index=None)
		
def paac(file,lambdaval, out, w=0.05):
    filename, file_extension = os.path.splitext(os.path.basename(file))
    paac_1(file,lambdaval,w=0.05)
    aac_comp(file,os.path.join(FEATURE_DIR, f"{filename}.aac"))
    data1 = pd.read_csv(os.path.join(FEATURE_DIR, f"{filename}.aac"))  
    data2 = pd.read_csv(os.path.join(FEATURE_DIR, f"{filename}.lam")) 
    header = ['PAAC'+str(lambdaval)+'_A','PAAC'+str(lambdaval)+'_C','PAAC'+str(lambdaval)+'_D','PAAC'+str(lambdaval)+'_E','PAAC'+str(lambdaval)+'_F','PAAC'+str(lambdaval)+'_G','PAAC'+str(lambdaval)+'_H','PAAC'+str(lambdaval)+'_I','PAAC'+str(lambdaval)+'_K','PAAC'+str(lambdaval)+'_L','PAAC'+str(lambdaval)+'_M','PAAC'+str(lambdaval)+'_N','PAAC'+str(lambdaval)+'_P','PAAC'+str(lambdaval)+'_Q','PAAC'+str(lambdaval)+'_R','PAAC'+str(lambdaval)+'_S','PAAC'+str(lambdaval)+'_T','PAAC'+str(lambdaval)+'_V','PAAC'+str(lambdaval)+'_W','PAAC'+str(lambdaval)+'_Y','Un']	
    data1.columns = header    
    data3 = pd.concat([data1.iloc[:,:-1],data2], axis = 1).reset_index(drop=True)
    data3.to_csv(out,index=None)
    # CRITICAL FIX
    del data1
    del data2
    del data3

    import gc
    gc.collect()

    import time
    time.sleep(0.5)

    # Safe delete
    for f in [os.path.join(FEATURE_DIR, filename + ".lam"), 
            os.path.join(FEATURE_DIR, filename + ".aac")]:
        try:
            if os.path.exists(f):
                 os.remove(f)
        except PermissionError:
            print(f"Warning: Could not delete {f}")

def paac_split(file,v,lg,out,w):
    paac(file,lg,os.path.join(FEATURE_DIR, 'tempfile_out'),w)
    df4 = pd.read_csv(os.path.join(FEATURE_DIR, 'tempfile_out'))
    head = []
    for j in range(1,v+1):
        for i in df4.columns:
            head.append(i+'_s'+str(j))
    ss = []
    for i in range(0,len(df4)):
        ss.extend(df4.loc[i])
    pp = []
    for i in range(0,len(ss),int(v*(len(ss)/len(df4)))):
        pp.append(ss[i:i+int((v*(len(ss)/len(df4))))])
    df5= pd.DataFrame(pp)
    df5.columns = head
    df5 = round(df5,2)
    df5.to_csv(out,index=None)
    os.remove(os.path.join(FEATURE_DIR, 'tempfile_out'))
"""
def run_pfeature_paac(fasta_file, output_file):
    cmd = ["python", "pfeature_comp.py", 
           "-i", fasta_file,
           "-o", output_file,
           "-j", "PAAC"]
    subprocess.run(cmd, check =False)
"""
def paac_pfeature(fasta_file):
    feature_dir = FEATURE_DIR
    temp_output = os.path.join(feature_dir,"PAAC.csv")
    
    paac(fasta_file,lambdaval=1, out=temp_output, w=0.05)
    
    if not os.path.exists(temp_output):
        raise FileNotFoundError("PAAC output not generated!")
    
    paac_df = pd.read_csv(temp_output)
    paac_df.to_csv(os.path.join(FEATURE_DIR, "PAAC.csv"), index=False)
        
    return paac_df

#fasta_file = "AOPP.test.2023.fasta"
paac_df = paac_pfeature(fasta_loc)
paac_df.index = ohe_df.index

#Combine both one hot encoded sequences + PAAC features dataframes

combined_df = pd.concat([ohe_df, paac_df], axis = 1)
combined_path =os.path.join(FEATURE_DIR, "ohe+paac.csv")

#reading the combined dataframe csv file
combined_df.to_csv(combined_path, index_label="Peptide")
combine_df = pd.read_csv(combined_path)

#276 features are selected from the combined dataframe.
Selected = pd.read_csv(os.path.join(directory_path, "Feature", "selected_features.csv"))
selected_features = Selected['Feature'].tolist()

combined_selected_df = combine_df[selected_features]
selected_path = os.path.join(FEATURE_DIR, "selected_ohe+paac.csv")
combined_selected_df.to_csv(selected_path, index = False)

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

