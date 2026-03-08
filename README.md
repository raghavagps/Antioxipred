ANTIOXIPRED
A method for predicting antioxidant potential in peptides

Introduction
AntioxiPred is a computational framework designed to predict the antioxidant potential of peptide sequences with high accuracy using machine learning approaches. The method integrates a similarity-based search using the Basic Local Alignment Search Tool (BLAST) with a CatBoost(Categorical Boosting) classifier. The predictive model is trained using composition-based features, including Pseudo Amino Acid Composition (PAAC) and one-hot encoded sequence profiles, enabling the capture of both global compositional characteristics and positional sequence information.

Requirements
scikit-learn=1.6.1
Pandas
Numpy
Joblib


You can set up the environment using either requirements.txt (for pip users) or environment.yml (for Conda users).

Using requirements.txt
pip install -r requirements.txt

Using environment.yml
conda env create -f env.yml

No additional package/tool is required for model = 1 (default model)

For the hybrid prediction mode (Model 2) in AntioxiPred requires the NCBI BLAST+ software.

A Windows-compatible BLAST executable is already included in this repository under:

ncbi_blast/bin/

Since the tool was developed and tested on Windows, Windows users can directly use the provided BLAST executable.

For Linux or macOS users

BLAST executables are platform-specific. Therefore, Linux and macOS users must download the appropriate version of NCBI BLAST+ from the official NCBI website:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Steps:

Download the BLAST+ package for your operating system.

Extract the downloaded archive.

Copy the bin folder from the extracted BLAST package.

Replace the existing folder inside the AntioxiPred directory:

ANTIOXIPRED/
│
├── antioxipred.py
├── ncbi_blast/
│   └── bin/   ← replace this folder with the one from the downloaded BLAST package

After replacing the bin folder, the program will automatically use the correct BLAST executable.



Standalone minimum usage
python3 antioxipred.py -f AOPP.test.2023.fasta -o output
Arguments description
Input File: It allow users to provide input in FASTA format.

Output File: Program will save the results to this folder

Model: User can pick which model to run, model = 1 runs only ML model (CATBOOT classifier), whereas model = 2 runs hybrid model (ML + BLAST), by default the tool runs model = 1

Threshold: User can provide threshold for classification (can be any value between 0-1 for model = 1 (by default = 0.50) and 0-2 for model = 2 (by default = 0.50))

ANTIOXIPRED Package Files
It contantain following files, brief description of these files given below

INSTALLATION : Installations instructions

README.md : This file provide information about this package

catboost_model_server.pkl : This file contains the pickled version of model

antioxipred.py : Main python program

AOPP.test.2023.fasta : Example file contain peptide sequences in FASTA format

blast_db : Database for BLAST search

