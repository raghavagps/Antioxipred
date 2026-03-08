## ANTIOXIPRED

A method for predicting antioxidant potential in peptides

# Introduction

AntioxiPred is a computational framework designed to predict the antioxidant potential of peptide sequences with high accuracy using machine learning approaches. The method integrates a similarity-based search using the Basic Local Alignment Search Tool (BLAST) with a CatBoost(Categorical Boosting) classifier. The predictive model is trained using composition-based features, including Pseudo Amino Acid Composition (PAAC) and one-hot encoded sequence profiles, enabling the capture of both global compositional characteristics and positional sequence information.

# Requirements
scikit-learn=1.6.1
Pandas
Numpy
Joblib


# You can set up the environment using either requirements.txt (for pip users) or environment.yml (for Conda users).

# Using requirements.txt

pip install -r requirements.txt

# Using environment.yml

conda env create -f env.yml

No additional package/tool is required for model = 1 (default model)

# For the hybrid prediction mode (Model 2) in AntioxiPred requires the NCBI BLAST+ software.

For Windows, Linux and macOS users

BLAST executables are platform-specific. Therefore, users must download the appropriate version of NCBI BLAST+ from the official NCBI website:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

## Steps:

1.Download the BLAST+ package for your operating system.
2.Extract the downloaded archive.
3.Copy the bin folder from the extracted BLAST package.
4.Replace the existing folder inside the AntioxiPred directory:

```
ANTIOXIPRED/
│
├── antioxipred.py
├── ncbi_blast/
│   └── bin/← replace this folder with the one from the downloaded BLAST package
``` 

After replacing the bin folder, the program will automatically use the correct BLAST executable.



## Standalone Minimum Usage

Run the program using:

```bash
python3 antioxipred.py -f AOPP.test.2023.fasta -o output
```

---

## Arguments Description

### Input File (`-f`)
Allows users to provide input peptide sequences in **FASTA format**.

### Output File (`-o`)
The program saves the prediction results in the specified **output folder**.

### Model (`-m`)
Users can choose which model to run:

- **model = 1** → Runs only the **Machine Learning model (CatBoost classifier)**
- **model = 2** → Runs the **Hybrid model (Machine Learning + BLAST)**

**Default:** `model = 1`

### Threshold (`-t`)
Users can provide a threshold for classification.

- For **model = 1** → threshold range: `0 – 1` (default = `0.50`)
- For **model = 2** → threshold range: `0 – 2` (default = `0.50`)

---

## ANTIOXIPRED Package Files

The package contains the following files:

| File | Description |
|-----|-------------|
| `INSTALLATION` | Installation instructions |
| `README.md` | Documentation and usage instructions |
| `catboost_model_server.pkl` | Pickled CatBoost prediction model |
| `antioxipred.py` | Main Python script used to run the prediction |
| `AOPP.test.2023.fasta` | Example FASTA file containing peptide sequences |
| `blast_db/` | Database used for BLAST similarity search |
