# AntiOxiPred

A method for predicting antioxidant potential in peptides

## Introduction

AntioxiPred is a computational framework designed to predict the antioxidant potential of peptide sequences with high accuracy using machine learning approaches. The method integrates a similarity-based search using the Basic Local Alignment Search Tool (BLAST) with a CatBoost(Categorical Boosting) classifier. The predictive model is trained using composition-based features, including Pseudo Amino Acid Composition (PAAC) and one-hot encoded sequence profiles, enabling the capture of both global compositional characteristics and positional sequence information.

## Requirements
scikit-learn=1.6.1

Pandas

Numpy

Joblib


You can set up the environment using either requirements.txt (for pip users) or environment.yml (for Conda users).

### Using requirements.txt
```
pip install -r requirements.txt
```

### Using environment.yml
```
conda env create -f env.yml
```
No additional package/tool is required for model = 1 (default model)

### For the hybrid prediction mode (Model 2) in AntioxiPred requires the NCBI BLAST+ software.

BLAST executables are platform-specific. Therefore, users must download the appropriate version of NCBI BLAST+ from the official NCBI website:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

## Minimum USAGE
To know about the available option for the standlone, type the following command:

```
antioxipred -h

```
To run the example, type the following command:

```
antioxipred -f AOPP.test.2023.fasta -o output

```
Here, -f argument is to enter the input file in Fasta format and -o argument is for giving the path to the output directory. By default, the package uses model (-m) = 1 which employs only ML algorithm (Categorical Boosting) to classify the peptide sequences, which generates a prediction file "classification_ml_(datetime).csv" in the specified output directory. If model (-m) = 2 is selected, then the hybrid model is employed (ML + BLAST) to classify the peptide sequences, which generates a prediction file "classification_hybrid(datetime).csv" in the specified output directory.

## Full Usage
```
usage: antioxipred [-h] --file FILE --output OUTPUT [--model MODEL] [--threshold THRESHOLD]
Please provide following arguments for successful run
required arguments:
  --file FILE, -f FILE                   Path to fasta file
  --output OUTPUT, -o OUTPUT             Path to output

optional arguments:

  --model MODEL, -m MODEL                Model selection: 1 for ML only, 2 for ML + BLAST (By default model = 1)
  --threshold THRESHOLD, -t THRESHOLD    Threshold for classification (can be any value between 0-1 for model = 1 (by default = 0.5) and 0-2 for model = 2 (by default = 0.52))

For help:
  -h, --help            show this help message and exit

```


## Standalone Minimum Usage

Run the program using:

```
python3 antioxipred.py -f AOPP.test.2023.fasta -o result
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

