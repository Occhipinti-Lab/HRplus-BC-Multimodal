
# Multiscale Flux Balance Analysis: Bridging Multi-Omics and Single-Cell Data for Therapy Response Prediction in Hormone Receptor-Positive Breast Cancer

This repository contains code and experiments for predicting treatment response in HR+ breast cancer using **single-cell RNA sequencing (scRNA-seq) data**, **fluxomics features inferred with scFBApy**, and a **multimodal integration of both**. The project evaluates predictive performance when models are trained on **transcriptomics**, **fluxomics**, or on their **combined multimodal features**, in order to assess whether integration provides improved accuracy over single-modality approaches.

---

## About

This project provides a reproducible pipeline for prediction of breast cancer treatment response at the single-cell level. It includes data acquisition, preprocessing, feature extraction, machine learning modeling, and explainability analyses.

**Key goals:**
- Evaluate predictive performance of transcriptomic and fluxomic features individually.
- Investigate the added value of multimodal integration.
- Provide a fully documented workflow for reproducibility and future extension.

---

## Repository Structure
```
project/
│── README.md
│── requirements.txt
│── .gitignore
│
│── notebooks/ 
│ ├── 0_data_acquistion.ipynb
│ ├── 1_preprocessing_genes.ipynb
│ ├── 2_generate_flux.ipynb
│ └── 3_ml_modality_evaluator.ipynb
│
|
│── models/model.xml (scFBA XML model) 
|
│── scripts
│ ├── pipeline_scFBApy.py
│ └── utils_scFBApy.py
│
└── dataset/
│ ├── csv
│ └── h5ad

```

**Notes:**
1. Notebook names are descriptive and numbered to indicate workflow order.  
2. `results/` stores experiment outputs, including csv for all modality and models.
3. `dataset/` contains all datasets used or generated throughout experiments for full reproducibility.
4. `models/` contains cobra model of scFBApy.
5. `scripts/` contains scFBApy scripts.
6. `requirements.txt` lists all Python dependencies for reproducibility.

---

## Notebooks

### 0_data_acquistion.ipynb
- Loads raw scRNA-seq data from GEO.  
- Prepares cell-level metadata, patient IDs, and response labels.  

### 1_preprocessing_genes.ipynb
- Performs quality control (QC): filters low-quality cells, removes outliers, and normalizes expression counts.
- Conducts Differential Expression (DE) analysis between responder and non-responder groups.
- Generates visualization plots including violin plots and gene expression distributions for QC validation.
- Saves the cleaned and preprocessed data for downstream flux generation.

### 2_generate_flux.ipynb
- Generates cell-level fluxes from scRNA-seq using `scFBApy`.    
- Outputs structured flux matrices for machine learning.

### 3_ml_modality_evaluator.ipynb 
- Uses gene expression features for ML modeling.  
- Uses fluxomic features for ML modeling.
- Combines transcriptomics and fluxomics into a unified dataset (based on highly variance gene and fluxomic features) 
- Models: Logistict Regression, Random Forest, XGBoost, ANN.
- Optuna for hyperparameter tuning
- SHAP contribution plots for interpretability in multimodal mode
- Compares predictive performance of:
  1. Transcriptomics
  2. Fluxomics
  3. Multimodal concatenated features

---

## Data Source

The primary dataset comes from the **NCBI Gene Expression Omnibus (GEO):**  

🔗 [GSE300475 – Single-cell RNA-seq of HR+ breast cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)  

- Includes single-cell transcriptomes and TCR repertoire data across multiple patients.   
- Patient response status (Responder vs Non-responder) was mapped using the associated publication:  
  🔗 [npj | Breast Cancer – Dynamic single-cell systemic immune responses in immunotherapy-treated early-stage HR+ breast cancer patients](https://www.nature.com/articles/s41523-025-00776-1).  
- Processed into feature matrices for downstream analysis in this project.


## How to Use

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Occhipinti-Lab/HRplus-BC-Multimodal.git
   cd HRplus-BC-Multimodal
---
2. **Install Dependencies**
    ```bash
    pip install -r requirements.txt
---
3. **Run notebooks in order:**
Execute the notebooks in the following order to reproduce the full pipeline:
1. ***`00_data_acquisition.ipynb`** → Download and organize the raw datasets.  
2. ***`01_preprocessing_genes.ipynb`** → Perform QC, outlier removal, normalization, and differential expression (DE) analysis on transcriptomic data.  
3. ***`02_generate_flux.ipynb`** → Infer single-cell metabolic fluxes using *scFBApy*.  
4. ***`03_ml_modality_evaluator.ipynb`** → Train, tune (Optuna), and evaluate machine-learning models across transcriptomic, fluxomic, and multimodal datasets.
---
