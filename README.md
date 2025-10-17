
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
│ ├── 3_experiment_1.ipynb
│ ├── 3_experiment_2.ipynb
│ └── 3_experiment_3.ipynb
│
|
│── models/model.xml (scFBA XML model) 
|
│── scripts (scFBA python scripts)
│ ├── pipeline_scFBApy.py
│ └── utils_scFBApy.py
│
└── dataset/
│ ├── csv
│ └── h5ad

```

**Notes:**
1. `notebooks/` contains Jupyter notebooks with descriptive, numbered filenames that reflect the logical order of the analysis workflow.
2. `dataset/` contains all datasets used or generated throughout experiments for full reproducibility.
3. `models/` contains cobra model of scFBApy.
4. `scripts/` contains scFBApy scripts.
5. `requirements.txt` lists all Python dependencies for reproducibility.

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

### 3_experiment_1.ipynb — Transcriptomic Features Only
- Uses **gene expression features** (HVGs) for ML modeling.  
- Models implemented: Logistic Regression, Random Forest, XGBoost, and ANN.  
- Hyperparameter tuning performed using **Optuna**.  
- Includes **SHAP contribution plots** for model interpretability.  
- Evaluates predictive performance based on transcriptomic data only.

### 3_experiment_2.ipynb — Fluxomic Features Only
- Uses **fluxomic features** (metabolic reaction–based features).  
- Applies the same modeling and tuning pipeline as in Experiment 1.  
- Evaluates how well flux-based data alone can predict treatment response.  
- Includes feature importance visualization and SHAP explanations.

### 3_experiment_3.ipynb — Multimodal Integration (Transcriptomics and Fluxomics)
- Integrates **transcriptomic** and **fluxomic** features into a unified dataset  
  (based on top highly variable genes and corresponding fluxomic features).  
- Implements the same ML pipeline and models (Logistic Regression, RF, XGBoost, ANN).  
- Performs **Optuna** hyperparameter optimization and **SHAP** interpretation in multimodal mode.  
- Visualizes contribution plots comparing top features across both modalities.  

📘 *Each notebook is self-contained and can be executed independently for modality-specific or multimodal analysis.*

---

## Data Source

The primary dataset comes from the **NCBI Gene Expression Omnibus (GEO):**  

🔗 [GSE300475 – Single-cell RNA-seq of HR+ breast cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)  

- Includes single-cell transcriptomes and TCR repertoire data across multiple patients.   
- Patient response status (Responder vs Non-responder) was mapped using the associated publication:  
  🔗 [npj | Breast Cancer – Dynamic single-cell systemic immune responses in immunotherapy-treated early-stage HR+ breast cancer patients](https://www.nature.com/articles/s41523-025-00776-1).  
- Processed into feature matrices for downstream analysis in this project.

---

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
- **`00_data_acquisition.ipynb`** → Download and organize the raw datasets.  
- **`01_preprocessing_genes.ipynb`** → Perform QC, outlier removal, normalization, and differential expression (DE) analysis on transcriptomic data.  
- **`02_generate_flux.ipynb`** → Infer single-cell metabolic fluxes using *scFBApy*.  
- **`3_experiment_1.ipynb`** → Train, tune (*Optuna*), and evaluate machine-learning models using **transcriptomic (gene expression)** features only, including **SHAP contribution plots** for top genes.
- **`3_experiment_2.ipynb`** → Train, tune (*Optuna*), and evaluate models using **fluxomic (metabolic reaction)** features only, with **feature importance** and **SHAP interpretability** visualizations.
- **`3_experiment_3.ipynb`** → Train, tune (*Optuna*), and evaluate **multimodal models** integrating **transcriptomic + fluxomic** features, and visualize **SHAP contribution plots** comparing top features across both modalities.
---
