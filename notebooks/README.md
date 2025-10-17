# 📓 Notebooks

This directory contains all Jupyter notebooks used in the **data processing**, **feature generation**, and **machine-learning evaluation** workflow.  
Each notebook is **numbered sequentially** to indicate the logical order of execution.

---

## 🧭 Recommended Execution Order
1. `0_data_acquisition.ipynb` → Load raw data and prepare metadata.  
2. `1_preprocessing_genes.ipynb` → Perform QC, normalization, and differential expression analysis.  
3. `2_generate_flux.ipynb` → Compute fluxomic features using `scFBApy`.  
4. `3_experiment_1.ipynb` → ML on transcriptomic features.  
5. `3_experiment_2.ipynb` → ML on fluxomic features.  
6. `3_experiment_3.ipynb` → ML on multimodal (combined) features.

---

### **0_data_acquisition.ipynb**
- Loads raw **scRNA-seq** data from GEO.  
- Prepares cell-level metadata, patient identifiers, and response labels.  
- Exports structured data for downstream analysis.

---

### **1_preprocessing_genes.ipynb**
- Performs **quality control (QC)** by filtering low-quality cells and removing outliers.  
- Normalizes expression counts and conducts **Differential Expression (DE)** analysis between responder and non-responder groups.  
- Generates QC visualizations (e.g., violin plots, expression distributions).  
- Saves the cleaned and preprocessed data for downstream flux generation.

---

### **2_generate_flux.ipynb**
- Generates **cell-level metabolic fluxes** from scRNA-seq data using `scFBApy`.  
- Produces structured flux matrices compatible with machine-learning models.  
- Outputs are saved for multimodal integration.

---

### **3_experiment_1.ipynb — Transcriptomic Features Only**
- Uses **gene expression (HVG)** features for machine learning modeling.  
- Implements models: Logistic Regression, Random Forest, XGBoost, and ANN.  
- Performs hyperparameter tuning using **Optuna**.  
- Includes **SHAP contribution plots** for model interpretability.  
- Evaluates predictive performance based solely on transcriptomic data.

---

### **3_experiment_2.ipynb — Fluxomic Features Only**
- Uses **fluxomic (metabolic reaction)** features for modeling.  
- Applies the same ML pipeline and tuning process as in Experiment 1.  
- Evaluates predictive power of flux-based data for treatment response.  
- Provides **feature importance** and **SHAP**-based interpretability visualizations.

---

### **3_experiment_3.ipynb — Multimodal Integration (Transcriptomics + Fluxomics)**
- Integrates **transcriptomic** and **fluxomic** features into a unified dataset  
  (based on top highly variable genes and corresponding fluxomic features).  
- Trains and tunes **multimodal models** using Logistic Regression, RF, XGBoost, and ANN.  
- Performs **Optuna** hyperparameter optimization and **SHAP-based** interpretation.  
- Visualizes **contribution plots** comparing top features across both modalities.

---

📘 *Each notebook is self-contained and can be executed independently for modality-specific or multimodal analysis.*