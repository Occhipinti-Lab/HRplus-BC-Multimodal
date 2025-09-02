
# Multimodal Prediction of Breast Cancer Treatment Response  
This repository contains code and experiments for predicting treatment response in HR+ breast cancer patients using **single-cell RNA sequencing (scRNA-seq)**, **fluxomics (scFBApy)**, and a **multimodal integration of both**. The project explores how transcriptomics and metabolic flux modeling can be combined to improve predictive performance.

## Repository Structure
The repository is organized into three main experimental workflows:

1. **Experiment 1 – Transcriptomics**  
   - Uses **scRNA-seq gene expression data** as input features.  
   - Preprocesses raw data, performs dimensionality reduction (PCA/UMAP), and trains predictive models such as **Random Forest**, **XGBoost**, and **ANNs**.  
   - Includes **explainability analysis (SHAP)** to identify top genes contributing to treatment response.  
   - Code: `experiment1_transcriptomics.ipynb`

2. **Experiment 2 – Fluxomics**  
   - Applies **scFBApy** to infer **cell-level metabolic fluxes** from single-cell gene expression data.  
   - Fixes the metabolic model to include boundary/exchange reactions and processes flux outputs into structured feature matrices.  
   - Predictive models are trained on fluxomic features to evaluate their ability to capture metabolic signatures of response.  
   - Code: `experiment2_fluxomics.ipynb`

3. **Experiment 3 – Multimodal Integration**  
   - Combines **transcriptomics** and **fluxomics** features into a **unified dataset**.  
   - Trains machine learning models on:  
     1. Transcriptomics only  
     2. Fluxomics only  
     3. Concatenated multimodal features  
   - Compares performance across modalities and evaluates whether multimodal integration improves predictive accuracy.  
   - Code: `experiment3_multimodal.ipynb`

## Data Source
The dataset used in this project comes from the **NCBI Gene Expression Omnibus (GEO)**:  

🔗 [GSE300475 – Single-cell RNA-seq of HR+ breast cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)

- Contains **single-cell transcriptomes and TCR repertoire data** across multiple patients.  
- Patients are annotated with **treatment response status** (Responder vs Non-responder).  
- These data are processed into feature matrices for downstream experiments.

---

## 🚀 How to Use
1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/your-repo-name.git
   cd your-repo-name
