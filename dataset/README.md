# 📂 Dataset Information

This folder contains all input data required to run the **Multi-Modal Biomedical Data Integration Pipeline**.  
The data are derived from the original raw dataset hosted on **NCBI GEO** and include preprocessed CSV and H5AD files.

---

## 🧬 Raw Data Source

- **GEO Accession:** [GSE300475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)  
- This dataset serves as the reference for the **original raw biological data** used in this project.

---

## 📁 Directory Structure

### 1️⃣ CSV Files  
All `.csv` files should be placed in dataset/csv/

These files contain preprocessed tabular data for model training and metadata reference.

| File | Description |
|------|--------------|
| `hvg_genes.csv` | List of highly variable genes used for feature selection |
| `transcriptomic_features.csv` | Transcriptomic feature matrix (samples × genes) |
| `fluxomic_features.csv` | Flux-based feature matrix (samples × reactions) |
| `metadata.csv` | Sample-level metadata with response labels (Responder / Non-responder) |

---

### 2️⃣ H5AD Files  
All `.h5ad` files should be placed in dataset/h5ad/

These contain the **AnnData objects** 


---

## ⚠️ Notes
- Ensure **all required `.csv` and `.h5ad` files** are downloaded and correctly placed before running the analysis scripts locally.  
- The `metadata.csv` file must have a `response` column containing binary class labels (`Responder` / `Non-responder`).  
- Sample indices across all files must be **consistent** for successful merging.

---

🧠 *This folder acts as the local data hub for all processed and intermediate datasets used throughout the project.*

