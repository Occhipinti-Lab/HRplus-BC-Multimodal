# Running the experiments on Google Colab

This folder contains Jupyter notebooks prepared for **Google Colaboratory**, enabling users to run the experiments directly in the cloud without local setup. The notebooks replicate the three main experiments from the primary repository, covering **transcriptomic**, **fluxomic**, and **multimodal** model training and evaluation.

## Folder Contents

| Notebook | Description |
|----------|-------------|
| `Colab_3_experiment_1.ipynb` | Experiment using **transcriptomic (gene expression)** features only.|
| `Colab_3_experiment_2.ipynb` | Experiment using **fluxomic (metabolic flux)** features only|
| `Colab_3_experiment_3.ipynb` | Experiment using **combined transcriptomic + fluxomic** features. Demonstrates multimodal model training, evaluation, and SHAP contribution comparison across modalities. |

---

## Usage on Google Colaboratory

1. Open the desired notebook in **Google Colab**:
   Click `File → Open notebook → Upload` or drag the notebook file into Colab.  
   
2. **Install dependencies**:
   Each notebook contains a dedicated cell to install required packages in the Colab environment. Simply run this cell to set up the environment.

3. **Mount Google Drive (optional)**:
   If you plan to read/write data or save model outputs, mount your Google Drive:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')

