#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel Fluxomics Pipeline using scFBApy
Author: Sadegh
"""

import sys
import os
import glob
import time
import numpy as np
import pandas as pd
import warnings
from tqdm.auto import tqdm
from anndata import AnnData
from pathlib import Path
from multiprocessing import Pool
import cobra
from scripts.utils_scFBApy import scFBApy, repairNeg

# ============================================================
# 🔸 Redirect all stdout to both terminal and log file
# ============================================================
class Tee:
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "a", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()  # write immediately to file

    def flush(self):
        pass


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))


# set log file path
base_dir = Path().resolve()
timestamp = time.strftime("%Y%m%d_%H%M%S")
log_file = base_dir / f"flux_pipeline_{timestamp}.log"
sys.stdout = Tee(log_file)

# ============================================================
# 🔸 Suppress unwanted warnings globally
# ============================================================
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message="ChainedAssignmentError")

# -------------------------
# Paths
# -------------------------
input_file = base_dir / "dataset/transcriptomic_hvg.csv"
meta_data = base_dir / "dataset/metadata.csv"
model_file = base_dir / "models/model.xml"
output_dir = base_dir / "flux_batch"
output_dir.mkdir(exist_ok=True)
final_output = base_dir / "dataset/fluxomics.csv"

# -------------------------
# Load transcriptomics data
# -------------------------
print("\n📥 Loading transcriptomics data...")
trans = pd.read_csv(input_file, index_col=0)
print(f"Expression matrix shape: {trans.shape}")

# -------------------------
# Load COBRA model
# -------------------------
print("\n🧬 Loading COBRA model...")
model = cobra.io.read_sbml_model(model_file)
model_genes = [g.id for g in model.genes]

# -------------------------
# Filter expression to model genes
# -------------------------
overlap_genes = [g for g in trans.columns if g in model_genes]
expr_filtered = trans[overlap_genes]
adata_expr = AnnData(expr_filtered)
print(f"Number of overlapping genes: {len(overlap_genes)}")
print(f"AnnData shape: {adata_expr.shape}")

# -------------------------
# Parameters
# -------------------------
batch_size = 100
objective = "Biomass"
n_cells = adata_expr.shape[0]
n_batches = (n_cells + batch_size - 1) // batch_size
print(f"\n⚙️ Total cells: {n_cells} | Batch size: {batch_size} → {n_batches} batches")

# -------------------------
# Function for one batch
# -------------------------
def run_batch(batch_idx):
    import pandas as pd
    import numpy as np
    import warnings
    from scripts.utils_scFBApy import scFBApy

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", message="ChainedAssignmentError")

    start = batch_idx * batch_size
    end = min((batch_idx + 1) * batch_size, n_cells)
    batch_file = output_dir / f"flux_batch_{batch_idx}.csv"

    if batch_file.exists():
        return f"🟢 Skipped batch {batch_idx+1} (already done)"

    adata_batch = adata_expr[start:end, :].copy()

    try:
        adata_flux_batch = scFBApy(
            model_orig=model,
            adata=adata_batch,
            objective=objective,
            cooperation=True,
            compute_fva=True,
            npop_fva=5,
            eps=0.001,             
            type_ras_normalization="max",
            and_expression=np.nanmin,
            or_expression=np.nansum,
            fraction_of_optimum=0,
            processes=1,
            round_c=10
        )

        flux_df = pd.DataFrame(
            adata_flux_batch.X,
            index=adata_flux_batch.obs.index,
            columns=adata_flux_batch.var.index
        )
        flux_df.to_csv(batch_file)
        return f"✅ Finished batch {batch_idx+1}/{n_batches}"

    except Exception as e:
        return f"⚠️ Error in batch {batch_idx}: {e}"

# -------------------------
# Parallel execution
# -------------------------
if __name__ == "__main__":
    print(f"\n🚀 Starting parallel flux computation using 6 processes...")
    start_time = time.time()

    with Pool(processes=15) as pool:
        results = list(tqdm(pool.imap(run_batch, range(n_batches)), total=n_batches))

    for r in results:
        print(r)

    duration = time.time() - start_time
    print(f"\n⏱️ Total runtime: {duration/60:.1f} minutes")

    # -------------------------
    # Combine all flux batches
    # -------------------------
    print("\n📦 Combining all batch files...")
    batch_files = sorted(output_dir.glob("flux_batch_*.csv"))
    combined_flux = pd.concat([pd.read_csv(f, index_col=0) for f in batch_files])
    print(f"Combined flux shape: {combined_flux.shape}")

    # -------------------------
    # Add metadata
    # -------------------------
    meta_df = pd.read_csv(meta_data, index_col=0)
    combined_flux.index = meta_df.index.values[: combined_flux.shape[0]]
    combined_flux["response"] = meta_df["response"].values[: combined_flux.shape[0]]

    combined_flux.to_csv(final_output, index=True)
    print(f"\n💾 Final fluxomics file saved to:\n{final_output}")
    print(f"🪶 Full log saved at: {log_file}")

    # -------------------------
    # Validation step (check duplicate values)
    # -------------------------
    print("\n🔎 Checking for duplicate flux value rows...")

    flux_columns = [c for c in combined_flux.columns if c != "response"]

    dup_mask = combined_flux.duplicated(subset=flux_columns, keep=False)
    dup_count = dup_mask.sum()
    dup_ratio = dup_count / combined_flux.shape[0]

    print(f"⚙️ Total rows: {combined_flux.shape[0]}")
    print(f"⚠️ Duplicate (identical) flux patterns: {dup_count} rows ({dup_ratio:.2%})")

    if dup_count == 0:
        print("✅ No duplicate flux rows — all cells have unique metabolic profiles.")
    elif dup_ratio < 0.05:
        print("✅ Minor duplication detected (<5%) — acceptable biological similarity.")
    elif dup_ratio < 0.2:
        print("⚠️ Moderate duplication (5–20%) — check normalization and rounding parameters.")
    else:
        print("🚨 High duplication (>20%) — model may be collapsing; recheck eps/round_c.")
