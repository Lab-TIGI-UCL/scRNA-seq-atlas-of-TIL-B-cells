########### Python script for scib. ###########
########### Created by Danwen Qian on 24-04-2024. ###########
########### Last modified by Danwen Qian on 24-04-2024. ###########

import numpy as np
import pandas as pd
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import matplotlib.pyplot as plt
import scib


adata = sc.read_h5ad("~/BCELL_integrated_SCT.h5ad")


adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
adata.obsm["CCA"] = adata.obsm["X_integrated.cca"]
adata.obsm["rPCA"] = adata.obsm["X_integrated.rpca"]
adata.obsm["Harmony"] = adata.obsm["X_harmony"]


 
bm = Benchmarker(
    adata,
    batch_key="study",
    label_key="cluster_anno",
    embedding_obsm_keys=["Unintegrated", "CCA","rPCA","Harmony"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    n_jobs=1,
)



bm.prepare()
bm.benchmark()

bm.plot_results_table(min_max_scale=False,save_dir="~/BCELL")


df = bm.get_results(min_max_scale=False)
df.to_csv("~/metrics_SCT.csv")



