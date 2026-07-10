def save_processed_adata(adata, output_path):
    adata.write_h5ad(output_path)

def export_obs_metadata(adata, output_path):
    adata.obs.to_csv(output_path)

def export_cluster_counts(adata, column, output_path):
    adata.obs[column].value_counts().to_csv(output_path)
