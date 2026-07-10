def validate_required_fields(adata):
    checks = {
        "has_X": adata.X is not None,
        "has_obs": adata.obs is not None,
        "has_var": adata.var is not None,
        "has_pca": "X_pca" in adata.obsm,
        "has_umap": "X_umap" in adata.obsm,
        "has_neighbors": "neighbors" in adata.uns,
    }
    return checks
