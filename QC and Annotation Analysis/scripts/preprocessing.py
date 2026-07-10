def normalize_log_transform(adata, target_sum=10000):
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=target_sum)
    sc.pp.log1p(adata)
    return adata

def select_highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.05
):
    sc.pp.highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp
    )
    return adata[:, adata.var["highly_variable"]].copy()

def run_pca_neighbors_umap(
    adata,
    n_comps=100,
    n_neighbors=50,
    n_pcs=99,
    random_state=42
):
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, random_state=random_state)
    return adata
