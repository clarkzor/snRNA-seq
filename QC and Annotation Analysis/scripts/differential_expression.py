def run_rank_genes_groups(adata, groupby, method="t-test"):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
    return adata

def export_rank_genes_groups(adata, output_path):
    ...
