"""Output helpers for the Xenopus gastrula snRNA-seq notebooks.

These functions provide:
- exact, descriptive figure filenames;
- CSV export of Scanpy rank_genes_groups results;
- annotation count and percentage tables;
- cluster-to-annotation cross-tabulations;
- label-refinement transition tables.

The module is written for the package versions recorded by this repository,
including Scanpy 1.9.6 and pandas 2.1.4.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Iterable, Optional

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def slugify(value: str) -> str:
    """Convert a label into a filesystem-safe lowercase slug."""
    value = str(value).strip().lower()
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^a-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("._-")


def save_current_figure(
    output_dir: Path | str,
    filename: str,
    *,
    dpi: int = 300,
    bbox_inches: str = "tight",
    close: bool = True,
) -> Path:
    """Save the active Matplotlib/Scanpy figure with an exact filename."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = filename if filename.lower().endswith(".png") else f"{filename}.png"
    path = output_dir / filename

    fig = plt.gcf()
    fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)

    if close:
        plt.close(fig)

    print(f"Saved figure: {path}")
    return path


def save_dataframe(
    dataframe: pd.DataFrame,
    output_dir: Path | str,
    filename: str,
    *,
    index: bool = False,
) -> Path:
    """Save a DataFrame as CSV and create the output directory if needed."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = filename if filename.lower().endswith(".csv") else f"{filename}.csv"
    path = output_dir / filename
    dataframe.to_csv(path, index=index)

    print(f"Saved table: {path}")
    return path


def _rank_groups(adata, key: str) -> list[str]:
    """Return group names stored in a Scanpy rank_genes_groups result."""
    result = adata.uns[key]
    names = result["names"]

    if isinstance(names, pd.DataFrame):
        return [str(column) for column in names.columns]

    dtype_names = getattr(getattr(names, "dtype", None), "names", None)
    if dtype_names:
        return [str(group) for group in dtype_names]

    params = result.get("params", {})
    groups = params.get("groups")
    if groups is not None:
        return [str(group) for group in groups]

    raise ValueError(f"Could not determine groups stored under adata.uns[{key!r}].")


def rank_genes_groups_to_long_df(adata, key: str) -> pd.DataFrame:
    """Convert all groups in a rank_genes_groups result into one long table."""
    groups = _rank_groups(adata, key)
    frames: list[pd.DataFrame] = []

    pts = adata.uns[key].get("pts")
    pts_rest = adata.uns[key].get("pts_rest")

    for group in groups:
        frame = sc.get.rank_genes_groups_df(
            adata,
            group=group,
            key=key,
        ).copy()

        frame.insert(0, "group", group)

        if pts is not None and group in pts.columns:
            frame["pct_nz_group"] = frame["names"].map(pts[group])

        if pts_rest is not None and group in pts_rest.columns:
            frame["pct_nz_reference"] = frame["names"].map(pts_rest[group])

        frames.append(frame)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def export_rank_genes_groups(
    adata,
    *,
    key: str,
    output_dir: Path | str,
    prefix: str,
    groupby: Optional[str] = None,
    fdr_max: float = 0.05,
    abs_log2fc_min: float = 1.0,
    top_n: int = 25,
) -> dict[str, Path]:
    """Export full, significant, and top-ranked DEG results plus metadata."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = slugify(prefix)

    results = rank_genes_groups_to_long_df(adata, key)
    paths: dict[str, Path] = {}

    paths["all_results"] = save_dataframe(
        results,
        output_dir,
        f"{prefix}_all_results.csv",
    )

    significant = results.copy()
    if "pvals_adj" in significant.columns:
        significant = significant[significant["pvals_adj"] <= fdr_max]
    if "logfoldchanges" in significant.columns:
        significant = significant[
            significant["logfoldchanges"].abs() >= abs_log2fc_min
        ]

    paths["significant_results"] = save_dataframe(
        significant,
        output_dir,
        f"{prefix}_significant_fdr{str(fdr_max).replace('.', 'p')}_"
        f"abslog2fc{str(abs_log2fc_min).replace('.', 'p')}.csv",
    )

    rank_column = "scores" if "scores" in results.columns else results.columns[-1]
    top_results = (
        results.sort_values(["group", rank_column], ascending=[True, False])
        .groupby("group", observed=True, group_keys=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    paths["top_results"] = save_dataframe(
        top_results,
        output_dir,
        f"{prefix}_top{top_n}_per_group.csv",
    )

    params = dict(adata.uns[key].get("params", {}))
    params.update(
        {
            "rank_genes_groups_key": key,
            "export_prefix": prefix,
            "fdr_max": fdr_max,
            "abs_log2fc_min": abs_log2fc_min,
            "top_n": top_n,
        }
    )

    parameters_path = output_dir / f"{prefix}_parameters.json"
    with parameters_path.open("w", encoding="utf-8") as handle:
        json.dump(params, handle, indent=2, default=str)
    paths["parameters"] = parameters_path
    print(f"Saved DEG parameters: {parameters_path}")

    if groupby is None:
        groupby = params.get("groupby")

    if groupby and groupby in adata.obs.columns:
        group_counts = (
            adata.obs[groupby]
            .astype(str)
            .value_counts(dropna=False)
            .rename_axis("group")
            .reset_index(name="cell_count")
        )
        group_counts["percent"] = (
            100 * group_counts["cell_count"] / group_counts["cell_count"].sum()
        )
        paths["group_counts"] = save_dataframe(
            group_counts,
            output_dir,
            f"{prefix}_group_counts.csv",
        )

    return paths


def annotation_summary(
    adata,
    annotation_col: str = "region_annotation",
) -> pd.DataFrame:
    """Return count and percentage for each annotation."""
    counts = (
        adata.obs[annotation_col]
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis(annotation_col)
        .reset_index(name="cell_count")
    )
    counts["percent"] = 100 * counts["cell_count"] / counts["cell_count"].sum()
    return counts


def save_annotation_summary(
    adata,
    *,
    output_dir: Path | str,
    filename: str,
    annotation_col: str = "region_annotation",
) -> Path:
    """Save annotation counts and percentages."""
    return save_dataframe(
        annotation_summary(adata, annotation_col=annotation_col),
        output_dir,
        filename,
    )


def save_cluster_annotation_crosstab(
    adata,
    *,
    cluster_col: str,
    annotation_col: str,
    output_dir: Path | str,
    filename: Optional[str] = None,
    normalize: bool = False,
) -> Path:
    """Save a cluster-by-annotation contingency table."""
    if cluster_col not in adata.obs.columns:
        raise KeyError(f"{cluster_col!r} is not present in adata.obs.")
    if annotation_col not in adata.obs.columns:
        raise KeyError(f"{annotation_col!r} is not present in adata.obs.")

    table = pd.crosstab(
        adata.obs[cluster_col].astype(str),
        adata.obs[annotation_col].astype(str),
        normalize="index" if normalize else False,
    )

    if normalize:
        table = table * 100

    table = table.reset_index()
    suffix = "_row_percent" if normalize else "_counts"
    filename = filename or (
        f"{slugify(cluster_col)}_by_{slugify(annotation_col)}{suffix}.csv"
    )

    return save_dataframe(table, output_dir, filename)


def save_label_transition_table(
    adata,
    *,
    before_col: str,
    after_col: str,
    output_dir: Path | str,
    filename: Optional[str] = None,
) -> Path:
    """Save a before-versus-after label transition table."""
    if before_col not in adata.obs.columns or after_col not in adata.obs.columns:
        missing = [
            column
            for column in [before_col, after_col]
            if column not in adata.obs.columns
        ]
        raise KeyError(f"Missing annotation columns: {missing}")

    table = pd.crosstab(
        adata.obs[before_col].astype(str),
        adata.obs[after_col].astype(str),
    ).reset_index()

    filename = filename or (
        f"{slugify(before_col)}_to_{slugify(after_col)}_transitions.csv"
    )
    return save_dataframe(table, output_dir, filename)
