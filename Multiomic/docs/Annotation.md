# Expression-based annotation

The repository preserves the original manually curated cluster map, but the primary
annotation workflow assigns labels at the individual-cell level from normalized RNA
expression using `config/stage11_expression_annotation_rules.csv`.

## Rule behavior

Rules are evaluated in ascending `priority`. Each rule can specify:

- `positive_all`: every listed marker must exceed the positive threshold.
- `positive_any`: at least one available marker must exceed the positive threshold.
- `negative_all`: none of the listed markers may exceed the negative threshold.

The first matching rule supplies `expression_annotation_seed`; all candidate rules and
labels are also retained so conflicts remain visible.

## Metadata produced

- `expression_annotation_seed`
- `expression_annotation_rule`
- `expression_annotation_candidates`
- `expression_annotation_candidate_rules`
- `expression_annotation_conflict`
- `expression_annotation_n_candidate_labels`
- `expression_annotation_cluster_filled`
- `final_annotation`

The rule engine explicitly accesses the RNA `data` layer with `LayerData()`. This avoids
ambiguity when RNA and SCT contain the same gene names. Metadata are added with explicit
cell-barcode row names to preserve cell alignment.

## Manual annotation provenance

`config/original_manual_cluster_map.csv` retains the original cluster-number-to-cell-type
mapping. It is **not automatically applied** after QC/reclustering because numeric cluster
IDs are not guaranteed to remain biologically equivalent after the pipeline changes.
