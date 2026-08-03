harmonize_annotation_seqlevels <- function(annotation, peak_names) {
  peak_chromosomes <- unique(sub(":.*$", "", peak_names))
  annotation_chromosomes <- GenomeInfoDb::seqlevels(annotation)
  common <- intersect(annotation_chromosomes, peak_chromosomes)

  if (length(common) == 0) {
    stripped <- sub("^chr", "", annotation_chromosomes)
    if (length(intersect(stripped, peak_chromosomes)) > 0) {
      mapping <- stats::setNames(stripped, annotation_chromosomes)
      annotation <- GenomeInfoDb::renameSeqlevels(annotation, mapping)
      annotation_chromosomes <- GenomeInfoDb::seqlevels(annotation)
      common <- intersect(annotation_chromosomes, peak_chromosomes)
    }
  }

  if (length(common) == 0) {
    stop(
      "No chromosome names overlap between EnsDb annotation and peak matrix. ",
      "Peak examples: ", paste(head(peak_chromosomes), collapse = ", "),
      "; annotation examples: ", paste(head(annotation_chromosomes), collapse = ", "),
      call. = FALSE
    )
  }

  GenomeInfoDb::keepSeqlevels(annotation, value = common, pruning.mode = "coarse")
}
