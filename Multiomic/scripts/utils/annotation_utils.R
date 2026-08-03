split_markers <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

append_token <- function(existing, new_value) {
  ifelse(is.na(existing) | existing == "", new_value, paste(existing, new_value, sep = ";"))
}

count_unique_tokens <- function(x) {
  if (is.na(x) || x == "") return(0L)
  tokens <- strsplit(x, ";", fixed = TRUE)[[1]]
  length(unique(tokens[tokens != ""]))
}

apply_expression_rules <- function(
    object,
    rules_file,
    assay = "RNA",
    layer = "data",
    unassigned_label = "Other") {

  rules <- read.csv(rules_file, stringsAsFactors = FALSE, check.names = FALSE)
  rules <- rules[order(rules$priority), , drop = FALSE]

  all_rule_genes <- unique(unlist(lapply(seq_len(nrow(rules)), function(i) {
    c(
      split_markers(rules$positive_all[i]),
      split_markers(rules$positive_any[i]),
      split_markers(rules$negative_all[i])
    )
  })))

  available_genes <- intersect(all_rule_genes, rownames(object[[assay]]))
  missing_genes <- setdiff(all_rule_genes, available_genes)
  marker_inventory <- data.frame(
    gene = all_rule_genes,
    present = all_rule_genes %in% available_genes
  )

  if (length(available_genes) == 0) {
    stop("None of the annotation marker genes are present in the requested assay.")
  }

  # Explicit assay/layer access avoids ambiguity when the same genes are present
  # in RNA and SCT assays.
  expression_matrix <- SeuratObject::LayerData(
    object = object,
    assay = assay,
    layer = layer,
    features = available_genes,
    cells = colnames(object)
  )

  # LayerData is genes x cells; annotation rules below use cells x genes.
  expression <- as.data.frame(as.matrix(t(expression_matrix)), check.names = FALSE)
  expression <- expression[colnames(object), , drop = FALSE]

  labels <- rep(unassigned_label, nrow(expression))
  names(labels) <- rownames(expression)
  assignment_rule <- rep(NA_character_, nrow(expression))
  names(assignment_rule) <- rownames(expression)
  candidate_rules <- rep("", nrow(expression))
  candidate_labels <- rep("", nrow(expression))
  rule_summaries <- vector("list", nrow(rules))

  for (i in seq_len(nrow(rules))) {
    positive_all <- split_markers(rules$positive_all[i])
    positive_any <- split_markers(rules$positive_any[i])
    negative_all <- split_markers(rules$negative_all[i])

    required_genes <- unique(c(positive_all, negative_all))
    missing_required <- setdiff(required_genes, available_genes)
    available_any <- intersect(positive_any, available_genes)

    if (length(missing_required) > 0 ||
        (length(positive_any) > 0 && length(available_any) == 0)) {
      rule_summaries[[i]] <- data.frame(
        priority = rules$priority[i], rule_name = rules$rule_name[i],
        label = rules$label[i], status = "skipped_missing_markers",
        missing_markers = paste(unique(c(missing_required, setdiff(positive_any, available_any))), collapse = ";"),
        matched_cells = 0, newly_assigned_cells = 0
      )
      next
    }

    positive_threshold <- as.numeric(rules$positive_threshold[i])
    negative_threshold <- as.numeric(rules$negative_threshold[i])
    mask <- rep(TRUE, nrow(expression))

    if (length(positive_all) > 0) {
      m <- as.matrix(expression[, positive_all, drop = FALSE])
      mask <- mask & rowSums(m > positive_threshold) == length(positive_all)
    }
    if (length(available_any) > 0) {
      m <- as.matrix(expression[, available_any, drop = FALSE])
      mask <- mask & rowSums(m > positive_threshold) >= 1
    }
    if (length(negative_all) > 0) {
      m <- as.matrix(expression[, negative_all, drop = FALSE])
      mask <- mask & rowSums(m > negative_threshold) == 0
    }

    candidate_rules[mask] <- append_token(candidate_rules[mask], rules$rule_name[i])
    candidate_labels[mask] <- append_token(candidate_labels[mask], rules$label[i])

    newly_assigned <- mask & labels == unassigned_label
    labels[newly_assigned] <- rules$label[i]
    assignment_rule[newly_assigned] <- rules$rule_name[i]

    rule_summaries[[i]] <- data.frame(
      priority = rules$priority[i], rule_name = rules$rule_name[i],
      label = rules$label[i], status = "evaluated", missing_markers = "",
      matched_cells = sum(mask), newly_assigned_cells = sum(newly_assigned)
    )
  }

  candidate_label_count <- vapply(candidate_labels, count_unique_tokens, integer(1))

  # Explicit barcode-indexed metadata avoids Seurat's "No cell overlap" error.
  annotation_metadata <- data.frame(
    expression_annotation_seed = labels,
    expression_annotation_rule = assignment_rule,
    expression_annotation_candidates = candidate_labels,
    expression_annotation_candidate_rules = candidate_rules,
    expression_annotation_conflict = candidate_label_count > 1,
    expression_annotation_n_candidate_labels = candidate_label_count,
    stringsAsFactors = FALSE,
    row.names = colnames(object)
  )
  annotation_metadata$expression_annotation_seed <- factor(
    annotation_metadata$expression_annotation_seed
  )

  object <- SeuratObject::AddMetaData(object = object, metadata = annotation_metadata)

  list(
    object = object,
    rules = rules,
    rule_summary = do.call(rbind, rule_summaries),
    marker_inventory = marker_inventory,
    missing_genes = missing_genes
  )
}

fill_unassigned_by_cluster <- function(
    object, cluster_col, seed_col, output_col,
    unassigned_label = "Other", min_labeled_cells = 10, min_fraction = 0.70) {

  clusters <- as.character(object[[cluster_col, drop = TRUE]])
  seed_labels <- as.character(object[[seed_col, drop = TRUE]])
  filled_labels <- seed_labels

  cluster_summaries <- lapply(sort(unique(clusters)), function(cluster_id) {
    in_cluster <- clusters == cluster_id
    labeled <- seed_labels[in_cluster]
    labeled <- labeled[labeled != unassigned_label]

    if (length(labeled) == 0) {
      return(data.frame(
        cluster = cluster_id, labeled_cells = 0, consensus_label = NA_character_,
        consensus_fraction = NA_real_, filled_cells = 0, status = "no_expression_seeds"
      ))
    }

    counts <- sort(table(labeled), decreasing = TRUE)
    consensus_label <- names(counts)[1]
    consensus_fraction <- as.numeric(counts[1]) / sum(counts)
    can_fill <- length(labeled) >= min_labeled_cells && consensus_fraction >= min_fraction
    fill_mask <- in_cluster & seed_labels == unassigned_label & can_fill
    filled_labels[fill_mask] <<- consensus_label

    data.frame(
      cluster = cluster_id, labeled_cells = length(labeled),
      consensus_label = consensus_label, consensus_fraction = consensus_fraction,
      filled_cells = sum(fill_mask),
      status = if (can_fill) "eligible" else "insufficient_consensus"
    )
  })

  metadata <- data.frame(value = factor(filled_labels), row.names = colnames(object))
  colnames(metadata) <- output_col
  object <- SeuratObject::AddMetaData(object, metadata = metadata)

  list(object = object, summary = do.call(rbind, cluster_summaries))
}
