# Static repository validation. Does not execute the biological analysis.
# Run from repository root:
#   Rscript scripts/validate_repository.R

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "README.md")) &&
        dir.exists(file.path(current, "scripts")) &&
        dir.exists(file.path(current, "config"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate repository root.")
}
root <- find_project_root()
validation_dir <- file.path(root, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "README.md", ".gitignore", "config/analysis_config.R",
  "config/paths_example.R", "config/stage11_expression_annotation_rules.csv",
  "config/original_manual_cluster_map.csv", "config/marker_panels.csv",
  "scripts/00_setup.R", "scripts/run_all.R", "scripts/07_expression_annotation.R",
  "scripts/08_annotation_validation.R", "docs/Pipeline.md", "docs/Annotation.md"
)

rows <- lapply(required_files, function(rel) {
  data.frame(check = paste0("file:", rel), status = if (file.exists(file.path(root, rel))) "PASS" else "FAIL", detail = "")
})

r_files <- list.files(file.path(root, "scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (path in r_files) {
  ok <- TRUE; detail <- ""
  tryCatch(parse(file = path), error = function(e) { ok <<- FALSE; detail <<- conditionMessage(e) })
  rows[[length(rows)+1]] <- data.frame(
    check = paste0("R_parse:", sub(paste0("^", root, "/?"), "", gsub("\\\\", "/", path))),
    status = if (ok) "PASS" else "FAIL", detail = detail
  )
}

rules <- read.csv(file.path(root, "config", "stage11_expression_annotation_rules.csv"), check.names = FALSE)
required_rule_cols <- c("priority", "rule_name", "label", "positive_all", "positive_any", "negative_all", "positive_threshold", "negative_threshold")
missing_cols <- setdiff(required_rule_cols, colnames(rules))
rows[[length(rows)+1]] <- data.frame(
  check = "annotation_rule_columns", status = if (length(missing_cols)==0) "PASS" else "FAIL",
  detail = paste(missing_cols, collapse = ";")
)

results <- do.call(rbind, rows)
write.csv(results, file.path(validation_dir, "repository_validation.csv"), row.names = FALSE)
summary_lines <- c(
  paste("Checks:", nrow(results)),
  paste("Passed:", sum(results$status == "PASS")),
  paste("Failed:", sum(results$status == "FAIL")),
  "",
  apply(results, 1, function(x) paste(x[["status"]], x[["check"]], x[["detail"]]))
)
writeLines(summary_lines, file.path(validation_dir, "repository_validation.txt"))
print(results)
if (any(results$status == "FAIL")) quit(status = 1)
