assert_file <- function(path, description = "Required file") {
  if (!file.exists(path)) {
    stop(description, " not found:\n", path, call. = FALSE)
  }
}

write_table <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = row.names)
  invisible(path)
}

save_checkpoint <- function(object, filename) {
  path <- file.path(OBJECT_DIR, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  message("Saved checkpoint: ", path)
  invisible(path)
}

load_checkpoint <- function(filename) {
  path <- file.path(OBJECT_DIR, filename)
  assert_file(path, "Checkpoint")
  readRDS(path)
}

fragment_files_available <- function() {
  file.exists(FRAGMENTS_FILE) && file.exists(FRAGMENTS_INDEX)
}

resolution_slug <- function(x) {
  gsub("\\.", "_", format(x, scientific = FALSE, trim = TRUE))
}
