invisible({
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Expand shortend inst/ pathnames (that were too long for tar)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Search-replace patterns
map <- c(
  "_BC_" = "BackgroundCorrection",
  "_CE_" = "ChromosomeExplorer",
  "_GCN_" = "GcContentNormalization",
  "_SN_" = "SmoothNormalization",
  "_DAR_" = "downloadAnnotationData"
);
for (kk in seq_along(map)) {
  pattern <- names(map)[kk];
  pathnames <- list.files("inst", pattern=pattern, full.names=TRUE, recursive=TRUE);
  if (length(pathnames) > 0L) {
    pathnamesD <- gsub(pattern, map[kk], pathnames);
    file.rename(pathnames, pathnamesD);
  }
}
})

