# Run from the project root: source("run_qc.R")
# Base R only; no installation or network access required.
options(stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
read_gz <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  read.delim(con, check.names = FALSE, quote = "", comment.char = "")
}
counts <- read_gz("data/raw/counts.tsv.gz")
annot <- read_gz("data/raw/annotation.tsv.gz")
con <- gzfile("data/raw/series_matrix.txt.gz", "rt")
meta_lines <- readLines(con, warn = FALSE)
close(con)
field <- function(key) {
  z <- meta_lines[startsWith(meta_lines, paste0(key, "\t"))]
  if (length(z) != 1L) stop("Missing or ambiguous metadata field: ", key)
  scan(text = z, what = character(), sep = "\t", quote = '"', quiet = TRUE)[-1]
}
ids <- field("!Sample_geo_accession")
titles <- field("!Sample_title")
conditions <- meta_lines[startsWith(meta_lines, "!Sample_characteristics_ch1\t")]
treatment <- conditions[grepl("treatment: day", conditions, fixed = TRUE)]
stopifnot(length(treatment) == 1L)
treatment <- scan(text = treatment, what = character(), sep = "\t", quote = '"', quiet = TRUE)[-1]
samples <- data.frame(sample_id = ids, title = titles, treatment = treatment)
samples <- samples[match(names(counts)[-1], samples$sample_id), ]
stopifnot(!anyNA(samples), all(startsWith(samples$title, "RNA-seq ")))
samples$condition <- ifelse(samples$treatment == "treatment: day 0 of osteogenic induction", "Day0",
  ifelse(samples$treatment == "treatment: day 21 of osteogenic induction", "Day21", NA))
stopifnot(!anyNA(samples$condition), all(table(samples$condition) == 3))
samples$donor <- NA_character_
samples$replicate_independence <- "Not established in GEO metadata"
write.csv(samples, "results/tables/sample_metadata.csv", row.names = FALSE, na = "")
stopifnot(!anyDuplicated(counts$GeneID), !anyDuplicated(annot$GeneID))
idx <- match(counts$GeneID, annot$GeneID)
stopifnot(!anyNA(idx))
annot <- annot[idx, ]
x <- as.matrix(counts[, -1, drop = FALSE])
stopifnot(is.numeric(x), all(is.finite(x)), all(x >= 0), all(x == floor(x)))
rownames(x) <- as.character(counts$GeneID)
lib <- colSums(x)
stopifnot(all(lib > 0))
qc <- data.frame(sample_id = colnames(x), condition = samples$condition,
  assigned_counts = lib, detected_genes = colSums(x > 0),
  zero_fraction = colMeans(x == 0), fraction_top10_genes = apply(x, 2, function(z) sum(sort(z, decreasing = TRUE)[1:10]) / sum(z)))
write.csv(qc, "results/tables/library_qc.csv", row.names = FALSE)
# This descriptive filter is predeclared, not tuned to separation of groups.
cpm <- sweep(x, 2, lib / 1e6, "/")
keep <- rowSums(cpm >= 1) >= 3
stopifnot(sum(keep) >= 2)
logcpm <- log2(cpm[keep, , drop = FALSE] + 1)
write.csv(data.frame(GeneID = counts$GeneID, Symbol = annot$Symbol, retained = keep),
  "results/tables/gene_filter.csv", row.names = FALSE)
colors <- ifelse(samples$condition == "Day0", "#2369A0", "#CF7440")
plot_file <- function(name, code) {
  png(file.path("results/figures", name), width = 1800, height = 1200, res = 180)
  on.exit(dev.off())
  par(mar = c(7, 5, 4, 2))
  force(code)
}
plot_file("01_library_sizes.png", {
  barplot(lib / 1e6, names.arg = samples$title, col = colors, las = 2,
    ylab = "Assigned gene counts (millions)", main = "Library size — not total sequenced reads")
})
plot_file("02_expression_distributions.png", {
  boxplot(logcpm, names = samples$title, col = colors, las = 2,
    ylab = "log2(CPM + 1)", main = "Expression distributions — filtered genes", outline = FALSE)
})
# PCA: 500 most variable filtered genes, centered but not variance-scaled.
vars <- apply(logcpm, 1, var)
top <- head(order(vars, decreasing = TRUE), min(500L, sum(vars > 0)))
stopifnot(length(top) >= 2)
pca <- prcomp(t(logcpm[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
write.csv(data.frame(sample_id = samples$sample_id, condition = samples$condition,
  pca$x), "results/tables/pca_scores.csv", row.names = FALSE)
write.csv(data.frame(component = seq_along(variance), proportion = variance),
  "results/tables/pca_variance.csv", row.names = FALSE)
plot_file("03_pca.png", {
  plot(pca$x[, 1], pca$x[, 2], pch = 19, col = colors, cex = 1.5,
    xlab = sprintf("PC1 (%.1f%%)", 100 * variance[1]),
    ylab = sprintf("PC2 (%.1f%%)", 100 * variance[2]),
    main = "Exploratory PCA — top 500 variable genes")
  text(pca$x[, 1], pca$x[, 2], labels = samples$title, pos = 3, cex = 0.7)
  legend("topright", c("Day 0", "Day 21"), col = c("#2369A0", "#CF7440"), pch = 19)
})
correlation <- cor(logcpm, method = "pearson")
write.csv(correlation, "results/tables/sample_correlations.csv")
plot_file("04_sample_correlations.png", {
  heatmap(correlation, symm = TRUE, scale = "none", margins = c(10, 10),
    labRow = samples$title, labCol = samples$title,
    col = colorRampPalette(c("white", "#2369A0"))(100),
    main = "Pearson correlations — log2(CPM + 1)")
})
# A preselected descriptive marker panel; not a significant-gene list.
markers <- c("RUNX2", "SP7", "ALPL", "COL1A1", "BGLAP", "SPP1", "IBSP")
mi <- which(annot$Symbol %in% markers)
marker <- data.frame(GeneID = counts$GeneID[mi], Symbol = annot$Symbol[mi],
  log2(cpm[mi, , drop = FALSE] + 1), check.names = FALSE)
write.csv(marker, "results/tables/marker_logcpm.csv", row.names = FALSE)
if (length(mi) > 0) plot_file("05_marker_panel.png", {
  par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
  for (i in mi) {
    y <- log2(cpm[i, ] + 1)
    plot(rep(c(1, 2), each = 3), y, pch = 19, col = colors, xaxt = "n",
      xlab = "", ylab = "log2(CPM + 1)", main = annot$Symbol[i], xlim = c(0.7, 2.3))
    axis(1, at = c(1, 2), labels = c("Day0", "Day21"))
  }
})
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
writeLines(c(
  "GSE178679: descriptive RNA-seq QC completed.",
  paste("Genes:", nrow(x), "| retained (CPM >= 1 in >= 3 samples):", sum(keep)),
  "No differential-expression tests were performed.",
  "No samples were automatically removed.",
  "Donor independence, pooling and pairing require clarification before confirmatory inference.",
  "Day0 vs Day21 also differs in culture time; not an isolated treatment effect.",
  "CPM plots and PCA are exploratory, not DESeq2 results.",
  "Count-level QC does not establish read quality, mapping rate or contamination status."
), "results/QC_NOTES.txt")
message("Done. See results/figures and results/tables. No statistical tests performed.")
