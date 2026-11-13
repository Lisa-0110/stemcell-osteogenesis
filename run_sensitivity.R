# Run from the project root. Uses base R only.
source("run_qc.R")
positive <- x[rowSums(x > 0) == ncol(x), , drop = FALSE]
stopifnot(nrow(positive) > 0)
sf <- apply(positive / exp(rowMeans(log(positive))), 2, median)
stopifnot(all(is.finite(sf)), all(sf > 0))
mor <- log2(sweep(x, 2, sf, "/") + 1)
day0 <- samples$condition == "Day0"
day21 <- samples$condition == "Day21"
delta <- function(z) rowMeans(z[, day21, drop = FALSE]) -
  rowMeans(z[, day0, drop = FALSE])
write.csv(data.frame(sample_id = colnames(x), size_factor = sf),
  "results/tables/size_factors.csv", row.names = FALSE)
write.csv(data.frame(Symbol = annot$Symbol[mi],
  delta_log2_cpm = delta(log2(cpm[mi, , drop = FALSE] + 1)),
  delta_log2_mor = delta(mor[mi, , drop = FALSE])),
  "results/tables/marker_sensitivity.csv", row.names = FALSE)
filtered <- mor[keep, , drop = FALSE]
top_mor <- head(order(apply(filtered, 1, var), decreasing = TRUE), 500)
p_mor <- prcomp(t(filtered[top_mor, , drop = FALSE]))
v_mor <- p_mor$sdev^2 / sum(p_mor$sdev^2)
write.csv(data.frame(method = rep(c("CPM", "Median-of-ratios"), each = 2),
  PC = rep(1:2, 2), percent = 100 * c(variance[1:2], v_mor[1:2])),
  "results/tables/pca_sensitivity.csv", row.names = FALSE)
plot_file("06_pca_sensitivity.png", {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  for (j in 1:2) {
    p <- list(pca, p_mor)[[j]]
    v <- list(variance, v_mor)[[j]]
    plot(p$x[, 1], p$x[, 2], col = colors, pch = 19, cex = 1.4,
      main = c("CPM", "Median-of-ratios")[j],
      xlab = sprintf("PC1 (%.1f%%)", 100 * v[1]),
      ylab = sprintf("PC2 (%.1f%%)", 100 * v[2]))
    legend("topright", c("Day 0", "Day 21"), col = c("#2369A0", "#CF7440"),
      pch = 19, bty = "n", horiz = TRUE)
  }
})
message("Sensitivity check saved. No statistical tests performed.")
