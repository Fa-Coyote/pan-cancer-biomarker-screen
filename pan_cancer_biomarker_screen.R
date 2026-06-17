# ============================================================
# Pan-Cancer Candidate Biomarker Screening (Demonstration)
# Author: Fatemeh Azizian Farsani
#
# This is a self-contained demonstration script using SIMULATED
# expression data. It illustrates the type of meta-analysis
# approach used to screen for candidate prognostic genes that are
# consistently dysregulated across multiple cancer types -- the
# same general logic behind pan-cancer biomarker discovery work
# I have contributed to (e.g. screening genes such as FAM72B,
# CDK5R1, SERPINE1, and CISH across twelve cancer types).
#
# No real patient, institutional, or third-party code is used
# here -- all data below is randomly generated for illustration.
# Uses base R only (stats package, loaded by default).
# ============================================================

set.seed(42)

# ---- 1. Simulate expression data for several "cancer types" ----
# Each cancer type has a Tumor and a Normal group, measured across
# a panel of genes. A few genes are given a built-in differential
# signal; the rest are pure noise -- mimicking a realistic
# screening scenario where most genes are not informative.

cancer_types <- paste0("Cancer_", 1:6)
genes        <- paste0("Gene_", 1:20)
n_samples    <- 25  # samples per group (Tumor / Normal), per cancer type

# Genes with a true differential signal baked in (for demonstration)
true_signal_genes <- c("Gene_3", "Gene_9", "Gene_14")

simulate_cancer_data <- function(cancer_type) {
  expr <- sapply(genes, function(g) {
    base <- rnorm(2 * n_samples, mean = 5, sd = 1)
    if (g %in% true_signal_genes) {
      base[1:n_samples] <- base[1:n_samples] + 1.5  # tumor up-shift
    }
    base
  })
  data.frame(
    cancer_type = cancer_type,
    sample_id   = paste0(cancer_type, "_S", 1:(2 * n_samples)),
    group       = rep(c("Tumor", "Normal"), each = n_samples),
    expr
  )
}

all_data <- do.call(rbind, lapply(cancer_types, simulate_cancer_data))

# ---- 2. Per-cancer, per-gene differential expression test ----
# Two-sample t-test comparing Tumor vs Normal expression for each
# gene, separately within each simulated cancer type.

run_de_test <- function(df, gene) {
  tumor  <- df[df$group == "Tumor", gene]
  normal <- df[df$group == "Normal", gene]
  test   <- t.test(tumor, normal)
  data.frame(
    cancer_type = unique(df$cancer_type),
    gene        = gene,
    log_fc      = mean(tumor) - mean(normal),
    p_value     = test$p.value
  )
}

de_results <- do.call(rbind, lapply(cancer_types, function(ct) {
  df_ct <- all_data[all_data$cancer_type == ct, ]
  do.call(rbind, lapply(genes, function(g) run_de_test(df_ct, g)))
}))

# ---- 3. Combine evidence across cancer types (meta-analysis) ----
# Fisher's method combines the per-cancer p-values for each gene
# into a single meta p-value, reflecting how consistently the gene
# is dysregulated across multiple independent cancer types.

fishers_method <- function(p_values) {
  p_values <- pmax(p_values, 1e-300)  # avoid log(0)
  stat <- -2 * sum(log(p_values))
  pchisq(stat, df = 2 * length(p_values), lower.tail = FALSE)
}

meta_results <- do.call(rbind, lapply(genes, function(g) {
  p_vals <- de_results$p_value[de_results$gene == g]
  data.frame(
    gene        = g,
    n_cancers   = length(p_vals),
    meta_p      = fishers_method(p_vals),
    mean_log_fc = mean(de_results$log_fc[de_results$gene == g])
  )
}))

meta_results$fdr <- p.adjust(meta_results$meta_p, method = "BH")
meta_results <- meta_results[order(meta_results$meta_p), ]

# ---- 4. Candidate biomarkers ----
cat("Top candidate pan-cancer biomarkers (lowest meta p-value):\n")
print(head(meta_results, 5))

candidates <- meta_results$gene[meta_results$fdr < 0.05]
cat("\nGenes flagged as candidate biomarkers (FDR < 0.05):\n")
print(candidates)
