# TMM normalization script for Python-R bridge
# Usage: Rscript tmm_normalize.R input_counts.csv output_logcpm.csv
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read counts matrix (genes x samples)
counts <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Convert to numeric matrix
counts_mat <- as.matrix(sapply(counts, as.numeric))
rownames(counts_mat) <- rownames(counts)

cat("Input:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")

# TMM normalization
dge <- DGEList(counts = counts_mat)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Write output
write.csv(as.data.frame(logcpm), output_file)
cat("TMM normalization complete. Output:", output_file, "\n")
