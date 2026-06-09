# Validation script: run full tAge pipeline on example data
# and save intermediate matrices for Python cross-check.

library(tAge)

# ---- Load example data ----
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()

cat("Example data loaded:\n")
cat("  Genes:", nrow(expr_data), "\n")
cat("  Samples:", ncol(expr_data), "\n")

# ---- Create ExpressionSet ----
eset <- make_ExpressionSet(expr_data, meta_data, verbose = FALSE)

# ---- Preprocessing (mouse, Ensembl IDs) ----
processed <- tAge_preprocessing(
  eset,
  species = "mouse",
  gene_mapping_type = "Ensembl",
  control_group_column = "Genotype",
  control_group_label = "WT",
  count_threshold = 10,
  percent_threshold = 20,
  verbose = TRUE
)

# Save intermediate matrices
cat("\nSaving intermediate matrices...\n")

# scaled_diff (the one we need for the EN_Mortality_Multispecies_Multitissue_scaleddiff model)
scaled_diff_expr <- as.data.frame(Biobase::exprs(processed$scaled_diff))
write.csv(scaled_diff_expr, "data/test/r_scaled_diff_matrix.csv")

# Also save scaled (without control subtraction) for reference
scaled_expr <- as.data.frame(Biobase::exprs(processed$scaled))
write.csv(scaled_expr, "data/test/r_scaled_matrix.csv")

# Save metadata with predictions
meta_out <- as.data.frame(Biobase::pData(eset))
write.csv(meta_out, "data/test/r_metadata.csv")

# ---- Prediction ----
model_paths <- list(
  scaled_diff = "models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
)

# Set Python path for reticulate
Sys.setenv(RETICULATE_PYTHON = Sys.which("python3"))

results <- predict_tAge(
  tAge_eset = processed,
  model_paths = model_paths,
  species = "mouse",
  mode = "EN"
)

write.csv(results, "data/test/r_predictions.csv")
cat("\nR predictions saved to data/test/r_predictions.csv\n")
print(head(results))
