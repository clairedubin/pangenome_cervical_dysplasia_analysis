#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(microSLAM)
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(locfdr)

})

# Define command line arguments
option_list <- list(
  make_option(c("--metadata"), type = "character", help = "Path to metadata CSV"),
  make_option(c("--gene_data"), type = "character", help = "Path to gene data CSV"),
  make_option(c("--covariate_equation"), type = "character", help = "GLM formula, e.g. 'y ~ A + B + 1'"),
  make_option(c("--output_dir"), type = "character", help = "Directory to save outputs"),
  make_option(c("--output_file_suffix"), type = "character", help = "Suffix for output filenames")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Read metadata and gene data
metadata <- read.csv(opt$metadata)
gene_data <- read.csv(opt$gene_data)

# Ensure sample_name is first column
metadata <- metadata %>% rename(sample_name = 1)
gene_data <- gene_data %>% rename(sample_name = 1)

# Filter to matching samples
metadata <- metadata %>% filter(sample_name %in% gene_data$sample_name)
gene_data <- gene_data %>% filter(sample_name %in% metadata$sample_name)

# Sort both by sample_name
metadata <- metadata %>% arrange(sample_name)
gene_data <- gene_data %>% arrange(sample_name)

# Set rownames and ensure column order
rownames(metadata) <- metadata$sample_name
rownames(gene_data) <- gene_data$sample_name

# Move sample_name to first column if needed
metadata <- metadata[, c("sample_name", setdiff(names(metadata), "sample_name"))]
sample_by_gene_matrix <- gene_data[, c("sample_name", setdiff(names(gene_data), "sample_name"))]

# Define y (binary outcome)
if (!"outcome" %in% names(metadata)) {
  stop("Metadata must include a column named 'outcome' to generate binary outcome y.")
}
metadata$y <- metadata$outcome

# Step 1: GRM
GRM <- calculate_grm(sample_by_gene_matrix)

# Step 2: Fit baseline GLM
glm_formula <- as.formula(opt$covariate_equation)

# Fit the model with error catching
glm_fit0 <- tryCatch({
  glm(glm_formula, data = metadata, family = "binomial")
}, error = function(e) {
  stop("GLM fitting failed: ", e$message)
})

# Print coefficients
# glm_summary <- base::summary(glm_fit0)

# cat("GLM Fit Coefficients:\n")
# print(glm_summary$coefficients)

# Step 3: Fit GLMM and tau test
glmm_fit <- fit_tau_test(glm_fit0, GRM, species_id = "Lcrispatus", verbose = TRUE)

# Step 4: Permutation test
n_tau <- 100
tautestfit <- run_tau_test(glm_fit0, GRM, n_tau, species_id = "Lcrispatus", tau0 = 1, phi0 = 1, seed = 63)

# Compute p-value
pvalue <- (sum(tautestfit$t >= glmm_fit$t) + 1) / n_tau
cat(sprintf("\nPermutation p-value for tau test: %.5g\n", pvalue))



# Step 5: β test
gene_test_df <- fit_beta(glmm_fit, glm_fit0, GRM, sample_by_gene_matrix, SPA = TRUE)


# glmm_fdr <- locfdr(gene_test_df$SPA_zvalue, pct0 = 0.1, plot = 4)
# gene_test_df$glmm_fdr <- glmm_fdr$fdr

# Sort by FDR
# gene_test_df <- gene_test_df %>% arrange(glmm_fdr)

# # Print selected columns for genes with SPA_zvalue < 0.05
# selected_cols <- c("gene_id", "tau", "cor_to_y", "cor_to_b", "z", "beta", "SPA_pvalue", "pvalue_noadj", "glmm_fdr")

# cat("\nGenes with glmm_fdr < 0.1:\n")
# print(gene_test_df %>%
#   filter(glmm_fdr < 0.1) %>%
#   select(all_of(selected_cols)))


# Calculate number of columns
num_cols <- ncol(sample_by_gene_matrix)

# Decrease bre if number of columns < 1000, per docs
if (num_cols > 1000) {
  bre <- 120
} else {
  bre <- 120 - 5 * ((1000 - num_cols) %/% 100)
}
# print(bre)

run_locfdr_with_retry <- function(zvalues, pct0 = 0.1, plot = 4, initial_df = 7, max_df = 50) {
  df <- initial_df
  repeat {
    warning_occurred <- FALSE
    last_warning <- NULL

    result <- withCallingHandlers(
      tryCatch(
        locfdr(zvalues, pct0 = pct0, plot = plot, df = df, bre = bre),
        error = function(e) {
          message("locfdr failed with error: ", e$message)
          return(NULL)
        }
      ),
      warning = function(w) {
        if (grepl("f\\(z\\) misfit", conditionMessage(w))) {
          warning_occurred <<- TRUE
          last_warning <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      }
    )

    if (!warning_occurred || df >= max_df) {
      converged <- !warning_occurred
      if (!converged) {
        warning("Max df reached (", df, ") but f(z) misfit warning persists.")
      }
      return(list(result = result, converged = converged))
    }

    message("f(z) misfit warning detected (df = ", df, "). Increasing df and retrying...")
    df <- df + 1
  }
}

# Run locfdr with convergence tracking
locfdr_result <- run_locfdr_with_retry(gene_test_df$SPA_zvalue, pct0 = 0.1, plot = 4)

if (!is.null(locfdr_result$result) && locfdr_result$converged) {
  gene_test_df$glmm_fdr <- locfdr_result$result$fdr
  gene_test_df <- gene_test_df %>% arrange(glmm_fdr)
} else {
  warning("locfdr did not converge. `glmm_fdr` column will not be added to gene_test_df.")
}



# Save outputs
output_base <- file.path(opt$output_dir, paste0("microSLAM_", opt$output_file_suffix))
write.csv(gene_test_df, paste0(output_base, ".results.csv"), row.names = FALSE)
write.csv(sample_by_gene_matrix, paste0(output_base, ".sample_by_gene_matrix.csv"), row.names = FALSE)
