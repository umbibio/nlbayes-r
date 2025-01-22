#!/usr/bin/env Rscript

#' TF Inference Simulation Example
#'
#' This script demonstrates how to use nlbayes to perform transcription factor (TF) inference
#' using simulated data. It generates a random network and evidence, then uses the OR-NOR model
#' to infer active TFs.

suppressPackageStartupMessages({
  library(nlbayes)
  library(pROC)
  library(dplyr)
  library(txtplot)
})

#' Generate a random network with specified parameters
#'
#' @param nx Number of transcription factors
#' @param ny Number of target genes
#' @param avg_ntf Average number of TFs that regulate a gene
#' @return A list representing the network structure
gen_network <- function(nx = 125, ny = 2500, avg_ntf = 5) {
  cat("Generating network with parameters:\n")
  cat("  Number of TFs:", nx, "\n")
  cat("  Number of target genes:", ny, "\n")
  cat("  Average TFs per gene:", avg_ntf, "\n")

  # Initialize empty network
  network <- list()
  total_edges <- 0

  # For each target gene, assign random TFs
  for (gene in 1:ny) {
    # Number of TFs follows Poisson distribution
    n_tfs <- rpois(1, avg_ntf)
    if (n_tfs > 0) {
      # Randomly select TFs
      tfs <- sample(1:nx, n_tfs, replace = FALSE)
      # For each selected TF, add this gene to its targets
      for (tf in tfs) {
        tf_name <- sprintf("TF%03d", tf)
        if (is.null(network[[tf_name]])) {
          network[[tf_name]] <- numeric(0)
        }
        network[[tf_name]] <- c(network[[tf_name]], sample(c(-1, 1), 1))
        names(network[[tf_name]])[length(network[[tf_name]])] <- as.character(gene)
        total_edges <- total_edges + 1
      }
    }
  }

  cat("Network generated:\n")
  cat("  Total edges:", total_edges, "\n")
  cat("  Average edges per TF:", round(total_edges/nx, 2), "\n\n")

  return(network)
}

#' Generate evidence based on network and active TFs
#'
#' @param network Network structure
#' @param active_tfs Vector of active TF indices
#' @param tf_target_fraction Fraction of TF targets that become differentially expressed
#' @return Vector of differentially expressed genes
gen_evidence <- function(network, active_tfs, tf_target_fraction = 0.2) {
  cat("Generating evidence:\n")
  cat("  Active TFs:", paste(active_tfs, collapse=", "), "\n")

  # Get all targets of active TFs
  all_targets <- unique(unlist(lapply(network[active_tfs], names)))
  cat("  Total unique targets of active TFs:", length(all_targets), "\n")

  # Randomly select a fraction of these targets
  n_targets <- round(length(all_targets) * tf_target_fraction)
  selected_targets <- sample(all_targets, n_targets, replace = FALSE)
  cat("  Selected targets:", length(selected_targets), "\n")

  # Create evidence vector with signs based on regulation mode
  evidence <- numeric(length(selected_targets))
  names(evidence) <- selected_targets

  for (tf in active_tfs) {
    tf_targets <- names(network[[tf]])
    common_targets <- intersect(tf_targets, selected_targets)
    if (length(common_targets) > 0) {
      evidence[common_targets] <- network[[tf]][common_targets]
    }
  }

  cat("Evidence generated:\n")
  cat("  Total DE genes:", length(evidence), "\n")
  cat("  Up-regulated:", sum(evidence > 0), "\n")
  cat("  Down-regulated:", sum(evidence < 0), "\n\n")

  return(evidence)
}

#' Plot ROC or PR curve using ASCII art
#'
#' @param metric_type Type of curve to plot ('roc' or 'pr')
#' @param true_labels True labels
#' @param scores Predicted scores
#' @param plot_title Title for the plot
plot_metric <- function(metric_type, true_labels, scores, plot_title) {
  if (metric_type == "roc") {
    # Calculate ROC curve
    roc_obj <- roc(true_labels, scores)
    auc_value <- auc(roc_obj)

    # Create ASCII plot
    cat(paste0("\n", plot_title, " (AUC = ", round(auc_value, 2), "):\n"))
    txtplot(1 - roc_obj$specificities,
        roc_obj$sensitivities,
        xlab = "False Positive Rate",
        ylab = "True Positive Rate")

  } else if (metric_type == "pr") {
    # Calculate precision-recall curve manually
    thresholds <- sort(unique(c(-Inf, scores, Inf)))
    precision <- numeric(length(thresholds))
    recall <- numeric(length(thresholds))

    for (i in seq_along(thresholds)) {
      pred_pos <- scores >= thresholds[i]
      tp <- sum(true_labels == 1 & pred_pos)
      fp <- sum(true_labels == 0 & pred_pos)
      fn <- sum(true_labels == 1 & !pred_pos)

      precision[i] <- if (tp + fp == 0) 1 else tp / (tp + fp)
      recall[i] <- if (tp + fn == 0) 0 else tp / (tp + fn)
    }

    # Sort points by recall for proper AUC calculation
    ord <- order(recall)
    recall <- recall[ord]
    precision <- precision[ord]

    # Remove duplicate recall points, keeping max precision
    unique_recalls <- unique(recall)
    unique_precisions <- sapply(unique_recalls, function(r) {
      max(precision[recall == r])
    })

    # Calculate AUC using trapezoidal rule
    auc_value <- 0
    for (i in 2:length(unique_recalls)) {
      delta_r <- unique_recalls[i] - unique_recalls[i-1]
      avg_p <- (unique_precisions[i] + unique_precisions[i-1]) / 2
      auc_value <- auc_value + delta_r * avg_p
    }

    # Create ASCII plot
    cat(paste0("\n", plot_title, " (AUC = ", round(auc_value, 2), "):\n"))
    txtplot(recall, precision,
        xlab = "Recall",
        ylab = "Precision")
  }

  return(auc_value)
}

# Main execution
main <- function() {
  # Set random seed for reproducibility
  set.seed(42)
  cat("Starting simulation with seed 42\n\n")

  # Generate a random network
  network <- gen_network(
    nx = 125,    # total number of transcription factors
    ny = 2500,   # total number of target genes
    avg_ntf = 5  # number of TFs that regulate a gene (average)
  )

  # Randomly select 5 TFs as active
  active_tfs <- sprintf("TF%03d", sample(1:125, 5))
  cat("Selected active TFs:", paste(active_tfs, collapse = ", "), "\n\n")

  # Generate evidence based on the network and active TFs
  evidence <- gen_evidence(
    network = network,
    active_tfs = active_tfs,  # known set of active TFs
    tf_target_fraction = 0.2  # only a fraction of a TF's targets will become diff. expr.
  )

  cat("Creating OR-NOR model...\n")
  # Create and run the OR-NOR model
  inference.model <- ORNOR.inference(network, evidence, n.graphs = 5, verbosity = 2)

  cat("\nSampling posterior distributions...\n")
  # Sample the posterior distributions
  inference.model <- sample.posterior(inference.model, N = 2000, gr.level = 1.1, burnin = TRUE)

  cat("\nPost-processing results...\n")
  # Get inference results
  inference.model <- postprocess.result(inference.model)
  tf.inference <- inference.model$result.info$tf.inference
  cat("Number of TFs in results:", nrow(tf.inference), "\n")

  # Add ground truth column
  cat("Adding ground truth information...\n")
  cat("Number of active TFs:", length(active_tfs), "\n")
  print(str(tf.inference))  # Print structure of tf.inference
  tf.inference$gt_act <- as.numeric(tf.inference$id %in% active_tfs)

  # Print top 10 TFs
  cat("\nTop 10 inferred active TFs:\n")
  print(head(tf.inference[order(-tf.inference$posterior.p), ], 10))

  # Plot ROC and PR curves
  cat("\nPlotting ROC curve...\n")
  roc_auc <- plot_metric("roc", tf.inference$gt_act, tf.inference$posterior.p,
              "Receiver Operating Characteristic (ROC) Curve")
  cat("AUC:", round(roc_auc, 3), "\n")

  cat("\nPlotting PR curve...\n")
  pr_auc <- plot_metric("pr", tf.inference$gt_act, tf.inference$posterior.p,
             "Precision-Recall Curve")
  cat("AUC:", round(pr_auc, 3), "\n")
}

# Run main function
if (!interactive()) {
  main()
}
