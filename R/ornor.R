#' Create an OR-NOR Bayesian Network Model
#'
#' Creates a model for inferring transcription factor activities from gene expression data
#' using an OR-NOR Bayesian network approach.
#'
#' @param network List where names are TF IDs and values are named numeric vectors
#'   mapping target gene names to regulation modes (-1 for repression, 1 for activation)
#' @param evidence Named numeric vector of gene expression states (-1 for down, 1 for up)
#' @param n_graphs Integer number of parallel graphs for convergence assessment (default: 3)
#' @param uniform_prior Logical, whether to use uniform prior for theta parameter (default: FALSE)
#' @param active_tfs Character vector of TFs known to be active (optional)
#'
#' @return An ornor model object
#' @export
ornor <- function(network, evidence = NULL, n_graphs = 3, uniform_prior = FALSE, active_tfs = NULL) {
    # Preprocess data
    if (is.null(evidence)) evidence <- numeric()
    if (is.null(active_tfs)) active_tfs <- character()
    
    # Create model using existing interface
    model <- ORNOR.inference(
        network = network,
        evidence = evidence,
        tf.active.set = active_tfs,
        uniform.t = uniform_prior,
        n.graphs = n_graphs
    )
    
    # Add class for S3 methods
    class(model) <- c("ornor", class(model))
    model
}

#' Fit an OR-NOR Model
#'
#' Fits the model using MCMC sampling to infer TF activities.
#'
#' @param model An ornor model object
#' @param n_samples Integer number of samples to draw (default: 2000)
#' @param gelman_rubin Numeric Gelman-Rubin convergence criterion (default: 1.1)
#' @param burnin Logical, whether to perform burn-in phase (default: TRUE)
#'
#' @return The fitted model
#' @export
fit.ornor <- function(model, n_samples = 2000, gelman_rubin = 1.1, burnin = TRUE) {
    # Sample posterior using existing interface
    model <- sample.posterior(
        model,
        N = n_samples,
        gr.level = gelman_rubin,
        burnin = burnin
    )
    
    # Mark as fitted
    attr(model, "fitted") <- TRUE
    model
}

#' Get OR-NOR Model Results
#'
#' Retrieves the inference results from a fitted OR-NOR model.
#'
#' @param model A fitted ornor model object
#'
#' @return A data frame containing TF inference results with columns:
#'   - TF_id: TF identifier
#'   - X: Inferred activity score
#'   - posterior_p: Posterior probability of activity
#' @export
get_results.ornor <- function(model) {
    if (!isTRUE(attr(model, "fitted")))
        stop("Model must be fitted before getting results")
    
    # Get results using existing interface and ensure it's a data frame
    model <- postprocess.result(model)
    results <- as.data.frame(model$result.info$tf.inference)
    
    # Rename columns to match Python API
    names(results)[names(results) == "id"] <- "TF_id"
    names(results)[names(results) == "posterior.p"] <- "posterior_p"
    
    results
}

#' @export
print.ornor <- function(x, ...) {
    cat("OR-NOR Bayesian Network Model\n")
    cat("Network size:", length(x$network), "TFs\n")
    if (!is.null(x$evidence))
        cat("Evidence size:", length(x$evidence), "genes\n")
    if (isTRUE(attr(x, "fitted")))
        cat("Status: Fitted\n")
    else
        cat("Status: Not fitted\n")
}

# Register S3 methods
#' @export
fit <- function(model, ...) UseMethod("fit")

#' @export
get_results <- function(model, ...) UseMethod("get_results")
