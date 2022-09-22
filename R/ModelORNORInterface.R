
compute.enrichment <- function(network, evidence) {

    non.zero.evid <- evidence[evidence != 0]
    n.genes <- length(unique(unlist(lapply(network, names))))
    all.ydeg <- length(non.zero.evid)
    all.ndeg <- n.genes - all.ydeg

    enrichment <- list()
    for (src in names(network)) {
        trg.ydeg <- sum(names(network[[src]]) %in% names(non.zero.evid))
        trg.ndeg <- length(network[[src]]) - trg.ydeg
        contingency <- matrix(
            c(trg.ydeg, trg.ndeg, all.ydeg - trg.ydeg, all.ndeg - trg.ndeg),
            ncol = 2)
        enrichment[[src]] <- fisher.test(
            contingency, alternative = "greater")$p.value
    }

    p.adjust(enrichment, method = "fdr")
}

preprocess.data <- function(
    network, evidence, tf.active.set = c(), verbosity=0) {

    all.tfs <- names(network)
    all.genes <- unique(unlist(lapply(network, names)))

    evidence.length.0 <- length(evidence)
    evidence <- evidence[names(evidence) %in% all.genes]
    evidence.length.1 <- length(evidence)
    removed.evidence <- evidence.length.0 - evidence.length.1

    if (removed.evidence > 0 & verbosity > 0) {
        cat("\nRemoved", removed.evidence,
        "genes from evidence not represented in network")
    }

    tf.active.set <- tf.active.set[tf.active.set %in% all.tfs]

    list(
        network = network,
        evidence = evidence,
        tf.active.set = tf.active.set,
        n.degs = length(evidence[evidence != 0]),
        enrichment = compute.enrichment(network, evidence)
    )
}

ORNOR.inference <- function(
    network, evidence = c(), tf.active.set = c(), N = 2000, s.leniency = 0.1,
    n.graphs = 3, gr.level = 1.10, uniform.t = TRUE, do.sample = FALSE,
    zy = 0., zn = 0., verbosity = 0, burnin = FALSE) {

    if (uniform.t) {
        t_alpha <- 1.
        t_beta <- 1.
    } else if (length(evidence) == 0) {
        t_alpha <- 18.
        t_beta <- 2.
    } else {
        t_alpha <- 2.
        t_beta <- 2.
    }

    input.data <- preprocess.data(network, evidence, tf.active.set, verbosity)

    model <- new(nlbayes::.__C__Rcpp_ModelORNOR)
    model$set.verbosity(verbosity)
    model$set.config(n.graphs, s.leniency, t_alpha, t_beta, zy, zn)

    model$load.data(
        input.data$network,
        input.data$evidence,
        input.data$tf.active.set)

    inference.model <- list(
        input.data = input.data,
        model = model,
        posterior = list(),
        result.info = list(),
        gr = Inf
    )
    if (do.sample)
        inference.model <- sample.posterior(
            inference.model, N = N, gr.level = gr.level, burnin = burnin)

    inference.model
}

transpose.net <- function(net) {
  net.t <- list()
  for (src in names(net)) {
    trg.lst <- net[[src]]
    for (trg in names(trg.lst)) {
      val <- trg.lst[[trg]]
      src.lst <- net.t[[trg]]
      if (is.null(src.lst)) src.lst <- c()

      val <- c(val)
      names(val) <- c(src)
      net.t[[trg]] <- c(src.lst, val)
    }
  }

  net.t
}

explain.ornor <- function(net.t, posterior.p, evidence, lvl=0.5) {
  x <- posterior.p > lvl
  de.explained <- 0
  for (trg in names(evidence)) {
    diff.exp <- evidence[[trg]]
    if (diff.exp == 0) next
    mask <- x[names(net.t[[trg]])]
    y <- net.t[[trg]][mask]
    y <- y[y != 0]
    if (length(y) == 0) next
    ornor.pred <- min(y)
    de.explained <- de.explained + as.integer(ornor.pred == diff.exp)
  }

  de.explained
}

postprocess.result <- function(inference.model, annotation=NULL) {

    if (!is.null(inference.model$posterior$X)) {
        df <- inference.model$posterior$X
        x.posterior.p <- df$posterior.p
        names(x.posterior.p) <- df$uid

        net.t <- transpose.net(inference.model$input.data$network)
        de.explained.80 <- explain.ornor(
            net.t, x.posterior.p, inference.model$input.data$evidence,
            lvl = 0.8)
        de.explained.50 <- explain.ornor(
            net.t, x.posterior.p, inference.model$input.data$evidence,
            lvl = 0.5)
        de.explained.20 <- explain.ornor(
            net.t, x.posterior.p, inference.model$input.data$evidence,
            lvl = 0.2)
        inference.model$result.info$de.explained.80 <- de.explained.80
        inference.model$result.info$de.explained.50 <- de.explained.50
        inference.model$result.info$de.explained.20 <- de.explained.20

        n.degs <- inference.model$input.data$n.degs
        cat("\nDE explained (lvl=0.8): ",
            round(de.explained.80 / n.degs * 100, 2), "%\n")
        cat("DE explained (lvl=0.5): ",
            round(de.explained.50 / n.degs * 100, 2), "%\n")
        cat("DE explained (lvl=0.2): ",
            round(de.explained.20 / n.degs * 100, 2), "%\n\n")

        dff <- data.frame(rank = seq_len(nrow(df)), id = df$uid)
        if (!is.null(annotation))
            dff$symbol <- unlist(annotation[df$uid], use.names = FALSE)
        dff$enrichment <- unlist(
            inference.model$input.data$enrichment[df$uid], use.names = FALSE)
        dff$posterior.p <- df$posterior.p

        dff <- dff[order(-dff$posterior.p), ]
        dff$rank <- seq_len(nrow(dff))
        rownames(dff) <- seq_len(nrow(dff))

        inference.model$result.info$tf.inference <- dff
    }

    inference.model
}

sample.posterior <- function(
    inference.model, N = 2000, gr.level = 1.10, burnin = FALSE) {

    if (burnin) {
        cat("\nInitializing model burn-in ...\n")
        converged <- FALSE
        while (!converged) {
            status <- inference.model$model$sample.posterior(200, 20, 5.0)
            inference.model$gr <- inference.model$model$get.gr()
            converged <- status == 0
            if (status == -1) stop("Interrupt signal received")
        }
        inference.model$model$burn.stats()
        cat("Burn-in complete ...\n")
    }

    n_sampled <- 0
    converged <- FALSE
    i <- 1
    while (n_sampled < N && !converged) {
        n <- min(200 * i, N - n_sampled)
        status <- inference.model$model$sample.posterior(n, 5, gr.level)
        inference.model$gr <- inference.model$model$get.gr()
        converged <- status == 0
        n_sampled <- n_sampled + n
        i <- i + 1
        if (status == -1) stop("Interrupt signal received")
    }

    if (!converged) {
        x <- inference.model$model$get.all.gr()
        n.vars <- length(x)
        x <- x[x > gr.level]
        n.did.not.converge <- length(x)
        cat(
            "\nThere are", n.vars, "random variables in the model. ",
            n.did.not.converge, " of them did not converge.\n")
        if (n.did.not.converge < 20) {
            print(x)
        }
    }

    if (length(inference.model$input.data$evidence) == 0) {
        # We are simulating expression and expect posterior for Y
        inference.model$posterior$H <- inference.model$model$get.posterior("H")
        inference.model$posterior$Y <- inference.model$model$get.posterior("Y")
    }

    if (length(inference.model$input.data$evidence) > 0) {
        df <- inference.model$model$get.posterior("X")
        df$posterior.p <- df$V2.mean
        inference.model$posterior$X <- df

        df <- inference.model$model$get.posterior("T")
        df$posterior.p <- df$V1.mean
        inference.model$posterior$T <- df

        df <- inference.model$model$get.posterior("S")
        tmp <- t(data.frame((strsplit(df$uid, "-->"))))
        rownames(tmp) <- NULL
        df[, c("src", "trg")] <- tmp
        df$posterior.mor <- rowSums(as.matrix(df[, 3:5]) %*% diag(c(-1, 0, 1)))
        inference.model$posterior$S <- df
    }

    inference.model
}
