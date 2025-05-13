mr2.cml.susie.step1 <- function(exposure.ids = NULL,
                                outcome.id = NULL,
                                sample.sizes,
                                beta.exposure.ls = NULL,
                                se.exposure.ls = NULL,
                                beta.outcome.ls = NULL,
                                se.outcome.ls = NULL,
                                use.openGWAS = TRUE) {

  # Check if OpenGWAS data extraction is required
  if (use.openGWAS) {

    # Ensure exposure.ids and outcome.id are provided
    if (is.null(exposure.ids) || is.null(outcome.id)) {
      stop("Please provide 'exposure.ids' and 'outcome.id' when GWAS data is not provided.")
    }

    L <- length(exposure.ids)

    # Check length of sample.sizes
    if (length(sample.sizes) != L) {
      stop("Length of 'sample.sizes' must be equal to the number of exposures (L) when GWAS data is not provided.")
    }

    pval.vec <- rep(NA, L)

    for (i in 1:L) {
      # Get instruments
      exposure.dat <- extract_instruments(exposure.ids[i])

      # Get effects of instruments on outcome
      outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id)

      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure.dat, outcome.dat)

      n <- min(sample.sizes[i], outcome.dat$samplesize.outcome)

      # Perform UVMR-cML
      cML.result <- mr_cML(dat$beta.exposure,
                           dat$beta.outcome,
                           dat$se.exposure,
                           dat$se.outcome,
                           n = n,
                           random_start = 100,
                           random_seed = 1)
      pval.vec[i] <- cML.result$MA_BIC_p
    }

  } else {
    # Ensure user-provided data is available
    if (is.null(beta.exposure.ls) || is.null(se.exposure.ls) || is.null(beta.outcome.ls) || is.null(se.outcome.ls)) {
      stop("Please provide 'beta.exposure.ls', 'se.exposure.ls', 'beta.outcome.ls', 'se.outcome.ls' when GWAS data is provided.")
    }

    # Check that all inputs are lists and of the same length
    if (!is.list(beta.exposure.ls) || !is.list(se.exposure.ls) || !is.list(beta.outcome.ls) || !is.list(se.outcome.ls)) {
      stop("All beta and se inputs must be lists when GWAS data is provided.")
    }

    L <- length(beta.exposure.ls)

    if (length(se.exposure.ls) != L || length(beta.outcome.ls) != L || length(se.outcome.ls) != L) {
      stop("All lists (beta.exposure.ls, se.exposure.ls, beta.outcome.ls, se.outcome.ls) must have the same length.")
    }

    # Check length of sample.sizes
    if (length(sample.sizes) != (L + 1)) {
      stop("Length of 'sample.sizes' must be equal to the number of exposures (L) + 1 (outcome) when GWAS data is provided.")
    }

    pval.vec <- rep(NA, L)

    for (i in 1:L) {
      print(i)
      # Extract vectors for the ith exposure
      beta.exposure <- beta.exposure.ls[[i]]
      se.exposure <- se.exposure.ls[[i]]
      beta.outcome <- beta.outcome.ls[[i]]
      se.outcome <- se.outcome.ls[[i]]

      # Calculate sample size for cML
      n <- min(sample.sizes[i], sample.sizes[L + 1])

      # Perform UVMR-cML using provided data
      cML.result <- mr_cML(beta.exposure,
                           beta.outcome,
                           se.exposure,
                           se.outcome,
                           n = n,
                           random_start = 100,
                           random_seed = 1)
      pval.vec[i] <- cML.result$MA_BIC_p
    }
  }

  return(pval.vec)
}

identify.exposure.subset.idx <- function(step1.res.list, cutoff = NULL) {
  # Ensure step1.res.list is provided
  if (is.null(step1.res.list)) {
    stop("Please provide a list of step 1 p-values (step1.res.list).")
  }

  L <- length(step1.res.list[[1]])

  num.traits <- length(step1.res.list)

  # Default cutoff is Bonferroni on the number of exposures x number of traits
  if (is.null(cutoff)) {
    cutoff = 0.05 / (num.traits * L)
  }

  # Initialize an empty vector to store indices of significant p-values
  sig.idx <- c()

  # Loop through each element in the list to find indices of significant p-values
  for (j in 1:num.traits) {
    # Find indices in the current sublist that are less than the cutoff
    subset.idx <- which(step1.res.list[[j]] < cutoff)

    # Append these indices to the vector of significant indices
    sig.idx <- union(sig.idx, subset.idx)
  }
  return(sort(sig.idx))
}

harmonize.mr2.data <- function(exposure.ids.subset,
                               outcome.id.list,
                               sample.sizes.subset) {

  # Ensure exposure.ids.subset and outcome.id.list are provided
  if (is.null(exposure.ids.subset) || is.null(outcome.id.list)) {
    stop("Please provide 'exposure.ids.subset' and 'outcome.id.list' when GWAS data is not provided.")
  }

  # Number of exposures in the subset
  L.star <- length(exposure.ids.subset)

  # Check length of sample.sizes.subset
  if (length(sample.sizes.subset) != L.star) {
    stop("Length of 'sample.sizes.subset' must be equal to the number of exposures (L.star).")
  }

  num.traits <- length(outcome.id.list)

  mvdat.list <- vector("list", length = num.traits)

  outcome.sample.sizes <- rep(NA, num.traits)

  for (j in 1:num.traits) {
    # Get instruments jointly for the subset of exposures
    exposure.dat <- mv_extract_exposures(exposure.ids.subset)

    # Get effects of instruments on an outcome
    outcome.dat <- extract_outcome_data(snps = exposure.dat$SNP, outcomes = outcome.id.list[[j]])

    outcome.sample.sizes[j] <- min(outcome.dat$samplesize.outcome)
    # Harmonise the exposure and outcome data
    mvdat.list[[j]] <- mv_harmonise_data(exposure.dat, outcome.dat)
  }
  # Find common SNPs across all outcome datasets
  common.snps <- Reduce(intersect, lapply(mvdat.list, function(mvdat) rownames(mvdat$exposure_beta)))

  # Filter each dataset to keep only common SNPs
  for (j in 1:num.traits) {
    mvdat <- mvdat.list[[j]]

    # Get indices of common SNPs in the jth mvdat
    common.snps.idx <- match(common.snps, rownames(mvdat$exposure_beta), nomatch = 0)

    # Filter exposure matrices and outcome arrays
    mvdat$exposure_beta <- mvdat$exposure_beta[common.snps.idx, , drop = FALSE]
    mvdat$exposure_pval <- mvdat$exposure_pval[common.snps.idx, , drop = FALSE]
    mvdat$exposure_se <- mvdat$exposure_se[common.snps.idx, , drop = FALSE]

    mvdat$outcome_beta <- mvdat$outcome_beta[common.snps.idx]
    mvdat$outcome_pval <- mvdat$outcome_pval[common.snps.idx]
    mvdat$outcome_se <- mvdat$outcome_se[common.snps.idx]

    # Update the list with filtered data
    mvdat.list[[j]] <- mvdat
  }
  return(list(mvdat.list = mvdat.list, exposure.sample.sizes = sample.sizes.subset,
              outcome.sample.sizes = outcome.sample.sizes))
}

mr2.cml.susie.step2 <- function(mr2dat, outcome.idx, cutoff = 5e-8) {

  # Number of exposures in the subset
  L.star <- dim(mr2dat$mvdat.list[[1]]$exposure_beta)[2]

  uvmr.ls <- list()

  mvdat <- mr2dat$mvdat.list[[outcome.idx]]

  # Focus only on QTL for each exposure based on p-values
  for (i in 1:L.star) {
    uvmr.ls[[i]] <- which(mvdat$exposure_pval[, i] < cutoff)
  }

  # Initialize result vectors and lists
  pval.vec <- rep(NA, L.star)
  theta.vec <- rep(NA, L.star)
  uvmr.invalid.ls <- list()

  for (i in 1:L.star) {
    print(i)
    n <- min(mr2dat$exposure.sample.sizes[i], mr2dat$outcome.sample.sizes[outcome.idx])

    # Perform UVMR-cML
    cML.result <- mr_cML(mvdat$exposure_beta[uvmr.ls[[i]], i],
                         mvdat$outcome_beta[uvmr.ls[[i]]],
                         mvdat$exposure_se[uvmr.ls[[i]], i],
                         mvdat$outcome_se[uvmr.ls[[i]]],
                         n = n,
                         random_start = 100,
                         random_seed = 1)

    pval.vec[i] <- cML.result$MA_BIC_p
    theta.vec[i] <- cML.result$MA_BIC_theta

    if (length(cML.result$BIC_invalid) != 0) {
      uvmr.invalid.ls[[i]] <- uvmr.ls[[i]][cML.result$BIC_invalid]
    } else {
      uvmr.invalid.ls[[i]] <- cML.result$BIC_invalid
    }
  }

  # IV is considered as invalid if invalid in any one of the exposures
  rm.idx <- sort(Reduce(union, uvmr.invalid.ls))
  return(list(mvdat = mvdat, invalid.idx = rm.idx, theta.vec = theta.vec))
}

# MVMR-cML-SuSiE
mr2.cml.susie.step3 <- function(mvdat, invalid.idx, theta.vec, rho.mat,
                                L = 10, max.iter = 200, tol = 1e-10) {

  m.star <- dim(mvdat$exposure_beta)[1]
  L.star <- dim(mvdat$exposure_beta)[2]

  # Check if the dimension of rho.mat is (L + 1) x (L + 1)
  if (ncol(rho.mat) != (L.star + 1) || nrow(rho.mat) != (L.star + 1)) {
    stop(paste("The dimension of 'rho.mat' should be", (L.star + 1), "x", (L.star + 1),
               "but it is", nrow(rho.mat), "x", ncol(rho.mat), "."))
  }

  X <- mvdat$exposure_beta
  Y <- mvdat$outcome_beta
  SigmaX <- mvdat$exposure_se
  seY <- mvdat$outcome_se

  # Create a list of genetic correlation corresponding to individual i
  sei.list <- list()

  for (i in 1:m.star) {
    sei.list[[i]] <- c(SigmaX[i,], seY[i])
  }

  # Create a list of variance corresponding to individual i
  Sigmai.list <- list()

  for (i in 1:m.star) {
    Sigmai.list[[i]] <- diag(sei.list[[i]]) %*% rho.mat %*% diag(sei.list[[i]])
  }

  # Create a list of profile likelihood variance (denominator) corresponding to individual i
  new.sigma.vec <- rep(NA, m.star)

  gamma.vec <- c(theta.vec, -1)
  for (i in 1:m.star) {
    new.sigma.vec[i] <- sqrt(t(gamma.vec) %*% Sigmai.list[[i]] %*% gamma.vec)
  }

  # Transform exposure and outcome using profile likehood weights
  Xnew <- X / new.sigma.vec
  Ynew <- Y / new.sigma.vec

  k <- length(invalid.idx)

  # Remove invalid IVs previously identified in step 2 if present
  if (k > 0) {
    X.sub <- Xnew[-invalid.idx,]
    Y.sub <- Ynew[-invalid.idx]
  } else {
    X.sub <- Xnew
    Y.sub <- Ynew
  }

  res <- susie(X.sub, Y.sub, L = 10, max = max.iter, intercept = FALSE)

  theta.vec <- unname(coef(res))[-1]
  prev.theta <- theta.vec
  iter <- 0
  diff <- 1e300
  thres <- tol
  while (diff > thres) {
    iter <- iter + 1

    gamma.vec <- c(unname(coef(res)[-1]), -1)

    for (i in 1:m.star) {
      new.sigma.vec[i] <- sqrt(t(gamma.vec) %*% Sigmai.list[[i]] %*% gamma.vec)
    }

    Xnew <- X / new.sigma.vec
    Ynew <- Y / new.sigma.vec

    if (k > 0) {
      X.sub <- Xnew[-invalid.idx,]
      Y.sub <- Ynew[-invalid.idx]
    } else {
      X.sub <- Xnew
      Y.sub <- Ynew
    }

    res <- susie(X.sub, Y.sub, L = L, max = max.iter, intercept = FALSE)
    theta.vec <- unname(coef(res))[-1]
    diff <- sum((prev.theta - theta.vec)^2)
    prev.norm <- diff
    if (diff > 0) {
      prev.theta <- theta.vec
    } else if (diff < 0) {
      iter <- iter - 1
      break
    }
    prev.norm <- diff
  }

  return(res)
}

# Function for obtaining the posterior mean of MESuSiE
get_posterior_mean_mesusie <- function(mesusie_object, prior_tol = 1e-9) {
  # Number of SNPs
  p <- mesusie_object$nSNP

  # Number of ancestries
  N_ancestry <- mesusie_object$nancestry

  # Initialize a matrix to store posterior means for each ancestry
  posterior_means <- matrix(0, nrow = p, ncol = N_ancestry)

  # Calculate posterior means for each ancestry
  for (l_index in seq_len(mesusie_object$L)) {
    if (all(abs(mesusie_object$V[[l_index]]) > prior_tol)) {
      alpha <- mesusie_object$alpha[[l_index]]
      mu1 <- mesusie_object$mu1[[l_index]]

      for (config_index in seq_along(mesusie_object$column_config)) {
        config <- mesusie_object$column_config[[config_index]]
        for (ancestry_position in seq_along(config)) {
          actual_ancestry <- config[ancestry_position]

          # Check if mu1[[config_index]] is a vector or a matrix
          if (is.vector(mu1[[config_index]])) {
            posterior_means[, actual_ancestry] <- posterior_means[, actual_ancestry] +
              alpha[, actual_ancestry] * mu1[[config_index]]
          } else if (is.matrix(mu1[[config_index]])) {
            posterior_means[, actual_ancestry] <- posterior_means[, actual_ancestry] +
              alpha[, actual_ancestry] * mu1[[config_index]][, ancestry_position]
          } else {
            stop("Unexpected structure in mu1.")
          }
        }
      }
    }
  }
  return(posterior_means)
}

mr2.cml.susie.step4 <- function(mvdat.list, invalid.idx, theta.vec.list, rho.mat,
                                L = 10, max.iter = 200, tol = 1e-10,
                                exposure.names = NULL, outcome.names = NULL) {
  m.star <- dim(mvdat.list[[1]]$exposure_beta)[1]
  L.star <- dim(mvdat.list[[1]]$exposure_beta)[2]

  num.traits <- length(mvdat.list)

  # Optional name checks
  if (!is.null(exposure.names)) {
    if (length(exposure.names) != L.star) {
      stop("Length of exposure.names (", length(exposure.names),
           ") must match the number of exposures (", L.star, ").")
    }
  }
  if (!is.null(outcome.names)) {
    if (length(outcome.names) != num.traits) {
      stop("Length of outcome.names (", length(outcome.names),
           ") must match the number of outcomes (", num.traits, ").")
    }
  }

  # Initialize storage for computed values
  X <- mvdat.list[[1]]$exposure_beta
  Y.list <- lapply(mvdat.list, function(mv) mv$outcome_beta)
  SigmaX.list <- lapply(mvdat.list, function(mv) mv$exposure_se)
  seY.list <- lapply(mvdat.list, function(mv) mv$outcome_se)

  # Calculate sei.list and Sigmai.list for each trait
  sei.list <- vector("list", length = num.traits)
  Sigmai.list <- vector("list", length = num.traits)
  for (j in 1:num.traits) {
    sei.list[[j]] <- lapply(1:m.star, function(i) c(SigmaX.list[[j]][i, ], seY.list[[j]][i]))
    Sigmai.list[[j]] <- lapply(1:m.star, function(i) {
      diag(sei.list[[j]][[i]]) %*% rho.mat %*% diag(sei.list[[j]][[i]])
    })
  }

  # Compute profile likelihood variance and transform X, Y for each trait
  new.sigma.vec.list <- vector("list", length = num.traits)
  Xnew.list <- vector("list", length = num.traits)
  Ynew.list <- vector("list", length = num.traits)

  # Initial theta.mat from first fit (e.g., from mvmr.cml.susie.step3 results or initialized values)
  theta.mat <- matrix(unlist(theta.vec.list), ncol = num.traits)

  prev.theta <- theta.mat
  iter <- 0
  diff <- 1e300
  
  # IMS

  while (diff > tol) {
    iter <- iter + 1

    new.sigma.vec.list <- list()
    for (j in 1:num.traits) {
      gamma.vec <- c(theta.mat[, j], -1)
      new.sigma.vec <- sapply(1:m.star, function(i) {
        sqrt(t(gamma.vec) %*% Sigmai.list[[j]][[i]] %*% gamma.vec)
      })
      new.sigma.vec.list[[j]] <- new.sigma.vec
    }

    Xnew.list <- list()
    Ynew.list <- list()
    for (j in 1:num.traits) {
      Xnew.list[[j]] <- X / new.sigma.vec.list[[j]]
      Ynew.list[[j]] <- Y.list[[j]] / new.sigma.vec.list[[j]]
    }
    
    # Process invalid indices for all traits
    Xsub.list <- list()
    Ysub.list <- list()
    for (j in 1:num.traits) {
      if (length(invalid.idx) > 0) {
        Xsub.list[[j]] <- Xnew.list[[j]][-invalid.idx, ]
        Ysub.list[[j]] <- Ynew.list[[j]][-invalid.idx]
      } else {
        Xsub.list[[j]] <- Xnew.list[[j]]
        Ysub.list[[j]] <- Ynew.list[[j]]
      }

      # Standardize X and Y
      Xsub.list[[j]] <- scale(Xsub.list[[j]])
      Ysub.list[[j]] <- scale(Ysub.list[[j]])
    }

    rho.list <- list()
    theta.list <- list()
    
    # Prepare rho.list and theta.list for meSuSie_core
    for (j in 1:num.traits) {

      # Use the outcome name if provided, else default
      outcome.name <- if (!is.null(outcome.names)) outcome.names[j] else paste0("Outcome", j)

      n <- dim(Xsub.list[[j]])[1]
      Z <- t(Xsub.list[[j]]) %*% Ysub.list[[j]] / sqrt(n)

      X_tX <- t(Xsub.list[[j]]) %*% Xsub.list[[j]]
      rho.list[[outcome.name]] <- cov2cor(X_tX)

      # Construct the SNP names based on whether exposure names are provided
      snp.names <- if (!is.null(exposure.names)) exposure.names else paste0("Exposure", 1:L.star)

      theta.list[[outcome.name]] <- data.frame(
        SNP = paste0("Exposure", 1:L.star),
        Beta = rep(0, L.star),
        Se = rep(0.001825742, L.star),
        Z = Z,
        N = rep(m.star - length(invalid.idx), L.star)
      )
      theta.list[[outcome.name]]$Beta <- theta.list[[outcome.name]]$Z * theta.list[[outcome.name]]$Se
    }

    res <- meSuSie_core(rho.list, theta.list, L = 10)
    theta.mat <- get_posterior_mean_mesusie(res)

    diff <- sum((prev.theta - theta.mat)^2, na.rm = TRUE)
    if (diff > 0) {
      prev.theta <- theta.mat
    } else {
      break
    }
  }

  for (l in 1:L) {
    rownames(res$alpha[[l]]) <- snp.names
  }

  config.names <- colnames(res$pip_config)

  for (l in 1:L) {
    colnames(res$alpha[[l]]) <- config.names
  }

  return(res)
}
