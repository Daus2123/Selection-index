# Advanced analysis extensions for Modules 2, 4, and 5
#
# This file collects additional app.R analysis families:
#   1) D2 genetic diversity / clustering / PCA / correlation
#   2) Expanded MET stability wrappers from metan
#   3) Multi-trait selection indices from metan
#
# The functions are intentionally UI-free. They can be called from app.R server
# logic, export logic, or future modules without changing the existing code paths.

# Mean comparisons retain fixed genotypes. This separate REML fit estimates
# generalized broad-sense H2 from joint prediction-error variances of differences.
si_lpsi_lmm_statistics <- function(data, trait, selection_proportion) {
  d <- data[is.finite(data[[trait]]), , drop = FALSE]
  d$ID <- droplevels(factor(d$ID)); d$Rep <- droplevels(factor(d$Rep))
  mean_value <- mean(d[[trait]])
  counts <- table(d$ID)
  harmonic <- if (length(counts)) 1 / mean(1 / as.numeric(counts)) else NA_real_
  intensity <- dnorm(qnorm(1 - selection_proportion)) / selection_proportion
  out <- data.frame(Trait = trait, Mean = mean_value, CV_pct = NA_real_,
    MS_genotype = NA_real_, MS_error = NA_real_, Harmonic_replication = harmonic,
    Genotypic_variance = NA_real_, Error_variance = NA_real_, Phenotypic_variance = NA_real_,
    Broad_sense_H2 = NA_real_, Selection_intensity_pct = selection_proportion * 100,
    Selection_intensity_k = intensity, Genetic_advance = NA_real_, Genetic_advance_pct_mean = NA_real_,
    CV_model_residual_variance = NA_real_, Mean_PEV_difference = NA_real_,
    Heritability_method = "Cullis genotype-difference H2 (separate REML model)",
    Heritability_model = "Trait ~ 1 + (1 | ID) + (1 | Rep)",
    CV_status = "Model failed", H2_status = "Model failed", Note = "", stringsAsFactors = FALSE)
  notes <- c("Residual CV uses Trait ~ ID + (1 | Rep). H2 uses a separate random-genotype model; mean comparisons remain fixed-genotype.",
    "Genetic advance is not calculated from generalized H2.")
  warnings <- character()
  capture_fit <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  check_fit <- function(fit) {
    messages <- fit@optinfo$conv$lme4$messages
    messages <- messages[!grepl("boundary.*singular", messages)]
    if (length(messages)) stop(paste(messages, collapse = "; "))
    code <- fit@optinfo$conv$opt
    if (length(code) && any(code != 0)) stop("Optimizer did not converge.")
  }
  # Diagnose exact boundary datasets before asking an optimizer to estimate
  # vanishing residual variance (common with rounded/replicated score values).
  if (nlevels(d$ID) >= 2 && nlevels(d$Rep) >= 2 && nrow(d) > nlevels(d$ID)) {
    boundary <- stats::lm(as.formula(paste0(backtick_name(trait), " ~ ID + Rep")), d)
    tolerance <- 1e-10 * max(1, stats::sd(d[[trait]]))
    if (df.residual(boundary) > 0 && max(abs(residuals(boundary))) <= tolerance) {
      out$CV_model_residual_variance <- 0
      out$CV_pct <- if (abs(mean_value) > .0001) 0 else NA_real_
      out$CV_status <- if (is.finite(out$CV_pct)) "Zero residual variation" else "Zero mean"
      out$H2_status <- if (diff(range(d[[trait]])) <= tolerance) "Constant trait" else "Zero residual variation"
      out$Note <- paste(c(notes, paste(out$H2_status,
        "makes the REML/Cullis variance ratio non-regular. H2 is not forced to 0 or 1; verify rounding and independent replication.")), collapse = " | ")
      return(out)
    }
  }
  # Fit on a standardized response; return all variances in original trait units.
  response_scale <- stats::sd(d[[trait]])
  if (!is.finite(response_scale) || response_scale <= 0) response_scale <- 1
  scaled <- d
  scaled[[trait]] <- (d[[trait]] - mean_value) / response_scale
  fit_with_retry <- function(comparison = FALSE) {
    last_error <- NULL
    for (optimizer in c("bobyqa", "nloptwrap", "Nelder_Mead")) {
      attempt <- tryCatch({
        form <- as.formula(paste0(backtick_name(trait),
          if (comparison) " ~ ID + (1 | Rep)" else " ~ 1 + (1 | ID) + (1 | Rep)"))
        fit <- capture_fit(lme4::lmer(form, scaled, REML = TRUE,
          control = lme4::lmerControl(optimizer = optimizer)))
        check_fit(fit)
        fit
      }, error = function(e) { last_error <<- conditionMessage(e); NULL })
      if (!is.null(attempt)) {
        if (optimizer != "bobyqa") notes <<- c(notes, paste("Converged after optimizer retry:", optimizer))
        return(attempt)
      }
    }
    stop(last_error %||% "Mixed model could not be fitted.")
  }
  tryCatch({
    comparison <- fit_with_retry(comparison = TRUE)
    residual_var <- stats::sigma(comparison)^2 * response_scale^2
    out$CV_model_residual_variance <- residual_var
    if (is.finite(mean_value) && abs(mean_value) > 0.0001) {
      out$CV_pct <- 100 * sqrt(residual_var) / abs(mean_value)
      out$CV_status <- "Estimated"
    } else {
      out$CV_status <- "Zero mean"
      notes <- c(notes, "CV unavailable: trait mean is zero or near zero.")
    }
  }, error = function(e) notes <<- c(notes, paste("CV unavailable:", conditionMessage(e))))
  tryCatch({
    if (nlevels(d$ID) < 2 || nlevels(d$Rep) < 2 || nrow(d) <= nlevels(d$ID)) stop("Insufficient replicated genotype data.")
    fit <- fit_with_retry()
    vg_scaled <- as.numeric(lme4::VarCorr(fit)$ID[1, 1])
    vg <- vg_scaled * response_scale^2
    out$Genotypic_variance <- vg
    out$Error_variance <- stats::sigma(fit)^2 * response_scale^2
    if (!is.finite(vg)) stop("Invalid genotype variance.")
    if (vg_scaled <= 1e-10 * max(1, stats::sigma(fit)^2)) {
      out$Genotypic_variance <- 0
      out$H2_status <- "Zero genetic variance"
      # Cullis' ratio is undefined at vg=0. Report the zero-reliability limit
      # explicitly as a boundary convention, rather than an interior estimate.
      out$Broad_sense_H2 <- 0
      out$Heritability_method <- "Zero-genetic-variance boundary (H2=0 limiting convention)"
      notes <- c(notes, "Genotype variance is at the zero boundary; H2=0 is the zero-reliability limiting convention.")
    } else {
    x <- lme4::getME(fit, "X"); z <- lme4::getME(fit, "Z")
    lambda <- lme4::getME(fit, "Lambda")
    design <- cbind(Matrix::Matrix(x, sparse = TRUE), z %*% lambda)
    p <- ncol(x); q <- ncol(z)
    equations <- Matrix::crossprod(design) + Matrix::Diagonal(p + q, c(rep(0, p), rep(1, q)))
    groups <- names(lme4::getME(fit, "cnms")); offsets <- lme4::getME(fit, "Gp")
    g <- match("ID", groups)
    indices <- seq.int(offsets[g] + 1L, offsets[g + 1L])
    rhs <- rbind(Matrix::Matrix(0, p, length(indices), sparse = TRUE),
      Matrix::t(lambda[indices, , drop = FALSE]))
    pev <- as.matrix(Matrix::crossprod(rhs, Matrix::solve(equations, rhs))) * stats::sigma(fit)^2 * response_scale^2
    n <- nrow(pev)
    mean_pev <- 2 * (n * sum(diag(pev)) - sum(pev)) / (n * (n - 1))
    h2 <- 1 - mean_pev / (2 * vg)
    if (!is.finite(h2) || h2 < -1e-6 || h2 > 1 + 1e-6) stop("Invalid genotype-difference PEV estimate.")
    out$Mean_PEV_difference <- mean_pev
    out$Broad_sense_H2 <- max(0, min(1, h2))
    out$H2_status <- "Estimated"
    }
    if (lme4::isSingular(fit)) notes <- c(notes, "The H2 model has boundary nuisance variance components.")
  }, error = function(e) notes <<- c(notes, paste("H2 not estimable:", conditionMessage(e))))
  out$CV_pct <- round(out$CV_pct, 2)
  out$Broad_sense_H2 <- round(out$Broad_sense_H2, 4)
  out$Note <- paste(unique(c(notes, warnings)), collapse = " | ")
  out
}

si_lpsi_final_diagnostics <- function(data, traits, model_type) {
  dplyr::bind_rows(lapply(traits, function(trait) {
    d <- data[is.finite(data[[trait]]), , drop = FALSE]
    fit <- fit_single_location_model(d, trait, model_type)
    residual <- as.numeric(residuals(fit))
    h <- as.numeric(hatvalues(fit))
    standardized <- residual / (stats::sigma(fit) * sqrt(pmax(1 - h, 0)))
    standardized[!is.finite(standardized)] <- NA_real_
    data.frame(Trait = trait, ID = as.character(d$ID), Original_ID = d$Original_ID,
      Rep = as.character(d$Rep), Observed = d[[trait]], Fitted = as.numeric(fitted(fit)),
      Residual = residual, Standardized_residual = standardized,
      Model = paste(deparse(formula(fit)), collapse = " "), stringsAsFactors = FALSE)
  }))
}

si_lpsi_diagnostic_plot <- function(result, trait) {
  d <- result$residual_diagnostics
  d <- d[d$Trait == trait & is.finite(d$Residual), , drop = FALSE]
  if (!nrow(d)) return(empty_plot("No fitted-model diagnostics are available for this trait."))
  qq <- ggplot2::ggplot(d, ggplot2::aes(sample = Residual)) +
    ggplot2::stat_qq(color = "#3498DB") + ggplot2::stat_qq_line(color = "#E74C3C") +
    ggplot2::labs(title = paste("Final model Q-Q:", trait), x = "Theoretical quantiles", y = "Residual quantiles") + ggplot2::theme_bw()
  rvf <- ggplot2::ggplot(d, ggplot2::aes(Fitted, Residual)) +
    ggplot2::geom_point(color = "#3498DB") + ggplot2::geom_hline(yintercept = 0, color = "#E74C3C", linetype = "dashed") +
    ggplot2::labs(title = "Final model residuals vs fitted") + ggplot2::theme_bw()
  qq + rvf
}

si_ext_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this analysis. Install it before running this extension.",
         call. = FALSE)
  }
  invisible(TRUE)
}

si_ext_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

si_ext_chr <- function(x) {
  trimws(as.character(x))
}

si_ext_safe <- function(expr) {
  tryCatch(expr, error = function(e) {
    structure(list(error = conditionMessage(e)), class = "si_ext_error")
  })
}

si_ext_is_error <- function(x) {
  inherits(x, "si_ext_error")
}

si_ext_pinv <- function(mat, tol = sqrt(.Machine$double.eps)) {
  s <- svd(mat)
  positive <- s$d > tol * max(s$d)
  if (!any(positive)) {
    stop("Covariance matrix is singular and has no positive singular values.", call. = FALSE)
  }
  s$v[, positive, drop = FALSE] %*%
    (t(s$u[, positive, drop = FALSE]) / s$d[positive])
}

si_ext_estimate_trait_values <- function(data, trait, genotype = "GEN", replication = NULL,
                                         method = c("mean", "blue", "blup")) {
  method <- match.arg(method)
  d <- data
  d$.GEN <- factor(si_ext_chr(d[[genotype]]))
  d$.y <- si_ext_num(d[[trait]])
  d <- d[!is.na(d$.GEN) & !is.na(d$.y), , drop = FALSE]
  if (!is.null(replication) && replication %in% names(d)) {
    d$.REP <- factor(si_ext_chr(d[[replication]]))
  }
  if (nrow(d) < 3 || length(unique(d$.GEN)) < 2) {
    return(data.frame(GEN = character(0), Trait = trait, Value = numeric(0), SE = numeric(0)))
  }

  if (method == "mean") {
    out <- stats::aggregate(d$.y, by = list(GEN = d$.GEN), FUN = mean, na.rm = TRUE)
    return(data.frame(GEN = as.character(out$GEN), Trait = trait, Value = out$x, SE = NA_real_))
  }
  if (is.null(replication) || !".REP" %in% names(d)) stop("Select a replication column for BLUE/BLUP estimation; no raw-mean fallback was used.")

  if (method == "blue") {
    fit <- stats::lm(.y ~ .GEN + .REP, data = d)
    if (requireNamespace("emmeans", quietly = TRUE)) {
      em <- as.data.frame(emmeans::emmeans(fit, specs = ".GEN"))
      return(data.frame(GEN = as.character(em$.GEN), Trait = trait, Value = em$emmean, SE = em$SE))
    }
    stop("BLUE estimation requires emmeans; no raw-mean fallback was used.", call. = FALSE)
  }

  si_ext_require("lme4")
  fit <- lme4::lmer(.y ~ .REP + (1 | .GEN), data = d, REML = TRUE)
  re <- lme4::ranef(fit, condVar = TRUE)[[".GEN"]]
  pv <- attr(re, "postVar")
  se <- if (!is.null(pv)) sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1, 1, i], numeric(1))) else NA_real_
  data.frame(
    GEN = rownames(re),
    Trait = trait,
    Value = unname(lme4::fixef(fit)[["(Intercept)"]]) + re[[1]],
    SE = se,
    stringsAsFactors = FALSE
  )
}

si_ext_estimate_genotype_matrix <- function(data, genotype, traits, replication = NULL,
                                            method = c("mean", "blue", "blup")) {
  method <- match.arg(method)
  long <- do.call(rbind, lapply(traits, function(trait) {
    si_ext_estimate_trait_values(data, trait, genotype, replication, method)
  }))
  if (is.null(long) || nrow(long) == 0) {
    stop("No genotype values could be estimated.", call. = FALSE)
  }
  wide <- stats::reshape(
    long[, c("GEN", "Trait", "Value"), drop = FALSE],
    idvar = "GEN",
    timevar = "Trait",
    direction = "wide"
  )
  names(wide) <- sub("^Value\\.", "", names(wide))
  wide[order(wide$GEN), , drop = FALSE]
}

si_ext_residual_covariance <- function(data, genotype, traits, replication = NULL) {
  d <- data
  d$.GEN <- factor(si_ext_chr(d[[genotype]]))
  if (!is.null(replication) && replication %in% names(d)) {
    d$.REP <- factor(si_ext_chr(d[[replication]]))
  }
  for (trait in traits) d[[trait]] <- si_ext_num(d[[trait]])
  d <- d[stats::complete.cases(d[, c(".GEN", traits), drop = FALSE]), , drop = FALSE]
  if (nrow(d) <= length(traits) + 2) {
    return(stats::cov(d[, traits, drop = FALSE], use = "pairwise.complete.obs"))
  }
  dv <- as.matrix(d[, traits, drop = FALSE])
  rhs <- ".GEN"
  if (".REP" %in% names(d) && length(unique(d$.REP)) > 1) rhs <- paste(rhs, "+ .REP")
  fit <- stats::manova(stats::as.formula(paste("dv ~", rhs)), data = d)
  ss <- stats::SSD(fit)
  ss$SSD / ss$df
}

si_ext_mahalanobis_matrix <- function(genotype_values, cov_mat = NULL, genotype_col = "GEN") {
  mat <- as.matrix(genotype_values[, setdiff(names(genotype_values), genotype_col), drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genotype_values[[genotype_col]]
  if (is.null(cov_mat)) cov_mat <- stats::cov(mat, use = "pairwise.complete.obs")
  inv_cov <- si_ext_pinv(cov_mat)
  n <- nrow(mat)
  dm <- matrix(0, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      delta <- mat[i, ] - mat[j, ]
      dm[i, j] <- as.numeric(t(delta) %*% inv_cov %*% delta)
    }
  }
  dm
}

si_ext_cluster_distance_summary <- function(distance_matrix, clusters) {
  cluster_ids <- sort(unique(clusters))
  out <- expand.grid(Cluster1 = cluster_ids, Cluster2 = cluster_ids, stringsAsFactors = FALSE)
  out$Mean_D2 <- vapply(seq_len(nrow(out)), function(i) {
    g1 <- names(clusters)[clusters == out$Cluster1[i]]
    g2 <- names(clusters)[clusters == out$Cluster2[i]]
    vals <- distance_matrix[g1, g2, drop = FALSE]
    if (out$Cluster1[i] == out$Cluster2[i]) vals <- vals[lower.tri(vals)]
    mean(vals, na.rm = TRUE)
  }, numeric(1))
  out
}

si_ext_auto_cluster_k <- function(hc, min_k = 2, max_k = NULL) {
  n <- length(hc$labels)
  if (n <= 2) return(2)
  if (is.null(max_k) || !is.finite(max_k)) {
    max_k <- min(6, n - 1)
  }
  min_k <- max(2, min(as.integer(min_k), n - 1))
  max_k <- max(min_k, min(as.integer(max_k), n - 1))
  heights <- as.numeric(hc$height)
  candidates <- data.frame(k = min_k:max_k)
  candidates$i <- n - candidates$k
  candidates <- candidates[candidates$i >= 1 & candidates$i < length(heights), , drop = FALSE]
  if (nrow(candidates) == 0) return(min_k)
  candidates$gap <- heights[candidates$i + 1] - heights[candidates$i]
  candidates$gap[!is.finite(candidates$gap)] <- -Inf
  candidates$k[which.max(candidates$gap)]
}

si_ext_run_genetic_diversity <- function(data, genotype, traits, replication = NULL,
                                         value_method = c("mean", "blue", "blup"),
                                         n_clusters = NULL, hc_method = "ward.D2",
                                         scale_pca = TRUE) {
  value_method <- match.arg(value_method)
  if (length(traits) < 2) stop("D2 genetic diversity needs at least two traits.", call. = FALSE)

  geno_values <- si_ext_estimate_genotype_matrix(data, genotype, traits, replication, value_method)
  geno_values <- geno_values[stats::complete.cases(geno_values[, traits, drop = FALSE]), , drop = FALSE]
  if (nrow(geno_values) < 2) {
    stop("Diversity analysis needs at least two genotypes with complete trait values.", call. = FALSE)
  }
  cov_mat <- si_ext_safe(si_ext_residual_covariance(data, genotype, traits, replication))
  if (si_ext_is_error(cov_mat)) {
    cov_mat <- stats::cov(geno_values[, traits, drop = FALSE], use = "pairwise.complete.obs")
  }
  d2_matrix <- si_ext_mahalanobis_matrix(geno_values, cov_mat)
  d2_for_clustering <- pmax(d2_matrix, 0)
  diag(d2_for_clustering) <- 0
  d2_dist <- stats::as.dist(sqrt(d2_for_clustering))
  hc <- stats::hclust(d2_dist, method = hc_method)
  if (is.null(n_clusters)) {
    n_clusters <- si_ext_auto_cluster_k(hc, min_k = 2, max_k = min(6, nrow(geno_values) - 1))
  }
  n_clusters <- suppressWarnings(as.integer(n_clusters))
  if (length(n_clusters) == 0 || !is.finite(n_clusters) || n_clusters < 2) {
    n_clusters <- si_ext_auto_cluster_k(hc, min_k = 2, max_k = min(6, nrow(geno_values) - 1))
  }
  n_clusters <- max(2, min(n_clusters, nrow(geno_values)))
  clusters <- stats::cutree(hc, k = n_clusters)
  clusters <- clusters[match(rownames(d2_matrix), names(clusters))]
  cluster_table <- data.frame(GEN = names(clusters), Cluster = as.integer(clusters), row.names = NULL)

  pca <- stats::prcomp(geno_values[, traits, drop = FALSE], center = TRUE, scale. = scale_pca)
  eigenvalues <- data.frame(
    PC = paste0("PC", seq_along(pca$sdev)),
    Eigenvalue = pca$sdev^2,
    Variance_pct = (pca$sdev^2 / sum(pca$sdev^2)) * 100
  )

  cor_mat <- stats::cor(geno_values[, traits, drop = FALSE], use = "pairwise.complete.obs")

  list(
    analysis = "D2 genetic diversity",
    genotype_values = geno_values,
    covariance = cov_mat,
    d2_matrix = d2_matrix,
    distance_matrix = sqrt(d2_for_clustering),
    hclust = hc,
    clusters = cluster_table,
    cluster_distances = si_ext_cluster_distance_summary(d2_matrix, clusters),
    pca = pca,
    pca_eigenvalues = eigenvalues,
    correlation = cor_mat
  )
}

si_ext_run_metan_met <- function(data, env, gen, rep, trait) {
  si_ext_require("metan")
  d <- data
  d$ENV <- factor(si_ext_chr(d[[env]]))
  d$GEN <- factor(si_ext_chr(d[[gen]]))
  d$REP <- factor(si_ext_chr(d[[rep]]))
  d$TRAIT <- si_ext_num(d[[trait]])
  d <- d[stats::complete.cases(d[, c("ENV", "GEN", "REP", "TRAIT")]), , drop = FALSE]
  if (nrow(d) == 0) stop("No complete ENV/GEN/REP/trait rows are available.", call. = FALSE)

  list(
    analysis = "Expanded MET stability via metan",
    trait = trait,
    descriptive = si_ext_safe(metan::desc_stat(d, stats = "all")),
    means_genotype = si_ext_safe(metan::means_by(d, GEN)),
    means_environment = si_ext_safe(metan::means_by(d, ENV)),
    anova_individual = si_ext_safe(metan::anova_ind(d, ENV, GEN, REP, resp = TRAIT)),
    anova_joint = si_ext_safe(metan::anova_joint(d, ENV, GEN, REP, TRAIT)),
    bartlett = si_ext_safe(stats::bartlett.test(d$TRAIT ~ d$ENV)),
    annicchiarico = si_ext_safe(metan::Annicchiarico(d, ENV, GEN, REP, TRAIT)),
    ecovalence = si_ext_safe(metan::ecovalence(d, ENV, GEN, REP, TRAIT)),
    shukla = si_ext_safe(metan::Shukla(d, ENV, GEN, REP, TRAIT)),
    regression_stability = si_ext_safe(metan::ge_reg(d, ENV, GEN, REP, TRAIT)),
    superiority = si_ext_safe(metan::superiority(d, ENV, GEN, TRAIT)),
    fox = si_ext_safe(metan::Fox(d, ENV, GEN, TRAIT)),
    factor_analysis = si_ext_safe(metan::ge_factanal(d, ENV, GEN, REP, TRAIT)),
    stability_wrap = si_ext_safe(metan::ge_stats(d, ENV, GEN, REP, TRAIT)),
    ammi = si_ext_safe(metan::performs_ammi(d, ENV, GEN, REP, TRAIT)),
    waas = si_ext_safe(metan::waas(d, ENV, GEN, REP, TRAIT)),
    gge_environment = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "environment")),
    gge_genotype = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "genotype")),
    gge_symmetrical = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "symmetrical"))
  )
}

si_ext_goal_codes <- function(goals, traits) {
  if (is.null(goals)) goals <- stats::setNames(rep("h", length(traits)), traits)
  goals <- goals[traits]
  goals[is.na(goals) | goals == ""] <- "h"
  goals <- tolower(as.character(goals))
  goals[goals %in% c("high", "higher", "increase", "max", "maximize", "+")] <- "h"
  goals[goals %in% c("low", "lower", "decrease", "min", "minimize", "-")] <- "l"
  goals[!goals %in% c("h", "l")] <- "h"
  stats::setNames(goals, traits)
}

si_ext_make_blup_matrix <- function(model_obj, traits) {
  out <- data.frame(GEN = as.character(model_obj[[traits[1]]]$BLUPgen$GEN), stringsAsFactors = FALSE)
  for (trait in traits) {
    tmp <- data.frame(
      GEN = as.character(model_obj[[trait]]$BLUPgen$GEN),
      value = model_obj[[trait]]$BLUPgen$Predicted,
      stringsAsFactors = FALSE
    )
    names(tmp)[2] <- trait
    out <- merge(out, tmp, by = "GEN", all = TRUE)
  }
  out
}

si_ext_selected_n <- function(total_n, selection_intensity) {
  max(1, round(total_n * selection_intensity / 100, 0))
}

si_ext_selection_differential <- function(blup_mat, selected_geno, goals, method_name) {
  selected <- blup_mat[blup_mat$GEN %in% selected_geno, , drop = FALSE]
  rows <- lapply(names(goals), function(trait) {
    overall <- mean(blup_mat[[trait]], na.rm = TRUE)
    selected_mean <- mean(selected[[trait]], na.rm = TRUE)
    gain <- selected_mean - overall
    data.frame(
      Method = method_name,
      Trait = trait,
      Goal = ifelse(goals[[trait]] == "l", "decrease", "increase"),
      OverallMean = overall,
      SelectedMean = selected_mean,
      SD = gain,
      SDpercent = gain / abs(overall) * 100,
      DesiredGain = ifelse((goals[[trait]] == "h" && gain > 0) ||
                             (goals[[trait]] == "l" && gain < 0), "Yes", "No"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

si_ext_selection_membership <- function(selected_list, all_genotypes) {
  out <- data.frame(GEN = all_genotypes, stringsAsFactors = FALSE)
  for (nm in names(selected_list)) out[[nm]] <- out$GEN %in% selected_list[[nm]]
  method_cols <- setdiff(names(out), "GEN")
  out$MethodsSelected <- rowSums(out[, method_cols, drop = FALSE])
  out$Pattern <- apply(out[, method_cols, drop = FALSE], 1, function(x) {
    hits <- names(x)[as.logical(x)]
    if (length(hits) == 0) "None" else paste(hits, collapse = " + ")
  })
  out[order(-out$MethodsSelected, out$GEN), , drop = FALSE]
}

si_ext_overlap_table <- function(selected_list, total_n) {
  methods <- names(selected_list)
  rows <- expand.grid(Method1 = methods, Method2 = methods, stringsAsFactors = FALSE)
  rows$CoincidenceIndex <- vapply(seq_len(nrow(rows)), function(i) {
    a <- selected_list[[rows$Method1[i]]]
    b <- selected_list[[rows$Method2[i]]]
    if (requireNamespace("metan", quietly = TRUE)) {
      return(tryCatch(
        as.numeric(metan::coincidence_index(sel1 = a, sel2 = b, total = total_n)),
        error = function(e) NA_real_
      ))
    }
    length(intersect(a, b)) / max(1, length(unique(c(a, b))))
  }, numeric(1))
  rows
}

si_ext_make_sh_index <- function(blup_mat, raw_data, traits, goals, selection_intensity,
                                 economic_weights = NULL) {
  si_ext_require("metan")
  pm <- stats::aggregate(raw_data[, traits, drop = FALSE],
                         by = list(GEN = raw_data$GEN), FUN = mean, na.rm = TRUE)
  pcov <- stats::cov(as.matrix(pm[, traits, drop = FALSE]), use = "pairwise.complete.obs")
  gcov <- stats::cov(as.matrix(blup_mat[, traits, drop = FALSE]), use = "pairwise.complete.obs")
  if (is.null(economic_weights)) economic_weights <- stats::setNames(rep(1, length(traits)), traits)
  weights <- ifelse(goals == "l", -1, 1) * as.numeric(economic_weights[traits])
  blup_for_sh <- as.matrix(blup_mat[, traits, drop = FALSE])
  rownames(blup_for_sh) <- blup_mat$GEN
  metan::Smith_Hazel(
    blup_for_sh,
    pcov = pcov,
    gcov = gcov,
    weights = weights,
    SI = selection_intensity
  )
}

si_ext_run_multitrait_selection <- function(data, env, gen, rep, traits, goals = NULL,
                                            selection_intensity = 15, yield_trait = NULL,
                                            economic_weights = NULL) {
  si_ext_require("metan")
  si_ext_require("rlang")
  d <- data
  d$ENV <- factor(si_ext_chr(d[[env]]))
  d$GEN <- factor(si_ext_chr(d[[gen]]))
  d$REP <- factor(si_ext_chr(d[[rep]]))
  for (trait in traits) d[[trait]] <- si_ext_num(d[[trait]])
  if (any(!vapply(d[, traits, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
    stop("Multi-trait selection requires complete trait records for this method. Missing/nonfinite measurements were found; no plots were silently removed. Review the data or select a supported trait set before continuing.", call. = FALSE)
  }
  d <- d[stats::complete.cases(d[, c("ENV", "GEN", "REP", traits), drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0) stop("No complete rows are available for multi-trait selection.", call. = FALSE)

  goals <- si_ext_goal_codes(goals, traits)
  ideotype_words <- ifelse(goals == "l", "min", "max")
  if (is.null(yield_trait) || !yield_trait %in% traits) yield_trait <- traits[length(traits)]

  gamem <- si_ext_safe(metan::gamem_met(d, ENV, GEN, REP, resp = dplyr::all_of(traits), verbose = FALSE))
  waasb <- si_ext_safe(metan::waasb(d, ENV, GEN, REP, resp = dplyr::all_of(traits),
                                    ideotype = paste(ideotype_words, collapse = ","),
                                    verbose = FALSE))
  if (si_ext_is_error(gamem)) stop("gamem_met failed: ", gamem$error, call. = FALSE)

  blup_mat <- si_ext_make_blup_matrix(gamem, traits)
  selected_n <- si_ext_selected_n(nrow(blup_mat), selection_intensity)

  direct <- blup_mat[order(if (goals[[yield_trait]] == "l") blup_mat[[yield_trait]] else -blup_mat[[yield_trait]]),
                     , drop = FALSE]
  direct <- direct[seq_len(min(selected_n, nrow(direct))), , drop = FALSE]

  mtsi <- if (!si_ext_is_error(waasb)) {
    si_ext_safe(metan::mtsi(waasb, SI = selection_intensity))
  } else waasb
  mgidi <- si_ext_safe(metan::mgidi(gamem, ideotype = unname(goals), SI = selection_intensity, verbose = FALSE))
  fai <- si_ext_safe(metan::fai_blup(gamem, DI = ideotype_words, SI = selection_intensity, verbose = FALSE))
  smith_hazel <- si_ext_safe(si_ext_make_sh_index(blup_mat, d, traits, goals, selection_intensity, economic_weights))

  selected_sets <- list(DirectSelection = as.character(direct$GEN))
  if (!si_ext_is_error(mtsi)) selected_sets$MTSI <- as.character(mtsi$sel_gen)
  if (!si_ext_is_error(mgidi)) selected_sets$MGIDI <- as.character(mgidi$sel_gen)
  if (!si_ext_is_error(fai)) selected_sets$FAIBLUP <- tryCatch(as.character(fai$sel_gen$ID1), error = function(e) character(0))
  if (!si_ext_is_error(smith_hazel)) selected_sets$SmithHazel <- tryCatch(as.character(smith_hazel$sel_gen), error = function(e) character(0))

  differentials <- do.call(rbind, lapply(names(selected_sets), function(nm) {
    si_ext_selection_differential(blup_mat, selected_sets[[nm]], goals, nm)
  }))

  list(
    analysis = "Multi-trait selection indices",
    gamem = gamem,
    waasb = waasb,
    blup_matrix = blup_mat,
    direct_selection = direct,
    mtsi = mtsi,
    mgidi = mgidi,
    fai_blup = fai,
    smith_hazel = smith_hazel,
    selected_sets = selected_sets,
    selection_differentials = differentials,
    coincidence = si_ext_overlap_table(selected_sets, nrow(blup_mat)),
    membership = si_ext_selection_membership(selected_sets, blup_mat$GEN),
    gt_biplot = si_ext_safe({
      gm_means <- stats::aggregate(d[, traits, drop = FALSE], by = list(GEN = d$GEN), mean)
      metan::gtb(gm_means, gen = GEN, resp = dplyr::all_of(traits))
    }),
    gyt_biplot = si_ext_safe({
      gm_means <- stats::aggregate(d[, traits, drop = FALSE], by = list(GEN = d$GEN), mean)
      rlang::inject(metan::gytb(
        gm_means,
        gen = GEN,
        yield = !!rlang::sym(yield_trait),
        traits = dplyr::all_of(setdiff(traits, yield_trait)),
        ideotype = unname(goals[setdiff(traits, yield_trait)])
      ))
    })
  )
}

si_ext_available_appR_extensions <- function() {
  data.frame(
    Extension = c(
      "D2 genetic diversity",
      "Expanded MET stability",
      "Multi-trait selection indices"
    ),
    Main_function = c(
      "si_ext_run_genetic_diversity()",
      "si_ext_run_metan_met()",
      "si_ext_run_multitrait_selection()"
    ),
    AppR_source_module = c(
      "Module 1 - D2 Genetic Diversity Analyser",
      "Module 2 - MET Analysis",
      "Module 3 - Multi-Trait Selection Suite"
    ),
    Adds_to_SIR = c(
      "Mahalanobis D2, clustering, PCA, trait correlation",
      "Annicchiarico, ecovalence, Shukla, regression, non-parametric, factor analysis, AMMI/WAAS/GGE via metan",
      "MTSI, MGIDI, FAI-BLUP, Smith-Hazel, direct selection, coincidence, membership, GT/GYT"
    ),
    stringsAsFactors = FALSE
  )
}
