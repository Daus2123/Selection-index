# mating_design_module.R
#
# Standalone statistical functions for mating-design analysis.
# This file intentionally contains no Shiny UI or server code.

# Analysis functions return lists of ANOVA, GCA, SCA, and related result tables.

add_significance_stars_robust <- function(df) {
  pval_col_name <- dplyr::case_when(
    "Pr(>F)" %in% names(df) ~ "Pr(>F)", "p_value" %in% names(df) ~ "p_value",
    "Pr(>|t|)" %in% names(df) ~ "Pr(>|t|)", TRUE ~ NA_character_
  )
  if (is.na(pval_col_name)) return(df)
  p_values <- df[[pval_col_name]]
  stars <- dplyr::case_when(
    is.na(p_values) ~ "", p_values < 0.001 ~ "***", p_values < 0.01 ~ "**",
    p_values < 0.05 ~ "*", p_values < 0.1 ~ ".", TRUE ~ ""
  )
  df$Signif <- stars
  return(df)
}

griffing_method1 <- function(df, rep_col = "Rep", male_col = "Male", female_col = "Female", trait_col = "Trait", blk_col = NULL) {
  data_ab1 <- df
  data_ab1$Rep    <- as.factor(data_ab1[[rep_col]])
  if (!is.null(blk_col)) data_ab1$Blk <- as.factor(data_ab1[[blk_col]])
  data_ab1$Male   <- as.character(data_ab1[[male_col]])
  data_ab1$Female <- as.character(data_ab1[[female_col]])
  data_ab1$YVAR   <- as.numeric(as.character(data_ab1[[trait_col]]))
  data_ab1 <- data_ab1[!is.na(data_ab1$Male) & !is.na(data_ab1$Female) & !is.na(data_ab1$YVAR), ]
  data_ab1$Male <- factor(data_ab1$Male)
  data_ab1$Female <- factor(data_ab1$Female)
  bc <- nlevels(data_ab1$Rep)
  ptypes <- sort(unique(c(as.character(data_ab1$Male), as.character(data_ab1$Female))))
  p <- length(ptypes)
  means_df <- aggregate(YVAR ~ Male + Female, data = data_ab1, mean)
  means_df <- merge(means_df, expand.grid(Male = ptypes, Female = ptypes), all.y = TRUE)
  myMatrix <- matrix(NA, nrow = p, ncol = p, dimnames = list(ptypes, ptypes))
  for (i in 1:nrow(means_df)) {
    myMatrix[as.character(means_df$Male[i]), as.character(means_df$Female[i])] <- means_df$YVAR[i]
  }
  if (any(is.na(myMatrix))) stop("Missing values in means matrix for Griffing I. Ensure all p^2 crosses are present in the data.")
  modelg1 <- lm(YVAR ~ factor(Rep) + factor(paste(Male, Female, sep="_x_")), data = data_ab1)
  anmodel <- anova(modelg1); rownames(anmodel)[nrow(anmodel)] <- "Residual"
  MSEAD <- as.numeric(anmodel[nrow(anmodel), "Mean Sq"]); error_DF <- as.numeric(anmodel[nrow(anmodel), "Df"])
  Xi.    <- rowSums(myMatrix)
  X.j    <- colSums(myMatrix)
  Xbar   <- sum(myMatrix)
  acon   <- sum((Xi. + X.j)^2) / (2*p)
  ssgca <- acon - (2/(p^2))*(Xbar^2)
  sssca <- sum(myMatrix * (myMatrix + t(myMatrix)))/2 - acon + (Xbar^2)/(p^2)
  ssrecp <- sum((myMatrix - t(myMatrix))^2)/4
  ssmat <- sum((Xi. - X.j)^2)/(2*p)
  ssnomat <- ssrecp - ssmat
  
  Df_display <- c(p-1, p*(p-1)/2, p*(p-1)/2)
  SSS_display <- c(ssgca, sssca, ssrecp) * bc
  MSSS_display <- SSS_display / Df_display
  names(SSS_display) <- names(MSSS_display) <- c("GCA", "SCA", "Reciprocal")
  FVAL_display <- MSSS_display / MSEAD
  pval_display <- 1 - pf(FVAL_display, Df_display, error_DF)
  
  anova_diallel <- data.frame(Df=Df_display, `Sum Sq`=SSS_display, `Mean Sq`=MSSS_display, `F value`=FVAL_display, `Pr(>F)`=pval_display, row.names=names(SSS_display), check.names=FALSE)
  anova_diallel <- add_significance_stars_robust(anova_diallel)
  anova_error <- data.frame(Df=error_DF, `Sum Sq`=MSEAD*error_DF, `Mean Sq`=MSEAD, `F value`=NA, `Pr(>F)`=NA, Signif="", row.names="Residual", check.names=FALSE)
  anova_final <- rbind(anova_diallel, anova_error)
  
  gca <- (Xi. + X.j) / (2*p) - Xbar/(p^2)
  gca_se <- sqrt(((p-1) * MSEAD) / (2*p*p*bc))
  parental_means <- diag(myMatrix)
  gca_df <- data.frame(Parent=ptypes, Parental_Mean = parental_means, GCA=gca, SE=gca_se, T_value = gca/gca_se)
  gca_df$p_value <- 2 * pt(-abs(gca_df$T_value), df = Df_display[1])
  gca_df <- add_significance_stars_robust(gca_df)
  
  sca <- (myMatrix + t(myMatrix))/2 - (matrix(Xi. + X.j, nrow=p, ncol=p, byrow=TRUE) + matrix(Xi. + X.j, nrow=p, ncol=p, byrow=FALSE))/(2*p) + Xbar/(p^2)
  sca[lower.tri(sca)] <- NA
  
  # --- FIX: Define sca_se before using it ---
  sca_se <- sqrt(((p-1) * MSEAD) / (2 * bc))
  
  sca_df <- data.frame(expand.grid(Female=ptypes, Male=ptypes), SCA=as.vector(sca))
  sca_df <- sca_df[!is.na(sca_df$SCA),]
  cross_means_vec <- as.vector(myMatrix)
  sca_df$Cross_Mean <- cross_means_vec[!is.na(as.vector(sca))]
  sca_df$SE <- sca_se
  sca_df$T_value <- sca_df$SCA / sca_df$SE
  sca_df$p_value <- 2 * pt(-abs(sca_df$T_value), df = Df_display[2])
  sca_df <- add_significance_stars_robust(sca_df)
  
  return(list(method="I", anova=anova_final, gca=gca_df, sca=sca_df))
}


griffing_method2 <- function(df, rep_col = "Rep", male_col = "Male", female_col = "Female", trait_col = "Trait", blk_col = NULL) {
  data_ab1 <- df
  data_ab1$Rep    <- as.factor(data_ab1[[rep_col]])
  if (!is.null(blk_col)) data_ab1$Blk <- as.factor(data_ab1[[blk_col]])
  data_ab1$Male   <- as.character(data_ab1[[male_col]])
  data_ab1$Female <- as.character(data_ab1[[female_col]])
  data_ab1$YVAR   <- as.numeric(as.character(data_ab1[[trait_col]]))
  data_ab1 <- data_ab1[!is.na(data_ab1$Male) & !is.na(data_ab1$Female) & !is.na(data_ab1$YVAR), ]
  data_ab1$Male <- factor(data_ab1$Male)
  data_ab1$Female <- factor(data_ab1$Female)
  bc <- nlevels(data_ab1$Rep)
  ptypes <- sort(unique(c(as.character(data_ab1$Male), as.character(data_ab1$Female))))
  p <- length(ptypes)
  data_ab1$Cross <- factor(paste(pmin(data_ab1$Female, data_ab1$Male), pmax(data_ab1$Female, data_ab1$Male), sep = "_x_"))
  means_df <- aggregate(YVAR ~ Male + Female, data = data_ab1, mean)
  means_df <- merge(means_df, expand.grid(Male = ptypes, Female = ptypes), all.y = TRUE)
  myMatrix <- matrix(NA, nrow = p, ncol = p, dimnames = list(ptypes, ptypes))
  for (i in 1:nrow(means_df)) {
    myMatrix[as.character(means_df$Male[i]), as.character(means_df$Female[i])] <- means_df$YVAR[i]
  }
  myMatrix[lower.tri(myMatrix, diag=FALSE)] <- t(myMatrix)[lower.tri(t(myMatrix), diag=FALSE)] # Symmetrize
  if (any(is.na(myMatrix[upper.tri(myMatrix, diag=TRUE)]))) stop("Missing values in means matrix for Griffing II. Ensure p(p+1)/2 crosses are present.")
  modelg <- lm(YVAR ~ Cross + Rep, data = data_ab1)
  anmodel <- anova(modelg)
  rownames(anmodel)[nrow(anmodel)] <- "Residual"
  MSEAD <- as.numeric(anmodel["Residual", "Mean Sq"])
  error_DF <- as.numeric(anmodel["Residual", "Df"])
  Xi. <- rowSums(myMatrix)
  Xbar <- sum(myMatrix[upper.tri(myMatrix, diag=TRUE)])
  acon <- sum((Xi. + diag(myMatrix))^2) / (p + 2)
  ssgca <- acon - (4 * Xbar^2) / (p * (p + 2))
  sssca <- sum((myMatrix[upper.tri(myMatrix, diag=TRUE)])^2) - acon + (2 * Xbar^2) / ((p + 1) * (p + 2))
  Df <- c(p - 1, p * (p - 1) / 2)
  SSS <- c(ssgca, sssca) * bc
  MSSS <- SSS / Df
  FVAL <- MSSS / MSEAD
  pval <- 1 - pf(FVAL, Df, error_DF)
  anova_diallel <- data.frame(Df=Df, `Sum Sq`=SSS, `Mean Sq`=MSSS, `F value`=FVAL, `Pr(>F)`=pval, row.names=c("GCA", "SCA"), check.names=FALSE)
  anova_diallel <- add_significance_stars_robust(anova_diallel)
  anova_error <- data.frame(Df=error_DF, `Sum Sq`=MSEAD*error_DF, `Mean Sq`=MSEAD, `F value`=NA, `Pr(>F)`=NA, Signif="", row.names="Residual", check.names=FALSE)
  anova_final <- rbind(anova_diallel, anova_error)
  gca <- (Xi. + diag(myMatrix) - 2 * Xbar / p) / (p + 2)
  gca_se <- sqrt(((p - 1) * MSEAD) / (p * (p + 2) * bc))
  parental_means <- diag(myMatrix)
  gca_df <- data.frame(Parent=ptypes, Parental_Mean = parental_means, GCA=gca, SE=gca_se, T_value=gca/gca_se)
  gca_df$p_value <- 2 * pt(-abs(gca_df$T_value), df = Df[1])
  gca_df <- add_significance_stars_robust(gca_df)
  sca_mat <- myMatrix - (matrix(Xi.+diag(myMatrix),nrow=p,ncol=p,byrow=TRUE) + matrix(Xi.+diag(myMatrix),nrow=p,ncol=p,byrow=FALSE))/(p+2) + 2*Xbar/((p+1)*(p+2))
  sca_mat[lower.tri(sca_mat)] <- NA
  sca_df <- expand.grid(Male=ptypes, Female=ptypes)
  sca_df$SCA <- as.vector(sca_mat)
  sca_df <- sca_df[!is.na(sca_df$SCA), ]
  cross_means_vec <- as.vector(myMatrix)
  sca_df$Cross_Mean <- cross_means_vec[!is.na(as.vector(sca_mat))]
  sca_se <- sqrt(((p*p + p + 2) * MSEAD) / ((p + 1) * (p + 2) * bc))
  sca_df$SE <- sca_se
  sca_df$T_value <- sca_df$SCA / sca_df$SE
  sca_df$p_value <- 2 * pt(-abs(sca_df$T_value), df = Df[2])
  sca_df <- add_significance_stars_robust(sca_df)
  return(list(method="II", anova=anova_final, gca=gca_df, sca=sca_df, griffing_anova = anova(modelg)))
}

griffing_method3 <- function(df, rep_col = "Rep", male_col = "Male", female_col = "Female", trait_col = "Trait", blk_col = NULL) {
  data_ab1 <- df
  data_ab1$Rep    <- as.factor(data_ab1[[rep_col]])
  if (!is.null(blk_col)) data_ab1$Blk <- as.factor(data_ab1[[blk_col]])
  data_ab1$Male   <- as.character(data_ab1[[male_col]])
  data_ab1$Female <- as.character(data_ab1[[female_col]])
  data_ab1$YVAR   <- as.numeric(as.character(data_ab1[[trait_col]]))
  data_ab1 <- data_ab1[!is.na(data_ab1$Male) & !is.na(data_ab1$Female) & !is.na(data_ab1$YVAR), ]
  data_ab1 <- data_ab1[data_ab1$Male != data_ab1$Female, ] # Method 3 excludes selfs
  data_ab1$Male <- factor(data_ab1$Male)
  data_ab1$Female <- factor(data_ab1$Female)
  bc <- nlevels(data_ab1$Rep)
  ptypes <- sort(unique(c(as.character(data_ab1$Male), as.character(data_ab1$Female))))
  p <- length(ptypes)
  means_df <- aggregate(YVAR ~ Male + Female, data = data_ab1, mean)
  means_df <- merge(means_df, expand.grid(Male = ptypes, Female = ptypes), all.y = TRUE)
  myMatrix <- matrix(NA, nrow = p, ncol = p, dimnames = list(ptypes, ptypes))
  for (i in 1:nrow(means_df)) {
    myMatrix[means_df$Male[i], means_df$Female[i]] <- means_df$YVAR[i]
  }
  diag(myMatrix) <- NA
  if (any(is.na(myMatrix[upper.tri(myMatrix)]) | is.na(myMatrix[lower.tri(myMatrix)]))) {
    stop("Griffing Method III: Incomplete matrix; missing F1 or reciprocal data (excluding selfs).")
  }
  modelg <- lm(YVAR ~ factor(Rep) + factor(paste(Male, Female, sep = "_x_")), data = data_ab1)
  anmodel <- anova(modelg); rownames(anmodel)[nrow(anmodel)] <- "Residual"
  MSEAD <- as.numeric(anmodel[nrow(anmodel), "Mean Sq"]); error_DF <- as.numeric(anmodel[nrow(anmodel), "Df"])
  Xi.    <- rowSums(myMatrix, na.rm = TRUE)
  X.j    <- colSums(myMatrix, na.rm = TRUE)
  Xbar   <- sum(myMatrix, na.rm = TRUE)
  acon <- sum((Xi. + X.j)^2) / (2 * (p - 2))
  ssgca <- acon - (2 / (p * (p - 2))) * (Xbar^2)
  sssca <- sum((myMatrix + t(myMatrix))^2, na.rm=TRUE)/4 - acon + (Xbar^2)/((p-1)*(p-2))
  ssrecp <- sum((myMatrix - t(myMatrix))^2, na.rm=TRUE)/4
  
  Df_display <- c((p - 1), (p * (p - 3) / 2), (p * (p - 1) / 2))
  SSS_display <- c(ssgca, sssca, ssrecp) * bc
  MSSS_display <- SSS_display / Df_display
  names(SSS_display) <- names(MSSS_display) <- c("GCA", "SCA", "Reciprocal")
  FVAL_display <- MSSS_display / MSEAD
  pval_display <- 1 - pf(FVAL_display, Df_display, error_DF)
  
  anova_diallel <- data.frame(Df=Df_display, `Sum Sq`=SSS_display, `Mean Sq`=MSSS_display, `F value`=FVAL_display, `Pr(>F)`=pval_display, row.names=names(SSS_display), check.names=FALSE)
  anova_diallel <- add_significance_stars_robust(anova_diallel)
  anova_error <- data.frame(Df=error_DF, `Sum Sq`=MSEAD*error_DF, `Mean Sq`=MSEAD, `F value`=NA, `Pr(>F)`=NA, Signif="", row.names="Residual", check.names=FALSE)
  anova_final <- rbind(anova_diallel, anova_error)
  
  gca <- (p * (Xi. + X.j) - 2 * Xbar) / (2 * p * (p - 2))
  gca_se <- sqrt(((p - 1) * MSEAD) / (2 * p * (p - 2) * bc))
  gca_df <- data.frame(Parent=ptypes, GCA=gca, SE=gca_se, T_value=gca/gca_se)
  gca_df$p_value <- 2 * pt(-abs(gca_df$T_value), df = Df_display[1])
  gca_df <- add_significance_stars_robust(gca_df)
  
  sca <- (myMatrix + t(myMatrix))/2 - (matrix(Xi.+X.j,nrow=p,ncol=p,byrow=TRUE) + matrix(Xi.+X.j,nrow=p,ncol=p,byrow=FALSE))/(2*(p-2)) + Xbar/((p-1)*(p-2))
  sca[lower.tri(sca, diag = TRUE)] <- NA
  sca_df <- expand.grid(Female=ptypes, Male=ptypes)
  sca_df$SCA <- as.vector(sca)
  sca_df <- sca_df[!is.na(sca_df$SCA), ]
  cross_means_vec <- as.vector(myMatrix)
  sca_df$Cross_Mean <- cross_means_vec[!is.na(as.vector(sca))]
  sca_se <- sqrt(((p - 3) * MSEAD) / (2 * (p - 1) * bc))
  sca_df$SE <- sca_se
  sca_df$T_value <- sca_df$SCA / sca_df$SE
  sca_df$p_value <- 2 * pt(-abs(sca_df$T_value), df = Df_display[2])
  sca_df <- add_significance_stars_robust(sca_df)
  
  return(list(method="III", anova=anova_final, gca=gca_df, sca=sca_df))
}


griffing_method4 <- function(df, rep_col, male_col, female_col, trait_col, blk_col = NULL) {
  dat <- df[, c(rep_col, male_col, female_col, trait_col)]
  names(dat) <- c("Rep", "Male", "Female", "YVAR")
  dat$YVAR <- as.numeric(as.character(dat$YVAR))
  dat <- dat[!is.na(dat$Male) & !is.na(dat$Female) & !is.na(dat$YVAR), ]
  dat <- dat[dat$Male != dat$Female, ] # Method 4 excludes selfs
  dat$Rep <- factor(dat$Rep)
  dat$Cross <- factor(paste(pmin(dat$Female, dat$Male), pmax(dat$Female, dat$Male), sep = "_x_"))
  bc <- nlevels(dat$Rep)
  ptypes <- sort(unique(c(dat$Male, dat$Female)))
  p <- length(ptypes)
  means_df <- aggregate(YVAR ~ Male + Female, data = dat, mean)
  myMatrix <- matrix(NA, nrow = p, ncol = p, dimnames = list(ptypes, ptypes))
  for (i in 1:nrow(means_df)) {
    myMatrix[as.character(means_df$Male[i]), as.character(means_df$Female[i])] <- means_df$YVAR[i]
  }
  myMatrix[lower.tri(myMatrix)] <- t(myMatrix)[lower.tri(t(myMatrix))]
  if (any(is.na(myMatrix[upper.tri(myMatrix)]))) {
    stop("Incomplete data: One or more cross combinations are missing for Method 4.")
  }
  model_g4 <- lm(YVAR ~ Cross + Rep, data = dat)
  anova_g4 <- anova(model_g4)
  MSEAD <- anova_g4["Residuals", "Mean Sq"]
  Randoms_Df <- anova_g4["Residuals", "Df"]
  diag(myMatrix) <- NA
  Xi. <- rowSums(myMatrix, na.rm = TRUE)
  Xbar <- sum(myMatrix, na.rm = TRUE) / 2
  ssgca <- (1/(p-2)) * sum(Xi.^2) - (4*Xbar^2)/(p*(p-2))
  sssca <- sum(myMatrix[upper.tri(myMatrix)]^2) - (1/(p-2))*sum(Xi.^2) + (2*Xbar^2)/((p-1)*(p-2))
  Df <- c(p - 1, p * (p - 3) / 2)
  SSS <- c(ssgca, sssca) * bc
  MSSS <- SSS / Df
  FVAL <- MSSS / MSEAD
  pval <- 1 - pf(FVAL, Df, Randoms_Df)
  anova_diallel <- data.frame(Df=Df, `Sum Sq`=SSS, `Mean Sq`=MSSS, `F value`=FVAL, `Pr(>F)`=pval, row.names=c("GCA", "SCA"), check.names=FALSE)
  anova_diallel <- add_significance_stars_robust(anova_diallel)
  anova_error <- data.frame(Df=Randoms_Df, `Sum Sq`=MSEAD*Randoms_Df, `Mean Sq`=MSEAD, `F value`=NA, `Pr(>F)`=NA, Signif="", row.names="Residual", check.names=FALSE)
  anova_final <- rbind(anova_diallel, anova_error)
  gcaeff <- (1/(p-2)) * (Xi. - (2 * Xbar / p))
  gca_se <- sqrt(((p-1)*MSEAD)/(bc*p*(p-2)))
  gca_tab <- data.frame(Parent=ptypes, GCA=gcaeff, SE=gca_se)
  gca_tab$T_value <- gca_tab$GCA / gca_tab$SE
  gca_tab$p_value <- 2 * pt(-abs(gca_tab$T_value), df = Df[1])
  gca_tab <- add_significance_stars_robust(gca_tab)
  mu_hat <- (2*Xbar)/(p*(p-1))
  gca_m1 <- matrix(gcaeff, nrow=p, ncol=p, byrow=TRUE)
  gca_m2 <- matrix(gcaeff, nrow=p, ncol=p, byrow=FALSE)
  scaeffmat <- myMatrix - gca_m1 - gca_m2 - mu_hat
  scaeffmat[lower.tri(scaeffmat, diag=TRUE)] <- NA
  sca_tab <- data.frame(expand.grid(Female=ptypes, Male=ptypes), SCA=as.vector(scaeffmat))
  sca_tab <- sca_tab[!is.na(sca_tab$SCA), ]
  cross_means_vec <- as.vector(myMatrix)
  sca_tab$Cross_Mean <- cross_means_vec[!is.na(as.vector(scaeffmat))]
  sca_tab$SE <- sqrt(((p-3)*MSEAD)/(bc*(p-1)))
  sca_tab$T_value <- sca_tab$SCA / sca_tab$SE
  sca_tab$p_value <- 2 * pt(-abs(sca_tab$T_value), df = Df[2])
  sca_tab <- add_significance_stars_robust(sca_tab)
  return(list(method="IV", anova=anova_final, gca=gca_tab, sca=sca_tab))
}

diallel_partial_manual <- function(df, trait, p1, p2, rep) {
  Y <- df[[trait]]
  P1 <- as.character(df[[p1]])
  P2 <- as.character(df[[p2]])
  Rep <- as.factor(df[[rep]])
  mask <- !is.na(Y)
  Y <- Y[mask]; P1 <- P1[mask]; P2 <- P2[mask]; Rep <- Rep[mask]
  parents <- sort(unique(c(P1, P2)))
  n_par <- length(parents)
  mu <- mean(Y, na.rm = TRUE)
  n_obs <- length(Y)
  s_tab <- table(c(P1, P2))
  s_vals <- as.numeric(s_tab)
  s <- mean(s_vals)
  is_balanced <- all(s_vals == s)
  if (!is_balanced) warning(
    paste0("[Partial Diallel] Unbalanced design detected: crosses per parent are not equal (s values: ",
           paste(unique(s_vals), collapse = ","), "). Results use average s = ", s, "."))
  cross_ids <- paste(pmin(P1, P2), pmax(P1, P2), sep = ":")
  reps_per_cross <- as.numeric(median(table(cross_ids)))
  if (any(table(cross_ids) != reps_per_cross)) {
    warning("[Partial Diallel] Unbalanced replication detected across crosses! Using modal # of reps.")
  }
  Xg <- matrix(0, nrow = n_obs, ncol = n_par)
  for (i in 1:n_obs) {
    Xg[i, which(parents == P1[i])] <- 1
    Xg[i, which(parents == P2[i])] <- 1
  }
  Xg_red <- Xg[, -n_par, drop = FALSE]
  Y_c <- Y - mu
  gca_hat_red <- solve(t(Xg_red) %*% Xg_red, t(Xg_red) %*% Y_c)
  gca_hat <- c(gca_hat_red, -sum(gca_hat_red))
  names(gca_hat) <- parents
  cross_names_unique <- unique(cross_ids)
  sca_cross <- numeric(length(cross_names_unique))
  names(sca_cross) <- cross_names_unique
  for (cr_name in cross_names_unique) {
    rows <- which(cross_ids == cr_name)
    ptemp <- unlist(strsplit(cr_name, ":"))
    expected <- mu + gca_hat[ptemp[1]] + gca_hat[ptemp[2]]
    sca_cross[cr_name] <- mean(Y[rows]) - expected
  }
  sca_fitted_for_each_obs <- sca_cross[cross_ids]
  rep_means <- tapply(Y, Rep, mean)
  rep_effects <- rep_means - mu
  rep_effect_for_each_obs <- rep_effects[as.character(Rep)]
  fitted_values <- mu + rep_effect_for_each_obs + gca_hat[P1] + gca_hat[P2] + sca_fitted_for_each_obs
  residuals <- Y - fitted_values
  ss_rep <- sum(table(Rep) * (rep_means - mu)^2)
  df_rep <- length(unique(Rep)) - 1
  ss_gca <- 2 * s * sum(gca_hat^2)
  df_gca <- n_par - 1
  ss_sca <- sum(sca_cross^2) * reps_per_cross
  df_sca <- length(unique(cross_ids)) - n_par
  ss_err <- sum(residuals^2)
  df_err <- n_obs - (df_rep + df_gca + df_sca + 1)
  ss_total <- sum((Y - mu)^2)
  df_total <- n_obs - 1
  ms_rep <- ss_rep / df_rep
  ms_gca <- ss_gca / df_gca
  ms_sca <- ss_sca / df_sca
  ms_err <- ss_err / df_err
  f_gca <- ms_gca / ms_err
  f_sca <- ms_sca / ms_err
  f_rep <- ms_rep / ms_err
  p_gca <- pf(f_gca, df_gca, df_err, lower.tail = FALSE)
  p_sca <- pf(f_sca, df_sca, df_err, lower.tail = FALSE)
  p_rep <- pf(f_rep, df_rep, df_err, lower.tail = FALSE)
  stars_func <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.1) "." else ""
  }
  anova_tab <- data.frame(
    Source = c("Replication", "GCA", "SCA", "Error", "Total"),
    Df = c(df_rep, df_gca, df_sca, df_err, df_total),
    `Sum Sq` = c(ss_rep, ss_gca, ss_sca, ss_err, ss_total),
    `Mean Sq` = c(ms_rep, ms_gca, ms_sca, ms_err, NA),
    `F value` = c(f_rep, f_gca, f_sca, NA, NA),
    `Pr(>F)` = c(p_rep, p_gca, p_sca, NA, NA),
    Signif = sapply(c(p_rep, p_gca, p_sca, NA, NA), stars_func),
    check.names = FALSE
  )
  r <- reps_per_cross
  n <- n_par
  sigma2_sca <- (ms_sca - ms_err) / r
  sigma2_gca <- ((ms_gca - ms_sca) * (n - 1)) / (r * s * (n - 2))
  var_tab <- data.frame(
    Variance_Component = c("GCA", "SCA"),
    Value = c(sigma2_gca, sigma2_sca)
  )
  gca_tab <- data.frame(Parent = parents, GCA = gca_hat)
  sca_tab <- data.frame(Cross = names(sca_cross), SCA = sca_cross)
  list( anova = anova_tab, gca = gca_tab, sca = sca_tab, var = var_tab)
}

line_tester_manual <- function(data, line_col, tester_col, rep_col, trait_col, type_col) {
  get_stars <- function(p_values) {
    dplyr::case_when(
      is.na(p_values)   ~ "", p_values < 0.001  ~ "***", p_values < 0.01   ~ "**",
      p_values < 0.05   ~ "*", p_values < 0.1    ~ ".", TRUE              ~ ""
    )
  }
  tryCatch({
    df <- data.frame(
      Rep    = as.factor(data[[rep_col]]), Line   = as.factor(data[[line_col]]),
      Tester = as.factor(data[[tester_col]]), Type   = as.factor(data[[type_col]]),
      Y      = as.numeric(data[[trait_col]])
    )
    df <- na.omit(df)
    if (!any(df$Type == "cross")) {
      stop("The 'Type' column must contain entries with the exact value 'cross'.")
    }
    parents <- df[df$Type != "cross", ]
    crosses <- df[df$Type == "cross", ]
    crosses$Line   <- droplevels(crosses$Line)
    crosses$Tester <- droplevels(crosses$Tester)
    has_parents <- nrow(parents) > 0
    l <- nlevels(crosses$Line)
    t <- nlevels(crosses$Tester)
    if (l < 2 || t < 1) {
      stop("Analysis requires at least two lines and one tester in the cross data.")
    }
    df$TreatmentID <- interaction(df$Line, df$Tester, df$Type, drop = TRUE)
    model_overall <- aov(Y ~ Rep + TreatmentID, data = df)
    anova_overall <- anova(model_overall)
    MS_Error <- anova_overall["Residuals", "Mean Sq"]
    DF_Error <- anova_overall["Residuals", "Df"]
    model_crosses <- aov(Y ~ Line * Tester, data = crosses)
    anova_crosses <- anova(model_crosses)
    DF_Parents <- 0; SS_Parents <- 0
    if (has_parents && nlevels(droplevels(parents$Line)) > 1) {
      parents$Parent <- droplevels(parents$Line)
      model_parents <- aov(Y ~ Parent, data = parents)
      anova_parents <- anova(model_parents)
      DF_Parents <- anova_parents["Parent", "Df"]
      SS_Parents <- anova_parents["Parent", "Sum Sq"]
    }
    DF_PvC <- 0; SS_PvC <- 0
    if (has_parents) {
      df$PvC <- factor(ifelse(df$Type == "cross", "Cross", "Parent"), levels = c("Parent", "Cross"))
      model_pvc <- aov(Y ~ PvC, data = df)
      anova_pvc <- anova(model_pvc)
      DF_PvC <- anova_pvc["PvC", "Df"]
      SS_PvC <- anova_pvc["PvC", "Sum Sq"]
    }
    ss_crosses_total <- sum(anova_crosses[c("Line", "Tester", "Line:Tester"), "Sum Sq"])
    df_crosses_total <- sum(anova_crosses[c("Line", "Tester", "Line:Tester"), "Df"])
    source_names <- c("Replications", "Treatments", "  Parents", "  Parents vs. Crosses", "  Crosses",
                      "    Lines", "    Testers", "    Lines X Testers", "Error", "Total")
    DF <- c(anova_overall["Rep", "Df"], anova_overall["TreatmentID", "Df"], DF_Parents, DF_PvC,
            df_crosses_total, anova_crosses["Line", "Df"], anova_crosses["Tester", "Df"],
            anova_crosses["Line:Tester", "Df"], DF_Error, sum(anova_overall[, "Df"], na.rm = TRUE))
    SS <- c(anova_overall["Rep", "Sum Sq"], anova_overall["TreatmentID", "Sum Sq"], SS_Parents, SS_PvC,
            ss_crosses_total, anova_crosses["Line", "Sum Sq"], anova_crosses["Tester", "Sum Sq"],
            anova_crosses["Line:Tester", "Sum Sq"], anova_overall["Residuals", "Sum Sq"],
            sum(anova_overall[, "Sum Sq"], na.rm = TRUE))
    MS <- ifelse(DF > 0, SS / DF, 0)
    F_value <- MS / MS_Error
    P_value <- pf(F_value, DF, DF_Error, lower.tail = FALSE)
    anova_final <- data.frame(Source = source_names, Df = DF, `Sum Sq` = SS, `Mean Sq` = MS,
                              `F value` = F_value, `Pr(>F)` = P_value, check.names = FALSE)
    anova_final[anova_final$Source == "Error", "Mean Sq"] <- MS_Error
    rows_to_blank <- c("Total")
    cols_to_blank <- c("Mean Sq", "F value", "Pr(>F)")
    anova_final[anova_final$Source %in% rows_to_blank, cols_to_blank] <- NA
    anova_final$Signif <- get_stars(anova_final$`Pr(>F)`)
    numeric_cols <- c("Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    anova_final[numeric_cols] <- lapply(anova_final[numeric_cols], function(x) sprintf("%.2f", x))
    anova_final[is.na(anova_final) | anova_final == "NA"] <- ""
    grand_mean <- mean(crosses$Y)
    
    parental_means_lines <- aggregate(Y ~ Line, data = parents, FUN = mean)
    names(parental_means_lines)[2] <- "Parental_Mean"
    
    parental_means_testers <- aggregate(Y ~ Tester, data = parents, FUN = mean)
    names(parental_means_testers)[2] <- "Parental_Mean"
    
    emm_lines <- emmeans::emmeans(model_crosses, ~ Line)
    summary_lines <- as.data.frame(summary(emm_lines))
    gca_lines_out <- data.frame(Line = summary_lines$Line, GCA = summary_lines$emmean - grand_mean, SE = summary_lines$SE)
    gca_lines_out <- merge(gca_lines_out, parental_means_lines, by = "Line", all.x = TRUE)
    gca_lines_out$`t value` <- gca_lines_out$GCA / gca_lines_out$SE
    gca_lines_out$`Pr(>|t|)` <- 2 * pt(-abs(gca_lines_out$`t value`), df = DF_Error)
    gca_lines_out$Signif <- get_stars(gca_lines_out$`Pr(>|t|)`)
    
    emm_testers <- emmeans::emmeans(model_crosses, ~ Tester)
    summary_testers <- as.data.frame(summary(emm_testers))
    gca_testers_out <- data.frame(Tester = summary_testers$Tester, GCA = summary_testers$emmean - grand_mean, SE = summary_testers$SE)
    gca_testers_out <- merge(gca_testers_out, parental_means_testers, by = "Tester", all.x = TRUE)
    gca_testers_out$`t value` <- gca_testers_out$GCA / gca_testers_out$SE
    gca_testers_out$`Pr(>|t|)` <- 2 * pt(-abs(gca_testers_out$`t value`), df = DF_Error)
    gca_testers_out$Signif <- get_stars(gca_testers_out$`Pr(>|t|)`)
    
    cross_means <- aggregate(Y ~ Line + Tester, data = crosses, FUN = mean)
    names(cross_means)[3] <- "Cross_Mean"
    
    emm_sca <- emmeans::emmeans(model_crosses, ~ Line:Tester)
    summary_sca <- as.data.frame(summary(emm_sca))
    summary_sca <- merge(summary_sca, gca_lines_out[, c("Line", "GCA")], by = "Line")
    names(summary_sca)[names(summary_sca) == "GCA"] <- "GCA_Line"
    summary_sca <- merge(summary_sca, gca_testers_out[, c("Tester", "GCA")], by = "Tester")
    names(summary_sca)[names(summary_sca) == "GCA"] <- "GCA_Tester"
    summary_sca$SCA <- summary_sca$emmean - summary_sca$GCA_Line - summary_sca$GCA_Tester - grand_mean
    sca_out <- data.frame(Line = summary_sca$Line, Tester = summary_sca$Tester, SCA = summary_sca$SCA, SE = summary_sca$SE)
    sca_out <- merge(sca_out, cross_means, by = c("Line", "Tester"))
    sca_out$`t value` <- sca_out$SCA / sca_out$SE
    sca_out$`Pr(>|t|)` <- 2 * pt(-abs(sca_out$`t value`), df = DF_Error)
    sca_out$Signif <- get_stars(sca_out$`Pr(>|t|)`)
    list(
      anova_full  = anova_final, gca_lines   = gca_lines_out,
      gca_testers = gca_testers_out, sca         = sca_out
    )
  }, error = function(e) {
    list(error = paste("Line x Tester analysis failed:", e$message))
  })
}

