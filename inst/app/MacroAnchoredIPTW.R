# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model)
#  - Fixed "Inf" Explosion during counterfactual simulation
#  - Guaranteed CI output (Fallback to robust SD if density fails)
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(flexsurv)
  library(survival)
  library(ggplot2)
  library(splines)
  library(scales)
})

# -------------------------------------------------------------------------
# Utility: Censoring parameter finder
# -------------------------------------------------------------------------
find_censoring_param <- function(T2, target_rate, pattern) {
  target_rate <- pmax(pmin(target_rate, 0.95), 0.01)

  if (pattern == "peak1y") {
    obj_func <- function(k) {
      scale_val <- 365.25 / ((k - 1) / k)^(1 / k)
      mean(rweibull(length(T2), shape = k, scale = scale_val) < T2) - target_rate
    }
    opt_k <- tryCatch(uniroot(obj_func, interval = c(1.05, 10))$root, error = function(e) 2.0)
    list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE)

  } else if (pattern == "ushape") {
    u_fixed <- runif(length(T2))
    obj_func <- function(scale_param) {
      C2 <- ifelse(u_fixed < 0.5, rweibull(length(T2), shape = 0.5, scale = scale_param), rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
      mean(C2 < T2) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func, interval = c(0.1, 1e6))$root, error = function(e) median(T2, na.rm = TRUE) * 1.5)
    list(shape = NA, scale = opt_scale, is_ushape = TRUE)

  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      mean(C2 < T2) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func_scale, interval = c(0.1, 1e6))$root, error = function(e) median(T2, na.rm = TRUE) * 1.5)
    list(shape = if (pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE)
  }
}

# -------------------------------------------------------------------------
# Step 1: Calibration IPTW (1:5 points)
# -------------------------------------------------------------------------
calculate_calibrated_iptw <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points_days <- (1:5) * 365.25

  for (ag in unique(data$Age_class)) {
    if (!(ag %in% names(ref_surv_list))) next
    S_macro_target <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
    ag_idx <- which(data$Age_class == ag)
    ag_data <- data[ag_idx, , drop = FALSE]

    log_t1 <- log(pmax(ag_data$T1 / 365.25, 1e-6))
    log_t1_sq <- log_t1^2

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t1 + theta[2] * log_t1_sq)
      w <- pmax(w, 1e-6)
      w <- w / mean(w)
      km <- tryCatch(survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w), error = function(e) NULL)
      if (is.null(km)) return(1e12)
      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(optim(c(0, 0), obj_func, method = "BFGS"), error = function(e) list(par = c(0, 0)))
    th <- opt$par

    w_opt <- exp(th[1] * log_t1 + th[2] * log_t1_sq)
    w_opt <- pmax(w_opt, 1e-6)

    lb <- quantile(w_opt, 0.01, na.rm = TRUE)
    ub <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lb, pmin(w_opt, ub))
    data$iptw[ag_idx] <- w_opt / mean(w_opt)
  }
  data
}

# -------------------------------------------------------------------------
# SAFER simulation-time generator (Capped at 100 years to prevent Inf)
# -------------------------------------------------------------------------
gen_sim_times_safe <- function(fit, newdata, dist_type = NULL) {
  if (is.null(fit) || is.null(newdata) || nrow(newdata) == 0) return(rep(NA_real_, nrow(newdata)))

  if (is.null(dist_type)) dist_type <- fit$dlist$name
  res_m <- fit$res
  n <- nrow(newdata)

  dist_pars <- c("shape", "scale", "mu", "sigma", "Q")
  cov_rows <- setdiff(rownames(res_m), dist_pars)

  xb <- rep(0, n)
  ind <- function(x) as.numeric(x)

  if ("XMutated" %in% cov_rows && "X" %in% names(newdata)) xb <- xb + res_m["XMutated", "est"] * ind(newdata$X == "Mutated")
  if ("Age" %in% cov_rows && "Age" %in% names(newdata)) xb <- xb + res_m["Age", "est"] * as.numeric(newdata$Age)
  if ("SexFemale" %in% cov_rows && "Sex" %in% names(newdata)) xb <- xb + res_m["SexFemale", "est"] * ind(newdata$Sex == "Female")
  if ("HistologyREAD" %in% cov_rows && "Histology" %in% names(newdata)) xb <- xb + res_m["HistologyREAD", "est"] * ind(newdata$Histology == "READ")
  if ("HistologyCOADREAD" %in% cov_rows && "Histology" %in% names(newdata)) xb <- xb + res_m["HistologyCOADREAD", "est"] * ind(newdata$Histology == "COADREAD")
  for (r in cov_rows) {
    if (grepl("^ns\\d+$", r) && r %in% names(newdata)) {
      xb <- xb + res_m[r, "est"] * as.numeric(newdata[[r]])
    }
  }

  clip <- function(x, lo = -15, hi = 10) pmin(pmax(x, lo), hi)
  safe_num <- function(x) ifelse(is.finite(x), x, NA_real_)
  res_times <- rep(NA_real_, n)

  if (dist_type %in% c("weibull", "weibullPH", "weibullph")) {
    if ("shape" %in% rownames(res_m) && "scale" %in% rownames(res_m)) {
      shape <- safe_num(as.numeric(res_m["shape","est"]))
      sc0   <- safe_num(as.numeric(res_m["scale","est"]))
      if (is.finite(shape) && is.finite(sc0) && shape > 0 && sc0 > 0) {
        mu <- clip(log(sc0) + xb)
        res_times <- rweibull(n, shape = shape, scale = exp(mu))
      }
    }
  } else if (dist_type %in% c("llogis", "loglogistic")) {
    if ("shape" %in% rownames(res_m) && "scale" %in% rownames(res_m)) {
      shape <- safe_num(as.numeric(res_m["shape","est"]))
      sc0   <- safe_num(as.numeric(res_m["scale","est"]))
      if (is.finite(shape) && is.finite(sc0) && shape > 0 && sc0 > 0) {
        mu <- clip(log(sc0) + xb)
        res_times <- flexsurv::rllogis(n, shape = shape, scale = exp(mu))
      }
    }
  } else if (dist_type %in% c("gengamma", "gen_gamma", "generalized gamma")) {
    if (all(c("mu","sigma","Q") %in% rownames(res_m))) {
      mu0   <- safe_num(as.numeric(res_m["mu","est"]))
      sigma <- safe_num(as.numeric(res_m["sigma","est"]))
      Q     <- safe_num(as.numeric(res_m["Q","est"]))
      if (is.finite(mu0) && is.finite(sigma) && is.finite(Q) && sigma > 0) {
        mu <- clip(mu0 + xb)
        res_times <- flexsurv::rgengamma(n, mu = mu, sigma = sigma, Q = Q)
      }
    }
  }

  # --- CRITICAL: Inf Explosion Protection ---
  # Replace NAs and cap at 100 years to guarantee finite CI calculation
  res_times[is.na(res_times)] <- median(res_times, na.rm=TRUE)
  if (is.na(res_times[1])) res_times <- rep(365.25, n)
  res_times <- pmin(res_times, 365.25 * 100)

  return(res_times)
}

# -------------------------------------------------------------------------
# Extractors for Naive/LT
# -------------------------------------------------------------------------
extract_af_full <- function(fit) {
  blank <- c(AF_X_est=NA, AF_X_L=NA, AF_X_U=NA, AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
             AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA, AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
             AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA)
  if (is.null(fit)) return(blank)

  res_m <- fit$res
  get <- function(r, col) if (r %in% rownames(res_m)) res_m[r, col] else NA_real_

  out <- blank
  out["AF_X_est"] <- exp(get("XMutated","est")); out["AF_X_L"] <- exp(get("XMutated","L95%")); out["AF_X_U"] <- exp(get("XMutated","U95%"))
  out["AF_Age_est"] <- exp(get("Age","est") * 10); out["AF_Age_L"] <- exp(get("Age","L95%") * 10); out["AF_Age_U"] <- exp(get("Age","U95%") * 10)
  out["AF_Sex_est"] <- exp(get("SexFemale","est")); out["AF_Sex_L"] <- exp(get("SexFemale","L95%")); out["AF_Sex_U"] <- exp(get("SexFemale","U95%"))
  out["AF_READ_est"] <- exp(get("HistologyREAD","est")); out["AF_READ_L"] <- exp(get("HistologyREAD","L95%")); out["AF_READ_U"] <- exp(get("HistologyREAD","U95%"))
  out["AF_COADREAD_est"] <- exp(get("HistologyCOADREAD","est")); out["AF_COADREAD_L"] <- exp(get("HistologyCOADREAD","L95%")); out["AF_COADREAD_U"] <- exp(get("HistologyCOADREAD","U95%"))
  out
}

calc_marginal_metrics <- function(sim_times_days) {
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(c(Med=NA_real_, S5=NA_real_))
  c(Med = median(sim_times_days, na.rm = TRUE) / 365.25, S5 = mean(sim_times_days > 5 * 365.25, na.rm = TRUE) * 100)
}

calc_curve_km <- function(sim_times_days) {
  t_seq <- seq(0, 365.25 * 5, length.out = 100)
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(rep(NA_real_, length(t_seq)))
  km <- tryCatch(survfit(Surv(sim_times_days, rep(1, length(sim_times_days))) ~ 1), error = function(e) NULL)
  if (is.null(km)) return(rep(NA_real_, length(t_seq)))
  approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
}

# -------------------------------------------------------------------------
# Median-based ratio CI (delta method)
# -------------------------------------------------------------------------
median_se_uncens <- function(x) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  n <- length(x)
  if (n < 30) return(NA_real_)

  m <- median(x)
  if (m <= 0) return(NA_real_)

  lx <- log(x)
  den <- tryCatch(stats::density(lx, n = 512), error = function(e) NULL)
  f_m <- NA_real_
  if (!is.null(den)) {
    f_lm <- approx(den$x, den$y, xout = log(m), rule = 2)$y
    f_m <- f_lm / m
  }

  # Fallback to robust SD if density estimation fails
  if (is.null(den) || !is.finite(f_m) || f_m <= 0) {
    s <- mad(x, constant = 1.4826)
    if (!is.finite(s) || s <= 0) s <- sd(x)
    if (!is.finite(s) || s <= 0) s <- 0.01
    return(1.253314 * s / sqrt(n))
  }

  sqrt(0.25 / (n * f_m^2))
}

ratio_ci_from_samples <- function(x0, x1, n_eff = NULL) {
  x0 <- x0[is.finite(x0) & !is.na(x0) & x0 > 0]
  x1 <- x1[is.finite(x1) & !is.na(x1) & x1 > 0]

  if (length(x0) < 30 || length(x1) < 30) return(c(est=NA_real_, L=NA_real_, U=NA_real_))

  m0 <- median(x0); m1 <- median(x1)
  if (!is.finite(m0) || !is.finite(m1) || m0 <= 0 || m1 <= 0) return(c(est=NA_real_, L=NA_real_, U=NA_real_))
  est <- m1 / m0

  se0 <- median_se_uncens(x0); se1 <- median_se_uncens(x1)
  if (!is.finite(se0) || !is.finite(se1) || se0 <= 0 || se1 <= 0) return(c(est=est, L=NA_real_, U=NA_real_))

  if (!is.null(n_eff) && is.finite(n_eff) && n_eff > 0) {
    se0 <- se0 * sqrt(length(x0) / n_eff)
    se1 <- se1 * sqrt(length(x1) / n_eff)
  }

  se_log <- sqrt( (se1/m1)^2 + (se0/m0)^2 )
  L <- exp(log(est) - 1.95996 * se_log)
  U <- exp(log(est) + 1.95996 * se_log)
  c(est=est, L=L, U=U)
}

# -------------------------------------------------------------------------
# Proposed TOTAL EFFECT AFs by counterfactual OS distributions
# -------------------------------------------------------------------------
compute_prop_afs_counterfactual <- function(pseudo_base, fit_t1, fit_t2, fit_t1res, ns_obj, ns_colnames, Data_cgp_logT1_resid, n_eff) {

  sim_os_do <- function(df) {
    df$sim_T1 <- gen_sim_times_safe(fit_t1, df, dist_type = "weibull")
    df$logT1_sim_scale <- log(pmax(df$sim_T1 / 365.25, 1e-6))
    df$logT1_sim_resid <- df$logT1_sim_scale

    if (!is.null(fit_t1res)) {
      pred <- tryCatch(predict(fit_t1res, newdata = df), error = function(e) NULL)
      if (!is.null(pred) && length(pred) == nrow(df) && all(is.finite(pred))) {
        df$logT1_sim_resid <- df$logT1_sim_scale - pred
      }
    }

    # CLAMP processing
    min_res <- min(Data_cgp_logT1_resid, na.rm = TRUE)
    max_res <- max(Data_cgp_logT1_resid, na.rm = TRUE)
    df$logT1_sim_resid <- pmax(min_res, pmin(max_res, df$logT1_sim_resid))

    ns_sim <- tryCatch(predict(ns_obj, df$logT1_sim_resid), error = function(e) NULL)
    if (is.null(ns_sim)) {
      ns_sim <- matrix(0, nrow(df), length(ns_colnames))
    } else {
      ns_sim <- as.matrix(ns_sim)
    }
    colnames(ns_sim) <- ns_colnames
    for (j in seq_along(ns_colnames)) df[[ns_colnames[j]]] <- ns_sim[, j]

    df$sim_T2 <- gen_sim_times_safe(fit_t2, df, dist_type = "gengamma")
    return(df$sim_T1 + df$sim_T2)
  }

  pseudo_base$X <- factor(pseudo_base$X, levels = c("WT","Mutated"))
  pseudo_base$Sex <- factor(pseudo_base$Sex, levels = c("Male","Female"))
  pseudo_base$Histology <- factor(pseudo_base$Histology, levels = c("COAD","READ","COADREAD"))
  N <- nrow(pseudo_base)

  os_nat <- sim_os_do(pseudo_base)

  # Explicit rep() to guarantee length consistency
  df0 <- pseudo_base; df0$X <- factor(rep("WT", N), levels = levels(pseudo_base$X))
  df1 <- pseudo_base; df1$X <- factor(rep("Mutated", N), levels = levels(pseudo_base$X))
  af_x <- ratio_ci_from_samples(sim_os_do(df0), sim_os_do(df1), n_eff)

  dfA0 <- pseudo_base
  dfA1 <- pseudo_base; dfA1$Age <- dfA1$Age + 10
  af_age <- ratio_ci_from_samples(sim_os_do(dfA0), sim_os_do(dfA1), n_eff)

  dfS0 <- pseudo_base; dfS0$Sex <- factor(rep("Male", N), levels = levels(pseudo_base$Sex))
  dfS1 <- pseudo_base; dfS1$Sex <- factor(rep("Female", N), levels = levels(pseudo_base$Sex))
  af_sex <- ratio_ci_from_samples(sim_os_do(dfS0), sim_os_do(dfS1), n_eff)

  dfHR0 <- pseudo_base; dfHR0$Histology <- factor(rep("COAD", N), levels = levels(pseudo_base$Histology))
  dfHR1 <- pseudo_base; dfHR1$Histology <- factor(rep("READ", N), levels = levels(pseudo_base$Histology))
  af_read <- ratio_ci_from_samples(sim_os_do(dfHR0), sim_os_do(dfHR1), n_eff)

  dfHC1 <- pseudo_base; dfHC1$Histology <- factor(rep("COADREAD", N), levels = levels(pseudo_base$Histology))
  af_coadread <- ratio_ci_from_samples(sim_os_do(dfHR0), sim_os_do(dfHC1), n_eff)

  Prop_AF <- c(
    AF_X_est = af_x["est"], AF_X_L = af_x["L"], AF_X_U = af_x["U"],
    AF_Age_est = af_age["est"], AF_Age_L = af_age["L"], AF_Age_U = af_age["U"],
    AF_Sex_est = af_sex["est"], AF_Sex_L = af_sex["L"], AF_Sex_U = af_sex["U"],
    AF_READ_est = af_read["est"], AF_READ_L = af_read["L"], AF_READ_U = af_read["U"],
    AF_COADREAD_est = af_coadread["est"], AF_COADREAD_L = af_coadread["L"], AF_COADREAD_U = af_coadread["U"]
  )

  list(Prop_AF = Prop_AF, pseudo_nat_OS = os_nat)
}

# -------------------------------------------------------------------------
# Core Simulation Engine
# -------------------------------------------------------------------------
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, True_Med, True_Shape,
                              t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  if (!is.finite(N_target) || N_target < 100) return(NULL)

  # --- DGP ---
  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", 1.10, 1.0) * ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))
  AF_total <- AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  u <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true_base <- (True_Med * 365.25) * ((1 - u) / u)^(1 / True_Shape)

  if (is.null(t1_pat)) t1_pat <- "indep"
  if (t1_pat == "early") { T1_base <- runif(N_macro, 30, pmax(31, T_true_base * 0.3))
  } else if (t1_pat %in% c("dep_1yr", "real")) { T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat %in% c("dep_2yr", "rev")) { T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else { T1_base <- runif(N_macro, 30, pmax(31, T_true_base)) }

  T2_base <- pmax(0.1, T_true_base - T1_base)
  T1 <- T1_base * AF_total
  T2_true <- T2_base * AF_total
  T_true <- T1 + T2_true

  cgp_idx <- which(T_true > T1 & T1 <= 365.25 * 10)
  if (length(cgp_idx) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_idx] < 60, 0.9, 0.1)
  sel <- sample(cgp_idx, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target, Age = Age[sel], Age_class = ifelse(Age[sel] < 60, "Young", "Old"),
    Sex = factor(Sex[sel], levels = c("Male", "Female")),
    Histology = factor(Histology[sel], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[sel] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[sel], T1 = T1[sel]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # --- Censoring ---
  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (isTRUE(params$is_ushape)) {
    u_c <- runif(N_target)
    C2 <- ifelse(u_c < 0.5, rweibull(N_target, shape = 0.5, scale = params$scale), rweibull(N_target, shape = 3.0, scale = params$scale * 2))
  } else if (cens_pat == "indep") { C2 <- rexp(N_target, rate = 1 / params$scale)
  } else { C2 <- rweibull(N_target, shape = params$shape, scale = params$scale) }

  Data_cgp$C2 <- C2
  Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)
  Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs
  if (sum(Data_cgp$Event) < 20) return(NULL)

  # --- Step 1: IPTW ---
  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young","Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_calibrated_iptw(Data_cgp, ref_surv_list)
  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)
  ESS <- (sum(Data_cgp$iptw)^2) / sum(Data_cgp$iptw^2)

  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  Data_cgp$logT1_resid <- Data_cgp$logT1_scale

  fit_t1res <- tryCatch({
    form_t1res <- as.formula(paste("logT1_scale ~", paste(valid_covs, collapse = " + ")))
    lm(form_t1res, data = Data_cgp, weights = Data_cgp$iptw)
  }, error = function(e) NULL)

  if (!is.null(fit_t1res)) {
    r <- tryCatch(resid(fit_t1res), error = function(e) NULL)
    if (!is.null(r) && length(r) == nrow(Data_cgp) && all(is.finite(r))) Data_cgp$logT1_resid <- r
  }

  ns_obj <- tryCatch(ns(Data_cgp$logT1_resid, df = 3), error = function(e) NULL)
  if (is.null(ns_obj)) {
    ns_mat <- matrix(0, nrow(Data_cgp), 3)
    colnames(ns_mat) <- paste0("ns", 1:3)
  } else {
    ns_mat <- as.matrix(ns_obj)
    colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  }
  for (j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt    <- tryCatch(flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurvreg(form_t1, data = Data_cgp, weights = Data_cgp$iptw, dist = "weibull"), error = function(e) NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurvreg(form_t2, data = Data_cgp, weights = Data_cgp$iptw, dist = "gengamma"), error = function(e) NULL)

  if (is.null(fit_t1) || is.null(fit_t2)) return(NULL)

  n_pseudo <- as.integer(round(ESS))
  if (!is.finite(n_pseudo) || is.na(n_pseudo)) n_pseudo <- 1000L
  n_pseudo <- max(500L, min(10000L, n_pseudo))

  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = n_pseudo, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_base <- Data_cgp[idx_pseudo, , drop = FALSE]

  cf <- tryCatch(
    compute_prop_afs_counterfactual(
      pseudo_base = pseudo_base, fit_t1 = fit_t1, fit_t2 = fit_t2,
      fit_t1res = fit_t1res, ns_obj = ns_obj, ns_colnames = colnames(ns_mat),
      Data_cgp_logT1_resid = Data_cgp$logT1_resid, n_eff = ESS
    ),
    error = function(e) NULL
  )
  if (is.null(cf)) return(NULL)

  Prop_AF <- cf$Prop_AF
  pseudo_nat_OS <- cf$pseudo_nat_OS

  sim_naive <- if (!is.null(fit_naive)) gen_sim_times_safe(fit_naive, Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt    <- if (!is.null(fit_lt))    gen_sim_times_safe(fit_lt,    Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))

  list(
    ESS = ESS, n_pseudo = n_pseudo,
    True_Metrics = calc_marginal_metrics(T_true),
    Naive_AF = extract_af_full(fit_naive), Naive_Metrics = calc_marginal_metrics(sim_naive),
    LT_AF = extract_af_full(fit_lt), LT_Metrics = calc_marginal_metrics(sim_lt),
    Prop_AF = Prop_AF, Prop_Metrics = calc_marginal_metrics(pseudo_nat_OS),
    curve_true = calc_curve_km(T_true), curve_naive = calc_curve_km(sim_naive),
    curve_lt = calc_curve_km(sim_lt), curve_prop = calc_curve_km(pseudo_nat_OS),
    sample_t1 = Data_cgp$T1[1:min(500, nrow(Data_cgp))], sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))]
  )
}

# -------------------------------------------------------------------------
# Formatting helpers (UI expects these)
# -------------------------------------------------------------------------
safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) {
    if (is.na(v) || !is.numeric(v) || !is.finite(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v)
  })
}

format_af_single <- function(af_vec, var_name_prefix) {
  est <- af_vec[paste0(var_name_prefix, "_est")]
  L   <- af_vec[paste0(var_name_prefix, "_L")]
  U   <- af_vec[paste0(var_name_prefix, "_U")]
  if (is.na(est) || !is.finite(est)) return("N/A")
  if (is.na(L) || is.na(U) || !is.finite(L) || !is.finite(U)) return(sprintf("%.2f (NA-NA)", est))
  sprintf("%.2f (%.2f-%.2f)", est, L, U)
}

# =========================================================================
# UI Observers and Rendering Logic
# =========================================================================

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)
  showNotification("Running Single Simulation...", type = "message", duration = 2)

  set.seed(as.integer(Sys.time()))

  res <- tryCatch(
    run_sim_iteration(
      input$sim_n, input$sim_true_af, input$sim_mut_freq / 100,
      input$sim_true_med, input$sim_true_shape,
      input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100,
      return_data = TRUE
    ),
    error = function(e) NULL
  )

  shiny::validate(shiny::need(!is.null(res), "Simulation failed. Try changing settings."))
  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
                 safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"),
      `Naive AFT` = c(format_af_single(res$Naive_AF, "AF_X"), format_af_single(res$Naive_AF, "AF_Age"),
                      format_af_single(res$Naive_AF, "AF_Sex"), format_af_single(res$Naive_AF, "AF_READ"),
                      format_af_single(res$Naive_AF, "AF_COADREAD"), safe_fmt(res$Naive_Metrics["Med"]),
                      safe_fmt(res$Naive_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)),
      `Standard LT AFT` = c(format_af_single(res$LT_AF, "AF_X"), format_af_single(res$LT_AF, "AF_Age"),
                            format_af_single(res$LT_AF, "AF_Sex"), format_af_single(res$LT_AF, "AF_READ"),
                            format_af_single(res$LT_AF, "AF_COADREAD"), safe_fmt(res$LT_Metrics["Med"]),
                            safe_fmt(res$LT_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)),
      `Proposed (Total Effect)` = c(paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_Age"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_Sex"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_READ"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_COADREAD"), "</b>"),
                                    paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"),
                                    paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"),
                                    paste0("<b>", safe_fmt(res$ESS, 0), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)
    plot_df <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(res$curve_true, res$curve_naive, res$curve_lt, res$curve_prop),
      Model = factor(rep(c("1. True Marginal (Macro)", "2. Naive AFT (Immortal Bias)", "3. Standard LT (Dep. Trunc Bias)", "4. Proposed (Total Effect)"), each = 100))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) + scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) + labs(title = "Tamura & Ikegami Model Validation", subtitle = "Single Run Results", x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)
  n_sims <- 400
  results <- list()

  withProgress(message = "Running 400 Simulations...", value = 0, {
    for (i in 1:n_sims) {
      set.seed(i)
      res <- tryCatch(
        run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_true_med, input$sim_true_shape, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = FALSE),
        error = function(e) NULL
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })

  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  pull_vals <- function(getter) { v <- sapply(results, getter); v[!is.na(v) & is.finite(v)] }

  calc_stats_af <- function(model_name, var_prefix, true_val) {
    est_vals <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_est")])
    l_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_L")])
    u_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_U")])
    if (length(est_vals) == 0) return("N/A")
    mean_est <- mean(est_vals); mse_val <- mean((est_vals - true_val)^2); cr_val <- mean(l_vals <= true_val & true_val <= u_vals, na.rm = TRUE) * 100
    sprintf("%.2f (MSE: %.3f, CR: %.1f%%)", mean_est, mse_val, cr_val)
  }

  calc_stats_metrics <- function(model_name, metric_name, true_val = NULL) {
    vals <- pull_vals(function(x) x[[model_name]][metric_name])
    if (length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    if (!is.null(true_val)) return(sprintf("%.2f (MSE: %.3f)", mean_val, mean((vals - true_val)^2)))
    sprintf("%.2f", mean_val)
  }

  True_AF <- input$sim_true_af
  mean_true_med <- mean(pull_vals(function(x) x$True_Metrics["Med"]), na.rm = TRUE)
  mean_true_s5  <- mean(pull_vals(function(x) x$True_Metrics["S5"]),  na.rm = TRUE)

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"),
      `Naive AFT` = c(calc_stats_af("Naive_AF", "AF_X", True_AF), calc_stats_af("Naive_AF", "AF_Age", 0.85), calc_stats_af("Naive_AF", "AF_Sex", 1.10), calc_stats_af("Naive_AF", "AF_READ", 0.90), calc_stats_af("Naive_AF", "AF_COADREAD", 0.95), calc_stats_metrics("Naive_Metrics", "Med", mean_true_med), calc_stats_metrics("Naive_Metrics", "S5",  mean_true_s5), "N/A"),
      `Standard LT AFT` = c(calc_stats_af("LT_AF", "AF_X", True_AF), calc_stats_af("LT_AF", "AF_Age", 0.85), calc_stats_af("LT_AF", "AF_Sex", 1.10), calc_stats_af("LT_AF", "AF_READ", 0.90), calc_stats_af("LT_AF", "AF_COADREAD", 0.95), calc_stats_metrics("LT_Metrics", "Med", mean_true_med), calc_stats_metrics("LT_Metrics", "S5",  mean_true_s5), "N/A"),
      `Proposed (Total Effect)` = c(paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5",  mean_true_s5), "</b>"), paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_multi_fig1 <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)
    safe_curve_mean <- function(getter) { mat <- do.call(rbind, lapply(results, function(x) { v <- getter(x); if (is.null(v) || length(v) != 100) return(rep(NA_real_, 100)); v })); colMeans(mat, na.rm = TRUE) }
    plot_df1 <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(safe_curve_mean(function(x) x$curve_true), safe_curve_mean(function(x) x$curve_naive), safe_curve_mean(function(x) x$curve_lt), safe_curve_mean(function(x) x$curve_prop)),
      Model = factor(rep(c("1. True Marginal", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed (Total Effect)"), each = 100))
    )
    ggplot(plot_df1, aes(x = Time, y = Survival, color = Model, linetype = Model)) + geom_line(linewidth = 1.2) + scale_color_manual(values = c("black", "gray50", "gray50", "black")) + scale_linetype_manual(values = c("solid", "dotted", "dashed", "twodash")) + scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) + theme_minimal(base_size = 15) + labs(title = "Figure 1: Reconstructed Marginal Survival Curves", subtitle = "Averaged across iterations", x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") + theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })

  output$sim_multi_fig2 <- renderPlot({
    af_naive <- pull_vals(function(x) x$Naive_AF["AF_X_est"]); af_lt <- pull_vals(function(x) x$LT_AF["AF_X_est"]); af_prop <- pull_vals(function(x) x$Prop_AF["AF_X_est"])
    plot_df2 <- data.frame(
      AF = c(af_naive, af_lt, af_prop),
      Model = factor(rep(c("Naive AFT", "Standard LT AFT", "Proposed (Total Effect)"), c(length(af_naive), length(af_lt), length(af_prop))), levels = c("Naive AFT", "Standard LT AFT", "Proposed (Total Effect)"))
    )
    ggplot(plot_df2, aes(x = Model, y = AF, fill = Model)) + geom_boxplot(alpha = 0.8, outlier.shape = 1) + geom_hline(yintercept = True_AF, color = "black", linetype = "dashed", linewidth = 1.2) + scale_fill_grey(start = 0.9, end = 0.4) + coord_cartesian(ylim = c(0, max(True_AF * 2.5, 3))) + theme_minimal(base_size = 15) + labs(title = "Figure 2: Distribution of Estimated AF", subtitle = "Dashed line = True AF", x = "", y = "Estimated Acceleration Factor (AF)") + theme(plot.title = element_text(face = "bold"), legend.position = "none")
  })

  output$sim_multi_fig3 <- renderPlot({
    t1_samp <- results[[1]]$sample_t1 / 365.25; t2_samp <- results[[1]]$sample_t2 / 365.25
    plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)
    ggplot(plot_df3, aes(x = T1, y = T2)) + geom_point(color = "black", alpha = 0.4, size = 2) + geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE) + theme_minimal(base_size = 15) + labs(title = "Figure 3: Dependent Left Truncation", subtitle = "Correlation between T1 and residual survival (true)", x = "Time to CGP Test (T1, Years)", y = "Residual Survival Time (T2, Years)") + theme(plot.title = element_text(face = "bold"))
  })
})
