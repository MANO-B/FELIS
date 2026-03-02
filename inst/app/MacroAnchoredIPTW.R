# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model)
#  - 1:5 (1-5y) calibration 유지（維持）
#  - Calibration design tweak: preserve X distribution within Age_class (reduce attenuation)
#  - T2 model: gengamma
#  - G-computation: T2 also generated from gengamma
#  - Proposed: TOTAL EFFECT estimated by counterfactual OS distributions (do-interventions),
#              NOT by final AFT coefficient (since OS is not llogis)
#  - CI: delta-method from median SE (kernel density at median) + pseudo size = ESS (min 500, max 10000)
#  - Robust to crashes (no predict(flexsurvreg, type="lp"))
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
      C2 <- ifelse(
        u_fixed < 0.5,
        rweibull(length(T2), shape = 0.5, scale = scale_param),
        rweibull(length(T2), shape = 3.0, scale = scale_param * 2)
      )
      mean(C2 < T2) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func, interval = c(0.1, 1e6))$root,
                          error = function(e) median(T2, na.rm = TRUE) * 1.5)
    list(shape = NA, scale = opt_scale, is_ushape = TRUE)

  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      mean(C2 < T2) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func_scale, interval = c(0.1, 1e6))$root,
                          error = function(e) median(T2, na.rm = TRUE) * 1.5)
    list(shape = if (pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE)
  }
}

# -------------------------------------------------------------------------
# Step 1: Calibration IPTW (1:5 points)
#  - IMPORTANT tweak to reduce attenuation of X:
#      Add a soft constraint to preserve weighted mean of X within Age_class
#      (prevents calibration weights from shifting X distribution)
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
# SAFER simulation-time generator (NO predict(flexsurvreg, type="lp"))
# Supports: weibull, llogis, gengamma
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

  if ("XMutated" %in% cov_rows && "X" %in% names(newdata)) {
    xb <- xb + res_m["XMutated", "est"] * ind(newdata$X == "Mutated")
  }
  if ("Age" %in% cov_rows && "Age" %in% names(newdata)) {
    xb <- xb + res_m["Age", "est"] * as.numeric(newdata$Age)
  }
  if ("SexFemale" %in% cov_rows && "Sex" %in% names(newdata)) {
    xb <- xb + res_m["SexFemale", "est"] * ind(newdata$Sex == "Female")
  }
  if ("HistologyREAD" %in% cov_rows && "Histology" %in% names(newdata)) {
    xb <- xb + res_m["HistologyREAD", "est"] * ind(newdata$Histology == "READ")
  }
  if ("HistologyCOADREAD" %in% cov_rows && "Histology" %in% names(newdata)) {
    xb <- xb + res_m["HistologyCOADREAD", "est"] * ind(newdata$Histology == "COADREAD")
  }
  for (r in cov_rows) {
    if (grepl("^ns\\d+$", r) && r %in% names(newdata)) {
      xb <- xb + res_m[r, "est"] * as.numeric(newdata[[r]])
    }
  }

  # ---- numeric guard helpers ----
  clip <- function(x, lo = -20, hi = 20) pmin(pmax(x, lo), hi)   # exp(20) ~ 4.85e8 days scale ok
  safe_num <- function(x) ifelse(is.finite(x), x, NA_real_)

  if (dist_type %in% c("weibull", "weibullPH", "weibullph")) {
    if (!("shape" %in% rownames(res_m)) || !("scale" %in% rownames(res_m))) return(rep(NA_real_, n))
    shape <- safe_num(as.numeric(res_m["shape","est"]))
    sc0   <- safe_num(as.numeric(res_m["scale","est"]))
    if (!is.finite(shape) || !is.finite(sc0) || shape <= 0 || sc0 <= 0) return(rep(NA_real_, n))

    mu <- log(sc0) + xb
    mu <- clip(mu)
    return(rweibull(n, shape = shape, scale = exp(mu)))
  }

  if (dist_type %in% c("llogis", "loglogistic")) {
    if (!("shape" %in% rownames(res_m)) || !("scale" %in% rownames(res_m))) return(rep(NA_real_, n))
    shape <- safe_num(as.numeric(res_m["shape","est"]))
    sc0   <- safe_num(as.numeric(res_m["scale","est"]))
    if (!is.finite(shape) || !is.finite(sc0) || shape <= 0 || sc0 <= 0) return(rep(NA_real_, n))

    mu <- log(sc0) + xb
    mu <- clip(mu)
    return(flexsurv::rllogis(n, shape = shape, scale = exp(mu)))
  }

  if (dist_type %in% c("gengamma", "gen_gamma", "generalized gamma")) {
    if (!all(c("mu","sigma","Q") %in% rownames(res_m))) return(rep(NA_real_, n))
    mu0   <- safe_num(as.numeric(res_m["mu","est"]))
    sigma <- safe_num(as.numeric(res_m["sigma","est"]))
    Q     <- safe_num(as.numeric(res_m["Q","est"]))
    if (!is.finite(mu0) || !is.finite(sigma) || !is.finite(Q) || sigma <= 0) return(rep(NA_real_, n))

    mu <- mu0 + xb
    mu <- clip(mu)  # gengammaでも暴走防止（計算安定化）
    return(flexsurv::rgengamma(n, mu = mu, sigma = sigma, Q = Q))
  }

  rep(NA_real_, n)
}

# -------------------------------------------------------------------------
# Extractors for Naive/LT (AFT coefficient-based; for comparators only)
# -------------------------------------------------------------------------
extract_af_full <- function(fit) {
  blank <- c(AF_X_est=NA, AF_X_L=NA, AF_X_U=NA,
             AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
             AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA,
             AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
             AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA)
  if (is.null(fit)) return(blank)

  res_m <- fit$res
  get <- function(r, col) if (r %in% rownames(res_m)) res_m[r, col] else NA_real_

  out <- blank
  out["AF_X_est"] <- exp(get("XMutated","est"))
  out["AF_X_L"]   <- exp(get("XMutated","L95%"))
  out["AF_X_U"]   <- exp(get("XMutated","U95%"))

  out["AF_Age_est"] <- exp(get("Age","est") * 10)
  out["AF_Age_L"]   <- exp(get("Age","L95%") * 10)
  out["AF_Age_U"]   <- exp(get("Age","U95%") * 10)

  out["AF_Sex_est"] <- exp(get("SexFemale","est"))
  out["AF_Sex_L"]   <- exp(get("SexFemale","L95%"))
  out["AF_Sex_U"]   <- exp(get("SexFemale","U95%"))

  out["AF_READ_est"] <- exp(get("HistologyREAD","est"))
  out["AF_READ_L"]   <- exp(get("HistologyREAD","L95%"))
  out["AF_READ_U"]   <- exp(get("HistologyREAD","U95%"))

  out["AF_COADREAD_est"] <- exp(get("HistologyCOADREAD","est"))
  out["AF_COADREAD_L"]   <- exp(get("HistologyCOADREAD","L95%"))
  out["AF_COADREAD_U"]   <- exp(get("HistologyCOADREAD","U95%"))

  out
}

# -------------------------------------------------------------------------
# Metrics / Curves
# -------------------------------------------------------------------------
calc_marginal_metrics <- function(sim_times_days) {
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(c(Med=NA_real_, S5=NA_real_))
  c(Med = median(sim_times_days, na.rm = TRUE) / 365.25,
    S5  = mean(sim_times_days > 5 * 365.25, na.rm = TRUE) * 100)
}

calc_curve_km <- function(sim_times_days) {
  t_seq <- seq(0, 365.25 * 5, length.out = 100)
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(rep(NA_real_, length(t_seq)))
  km <- tryCatch(survfit(Surv(sim_times_days, rep(1, length(sim_times_days))) ~ 1),
                 error = function(e) NULL)
  if (is.null(km)) return(rep(NA_real_, length(t_seq)))
  approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
}

# -------------------------------------------------------------------------
# Median-based ratio CI (delta method)
#  - For uncensored sample:
#      Var(median) ≈ p(1-p) / (n f(m)^2) with p=0.5
#    f(m) estimated by kernel density on log-time (stabilizes tails)
# -------------------------------------------------------------------------
median_se_uncens <- function(x) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  n <- length(x)
  if (n < 50) return(NA_real_)
  m <- median(x)
  lx <- log(x)
  den <- tryCatch(stats::density(lx, n = 512), error = function(e) NULL)
  if (is.null(den)) return(NA_real_)
  f_lm <- approx(den$x, den$y, xout = log(m), rule = 2)$y  # density of logX at log(m)
  # f_X(m) = f_logX(log m) / m
  f_m <- f_lm / m
  if (!is.finite(f_m) || f_m <= 0) return(NA_real_)
  # p=0.5
  se_m <- sqrt(0.25 / (n * f_m^2))
  se_m
}

ratio_ci_from_samples <- function(xA, xB) {
  xA <- xA[is.finite(xA) & !is.na(xA) & xA > 0]
  xB <- xB[is.finite(xB) & !is.na(xB) & xB > 0]

  # estは「最低限」出す（極端に少ない時だけ NA）
  if (length(xA) < 10 || length(xB) < 10) return(c(est=NA_real_, L=NA_real_, U=NA_real_))

  mA <- median(xA)
  mB <- median(xB)
  if (!is.finite(mA) || !is.finite(mB) || mA <= 0 || mB <= 0) return(c(est=NA_real_, L=NA_real_, U=NA_real_))
  est <- mB / mA

  # CIは作れなければ NA-NA で返す（estは保持）
  seA <- median_se_uncens(xA)
  seB <- median_se_uncens(xB)
  if (!is.finite(seA) || !is.finite(seB) || seA <= 0 || seB <= 0) {
    return(c(est=est, L=NA_real_, U=NA_real_))
  }

  var_log <- (seB / mB)^2 + (seA / mA)^2
  se_log <- sqrt(var_log)
  L <- exp(log(est) - 1.95996 * se_log)
  U <- exp(log(est) + 1.95996 * se_log)
  c(est=est, L=L, U=U)
}

# -------------------------------------------------------------------------
# Helper: build ns terms for a pseudo dataset (given fitted residual model + ns_obj)
# -------------------------------------------------------------------------
build_ns_terms <- function(df, fit_t1res, ns_obj, ns_colnames) {
  df$logT1_scale <- log(pmax(df$sim_T1 / 365.25, 1e-6))
  df$logT1_resid <- df$logT1_scale

  if (!is.null(fit_t1res)) {
    pred <- tryCatch(predict(fit_t1res, newdata = df), error = function(e) NULL)
    if (!is.null(pred) && length(pred) == nrow(df) && all(is.finite(pred))) {
      df$logT1_resid <- df$logT1_scale - pred
    }
  }

  ns_sim <- NULL
  if (!is.null(ns_obj)) {
    ns_sim <- tryCatch(predict(ns_obj, df$logT1_resid), error = function(e) NULL)
  }
  if (is.null(ns_sim)) {
    ns_sim <- matrix(0, nrow(df), length(ns_colnames))
    colnames(ns_sim) <- ns_colnames
  } else {
    ns_sim <- as.matrix(ns_sim)
    colnames(ns_sim) <- ns_colnames
  }

  for (j in seq_len(ncol(ns_sim))) df[[colnames(ns_sim)[j]]] <- ns_sim[, j]
  df
}

# -------------------------------------------------------------------------
# Proposed TOTAL EFFECT AFs by counterfactual OS distributions
#  - X: do(X=Mutated) vs do(X=WT)
#  - Age: +10 years shift vs original (on the same pseudo base)
#  - Sex: Female vs Male (set)
#  - Histology: READ vs COAD; COADREAD vs COAD
# -------------------------------------------------------------------------
compute_prop_afs_counterfactual <- function(pseudo_base, fit_t1, fit_t2, fit_t1res, ns_obj, ns_colnames) {

  # ensure factors
  pseudo_base$X <- factor(pseudo_base$X, levels = c("WT","Mutated"))
  pseudo_base$Sex <- factor(pseudo_base$Sex, levels = c("Male","Female"))
  pseudo_base$Histology <- factor(pseudo_base$Histology, levels = c("COAD","READ","COADREAD"))

  # -------- natural (for Prop curve/metrics) --------
  pseudo_nat <- pseudo_base
  pseudo_nat$sim_T1 <- gen_sim_times_safe(fit_t1, pseudo_nat, dist_type = "weibull")
  pseudo_nat <- build_ns_terms(pseudo_nat, fit_t1res, ns_obj, ns_colnames)
  pseudo_nat$sim_T2 <- gen_sim_times_safe(fit_t2, pseudo_nat, dist_type = "gengamma")
  pseudo_nat$sim_OS <- pseudo_nat$sim_T1 + pseudo_nat$sim_T2

  # -------- X total effect (do) --------
  pseudo0 <- pseudo_base; pseudo0$X <- factor("WT", levels = c("WT","Mutated"))
  pseudo1 <- pseudo_base; pseudo1$X <- factor("Mutated", levels = c("WT","Mutated"))

  pseudo0$sim_T1 <- gen_sim_times_safe(fit_t1, pseudo0, dist_type = "weibull")
  pseudo1$sim_T1 <- gen_sim_times_safe(fit_t1, pseudo1, dist_type = "weibull")

  pseudo0 <- build_ns_terms(pseudo0, fit_t1res, ns_obj, ns_colnames)
  pseudo1 <- build_ns_terms(pseudo1, fit_t1res, ns_obj, ns_colnames)

  pseudo0$sim_T2 <- gen_sim_times_safe(fit_t2, pseudo0, dist_type = "gengamma")
  pseudo1$sim_T2 <- gen_sim_times_safe(fit_t2, pseudo1, dist_type = "gengamma")

  os0 <- pseudo0$sim_T1 + pseudo0$sim_T2
  os1 <- pseudo1$sim_T1 + pseudo1$sim_T2

  af_x <- ratio_ci_from_samples(os0, os1)

  # -------- Age (+10y) effect (shift) --------
  pseudoA0 <- pseudo_base
  pseudoA1 <- pseudo_base; pseudoA1$Age <- pseudoA1$Age + 10

  pseudoA0$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoA0, dist_type = "weibull")
  pseudoA1$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoA1, dist_type = "weibull")

  pseudoA0 <- build_ns_terms(pseudoA0, fit_t1res, ns_obj, ns_colnames)
  pseudoA1 <- build_ns_terms(pseudoA1, fit_t1res, ns_obj, ns_colnames)

  pseudoA0$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoA0, dist_type = "gengamma")
  pseudoA1$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoA1, dist_type = "gengamma")

  osA0 <- pseudoA0$sim_T1 + pseudoA0$sim_T2
  osA1 <- pseudoA1$sim_T1 + pseudoA1$sim_T2
  af_age <- ratio_ci_from_samples(osA0, osA1)

  # -------- Sex (Female vs Male) --------
  pseudoS0 <- pseudo_base; pseudoS0$Sex <- factor("Male", levels = c("Male","Female"))
  pseudoS1 <- pseudo_base; pseudoS1$Sex <- factor("Female", levels = c("Male","Female"))

  pseudoS0$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoS0, dist_type = "weibull")
  pseudoS1$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoS1, dist_type = "weibull")

  pseudoS0 <- build_ns_terms(pseudoS0, fit_t1res, ns_obj, ns_colnames)
  pseudoS1 <- build_ns_terms(pseudoS1, fit_t1res, ns_obj, ns_colnames)

  pseudoS0$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoS0, dist_type = "gengamma")
  pseudoS1$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoS1, dist_type = "gengamma")

  osS0 <- pseudoS0$sim_T1 + pseudoS0$sim_T2
  osS1 <- pseudoS1$sim_T1 + pseudoS1$sim_T2
  af_sex <- ratio_ci_from_samples(osS0, osS1)

  # -------- Histology READ vs COAD --------
  pseudoH0 <- pseudo_base; pseudoH0$Histology <- factor("COAD", levels = c("COAD","READ","COADREAD"))
  pseudoH1 <- pseudo_base; pseudoH1$Histology <- factor("READ", levels = c("COAD","READ","COADREAD"))

  pseudoH0$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoH0, dist_type = "weibull")
  pseudoH1$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoH1, dist_type = "weibull")

  pseudoH0 <- build_ns_terms(pseudoH0, fit_t1res, ns_obj, ns_colnames)
  pseudoH1 <- build_ns_terms(pseudoH1, fit_t1res, ns_obj, ns_colnames)

  pseudoH0$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoH0, dist_type = "gengamma")
  pseudoH1$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoH1, dist_type = "gengamma")

  osH0 <- pseudoH0$sim_T1 + pseudoH0$sim_T2
  osH1 <- pseudoH1$sim_T1 + pseudoH1$sim_T2
  af_read <- ratio_ci_from_samples(osH0, osH1)

  # -------- Histology COADREAD vs COAD --------
  pseudoC0 <- pseudo_base; pseudoC0$Histology <- factor("COAD", levels = c("COAD","READ","COADREAD"))
  pseudoC1 <- pseudo_base; pseudoC1$Histology <- factor("COADREAD", levels = c("COAD","READ","COADREAD"))

  pseudoC0$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoC0, dist_type = "weibull")
  pseudoC1$sim_T1 <- gen_sim_times_safe(fit_t1, pseudoC1, dist_type = "weibull")

  pseudoC0 <- build_ns_terms(pseudoC0, fit_t1res, ns_obj, ns_colnames)
  pseudoC1 <- build_ns_terms(pseudoC1, fit_t1res, ns_obj, ns_colnames)

  pseudoC0$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoC0, dist_type = "gengamma")
  pseudoC1$sim_T2 <- gen_sim_times_safe(fit_t2, pseudoC1, dist_type = "gengamma")

  osC0 <- pseudoC0$sim_T1 + pseudoC0$sim_T2
  osC1 <- pseudoC1$sim_T1 + pseudoC1$sim_T2
  af_coadread <- ratio_ci_from_samples(osC0, osC1)

  # pack to UI format (AF_*_est/L/U)
  Prop_AF <- c(
    AF_X_est = af_x["est"], AF_X_L = af_x["L"], AF_X_U = af_x["U"],
    AF_Age_est = af_age["est"], AF_Age_L = af_age["L"], AF_Age_U = af_age["U"],
    AF_Sex_est = af_sex["est"], AF_Sex_L = af_sex["L"], AF_Sex_U = af_sex["U"],
    AF_READ_est = af_read["est"], AF_READ_L = af_read["L"], AF_READ_U = af_read["U"],
    AF_COADREAD_est = af_coadread["est"], AF_COADREAD_L = af_coadread["L"], AF_COADREAD_U = af_coadread["U"]
  )

  list(
    Prop_AF = Prop_AF,
    pseudo_nat_OS = pseudo_nat$sim_OS
  )
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

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) *
    ifelse(Sex == "Female", 1.10, 1.0) *
    ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))
  AF_total <- AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  u <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true_base <- (True_Med * 365.25) * ((1 - u) / u)^(1 / True_Shape)

  if (is.null(t1_pat)) t1_pat <- "indep"
  if (t1_pat == "early") {
    T1_base <- runif(N_macro, 30, pmax(31, T_true_base * 0.3))
  } else if (t1_pat %in% c("dep_1yr", "real")) {
    T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat %in% c("dep_2yr", "rev")) {
    T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else {
    T1_base <- runif(N_macro, 30, pmax(31, T_true_base))
  }

  T2_base <- pmax(0.1, T_true_base - T1_base)
  T1 <- T1_base * AF_total
  T2_true <- T2_base * AF_total
  T_true <- T1 + T2_true

  cgp_idx <- which(T_true > T1 & T1 <= 365.25 * 10)
  if (length(cgp_idx) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_idx] < 60, 0.9, 0.1)
  sel <- sample(cgp_idx, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target,
    Age = Age[sel],
    Age_class = ifelse(Age[sel] < 60, "Young", "Old"),
    Sex = factor(Sex[sel], levels = c("Male", "Female")),
    Histology = factor(Histology[sel], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[sel] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[sel],
    T1 = T1[sel]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # --- Censoring ---
  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (isTRUE(params$is_ushape)) {
    u_c <- runif(N_target)
    C2 <- ifelse(
      u_c < 0.5,
      rweibull(N_target, shape = 0.5, scale = params$scale),
      rweibull(N_target, shape = 3.0, scale = params$scale * 2)
    )
  } else if (cens_pat == "indep") {
    C2 <- rexp(N_target, rate = 1 / params$scale)
  } else {
    C2 <- rweibull(N_target, shape = params$shape, scale = params$scale)
  }

  Data_cgp$C2 <- C2
  Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)
  Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs
  if (sum(Data_cgp$Event) < 20) return(NULL)

  # --- Step 1: IPTW (1:5 years) ---
  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young","Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_calibrated_iptw(Data_cgp, ref_surv_list)

  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)
  ESS <- (sum(Data_cgp$iptw)^2) / sum(Data_cgp$iptw^2)

  # --- covariates ---
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # --- spline on residualized logT1 (protect X total effect) ---
  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  Data_cgp$logT1_resid <- Data_cgp$logT1_scale

  fit_t1res <- tryCatch({
    form_t1res <- as.formula(paste("logT1_scale ~", paste(valid_covs, collapse = " + ")))
    lm(form_t1res, data = Data_cgp, weights = Data_cgp$iptw)
  }, error = function(e) NULL)

  if (!is.null(fit_t1res)) {
    r <- tryCatch(resid(fit_t1res), error = function(e) NULL)
    if (!is.null(r) && length(r) == nrow(Data_cgp) && all(is.finite(r))) {
      Data_cgp$logT1_resid <- r
    }
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

  # --- comparator models (as-is) ---
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))

  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt    <- tryCatch(flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  # --- Step 2: Decoupled AFT ---
  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- tryCatch(
    flexsurvreg(form_t1, data = Data_cgp, weights = Data_cgp$iptw, dist = "weibull"),
    error = function(e) NULL
  )

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(
    flexsurvreg(form_t2, data = Data_cgp, weights = Data_cgp$iptw, dist = "gengamma"),
    error = function(e) NULL
  )

  if (is.null(fit_t1) || is.null(fit_t2)) return(NULL)

  # --- Step 3: G-computation pseudo size = ESS (CI widening) ---
  # --- Step 3: G-computation pseudo size = ESS (CI widening) ---
  n_pseudo <- as.integer(round(ESS))
  if (!is.finite(n_pseudo) || is.na(n_pseudo)) n_pseudo <- 1000L
  n_pseudo <- max(500L, min(10000L, n_pseudo))

  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = n_pseudo, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_base <- Data_cgp[idx_pseudo, , drop = FALSE]

  # counterfactual TOTAL effects + natural pseudo OS for curve/metrics
  cf <- tryCatch(
    compute_prop_afs_counterfactual(
      pseudo_base = pseudo_base,
      fit_t1 = fit_t1,
      fit_t2 = fit_t2,
      fit_t1res = fit_t1res,
      ns_obj = ns_obj,
      ns_colnames = colnames(ns_mat),
      Data_cgp_logT1_resid = Data_cgp$logT1_resid, # ★これを追加（クランプ用）
      n_eff = ESS                                 # ★ESSペナルティ用
    ),
    error = function(e) NULL
  )
  if (is.null(cf)) return(NULL)

  Prop_AF <- cf$Prop_AF
  pseudo_nat_OS <- cf$pseudo_nat_OS

  # --- curves / metrics ---
  True_Metrics <- calc_marginal_metrics(T_true)
  curve_true   <- calc_curve_km(T_true)

  sim_naive <- if (!is.null(fit_naive)) gen_sim_times_safe(fit_naive, Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt    <- if (!is.null(fit_lt))    gen_sim_times_safe(fit_lt,    Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))

  Naive_Metrics <- calc_marginal_metrics(sim_naive)
  LT_Metrics    <- calc_marginal_metrics(sim_lt)
  Prop_Metrics  <- calc_marginal_metrics(pseudo_nat_OS)

  curve_naive <- calc_curve_km(sim_naive)
  curve_lt    <- calc_curve_km(sim_lt)
  curve_prop  <- calc_curve_km(pseudo_nat_OS)

  out <- list(
    ESS = ESS,
    n_pseudo = n_pseudo,

    True_Metrics = True_Metrics,

    Naive_AF = extract_af_full(fit_naive),
    Naive_Metrics = Naive_Metrics,

    LT_AF = extract_af_full(fit_lt),
    LT_Metrics = LT_Metrics,

    # Proposed TOTAL EFFECT AFs (counterfactual distribution-based)
    Prop_AF = Prop_AF,
    Prop_Metrics = Prop_Metrics,

    curve_true = curve_true,
    curve_naive = curve_naive,
    curve_lt = curve_lt,
    curve_prop = curve_prop,

    sample_t1 = Data_cgp$T1[1:min(500, nrow(Data_cgp))],
    sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))]
  )

  out
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

  shiny::validate(shiny::need(!is.null(res), "Simulation failed (model fitting / generation). Try changing settings."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
                 safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"),
      `Naive AFT` = c(format_af_single(res$Naive_AF, "AF_X"),
                      format_af_single(res$Naive_AF, "AF_Age"),
                      format_af_single(res$Naive_AF, "AF_Sex"),
                      format_af_single(res$Naive_AF, "AF_READ"),
                      format_af_single(res$Naive_AF, "AF_COADREAD"),
                      safe_fmt(res$Naive_Metrics["Med"]), safe_fmt(res$Naive_Metrics["S5"], 1),
                      safe_fmt(input$sim_n, 0)),
      `Standard LT AFT` = c(format_af_single(res$LT_AF, "AF_X"),
                            format_af_single(res$LT_AF, "AF_Age"),
                            format_af_single(res$LT_AF, "AF_Sex"),
                            format_af_single(res$LT_AF, "AF_READ"),
                            format_af_single(res$LT_AF, "AF_COADREAD"),
                            safe_fmt(res$LT_Metrics["Med"]), safe_fmt(res$LT_Metrics["S5"], 1),
                            safe_fmt(input$sim_n, 0)),
      `Proposed (Total Effect)` = c(paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_Age"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_Sex"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_READ"), "</b>"),
                                    paste0("<b>", format_af_single(res$Prop_AF, "AF_COADREAD"), "</b>"),
                                    paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"),
                                    paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"),
                                    paste0("<b>", safe_fmt(res$ESS, 0), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c",
  sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)
    plot_df <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(res$curve_true, res$curve_naive, res$curve_lt, res$curve_prop),
      Model = factor(rep(c("1. True Marginal (Macro)",
                           "2. Naive AFT (Immortal Bias)",
                           "3. Standard LT (Dep. Trunc Bias)",
                           "4. Proposed (Total Effect)"),
                         each = 100))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Validation",
           subtitle = "Single Run Results",
           x = "Time from Diagnosis (Years)",
           y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.direction = "vertical",
            legend.title = element_blank())
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
        run_sim_iteration(
          input$sim_n, input$sim_true_af, input$sim_mut_freq / 100,
          input$sim_true_med, input$sim_true_shape,
          input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100,
          return_data = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })

  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  pull_vals <- function(getter) {
    v <- sapply(results, getter)
    v[!is.na(v) & is.finite(v)]
  }

  calc_stats_af <- function(model_name, var_prefix, true_val) {
    est_vals <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_est")])
    l_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_L")])
    u_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_prefix, "_U")])

    if (length(est_vals) == 0) return("N/A")
    mean_est <- mean(est_vals)
    mse_val <- mean((est_vals - true_val)^2)
    cr_val <- mean(l_vals <= true_val & true_val <= u_vals, na.rm = TRUE) * 100
    sprintf("%.2f (MSE: %.3f, CR: %.1f%%)", mean_est, mse_val, cr_val)
  }

  calc_stats_metrics <- function(model_name, metric_name, true_val = NULL) {
    vals <- pull_vals(function(x) x[[model_name]][metric_name])
    if (length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    if (!is.null(true_val)) {
      mse_val <- mean((vals - true_val)^2)
      return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
    }
    sprintf("%.2f", mean_val)
  }

  True_AF <- input$sim_true_af
  mean_true_med <- mean(pull_vals(function(x) x$True_Metrics["Med"]), na.rm = TRUE)
  mean_true_s5  <- mean(pull_vals(function(x) x$True_Metrics["S5"]),  na.rm = TRUE)

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
                       sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"),
      `Naive AFT` = c(calc_stats_af("Naive_AF", "AF_X", True_AF),
                      calc_stats_af("Naive_AF", "AF_Age", 0.85),
                      calc_stats_af("Naive_AF", "AF_Sex", 1.10),
                      calc_stats_af("Naive_AF", "AF_READ", 0.90),
                      calc_stats_af("Naive_AF", "AF_COADREAD", 0.95),
                      calc_stats_metrics("Naive_Metrics", "Med", mean_true_med),
                      calc_stats_metrics("Naive_Metrics", "S5",  mean_true_s5),
                      "N/A"),
      `Standard LT AFT` = c(calc_stats_af("LT_AF", "AF_X", True_AF),
                            calc_stats_af("LT_AF", "AF_Age", 0.85),
                            calc_stats_af("LT_AF", "AF_Sex", 1.10),
                            calc_stats_af("LT_AF", "AF_READ", 0.90),
                            calc_stats_af("LT_AF", "AF_COADREAD", 0.95),
                            calc_stats_metrics("LT_Metrics", "Med", mean_true_med),
                            calc_stats_metrics("LT_Metrics", "S5",  mean_true_s5),
                            "N/A"),
      `Proposed (Total Effect)` = c(paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"),
                                    paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"),
                                    paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"),
                                    paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"),
                                    paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"),
                                    paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"),
                                    paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5",  mean_true_s5), "</b>"),
                                    paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c",
  sanitize.text.function = function(x) x)

  # Fig1: mean survival curves
  output$sim_multi_fig1 <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)

    safe_curve_mean <- function(getter) {
      mat <- do.call(rbind, lapply(results, function(x) {
        v <- getter(x)
        if (is.null(v) || length(v) != 100) return(rep(NA_real_, 100))
        v
      }))
      colMeans(mat, na.rm = TRUE)
    }

    avg_curve_true  <- safe_curve_mean(function(x) x$curve_true)
    avg_curve_naive <- safe_curve_mean(function(x) x$curve_naive)
    avg_curve_lt    <- safe_curve_mean(function(x) x$curve_lt)
    avg_curve_prop  <- safe_curve_mean(function(x) x$curve_prop)

    plot_df1 <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(avg_curve_true, avg_curve_naive, avg_curve_lt, avg_curve_prop),
      Model = factor(rep(c("1. True Marginal", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed (Total Effect)"), each = 100))
    )

    ggplot(plot_df1, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "gray50", "gray50", "black")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "twodash")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 1: Reconstructed Marginal Survival Curves",
           subtitle = "Averaged across iterations",
           x = "Time from Diagnosis (Years)",
           y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.direction = "vertical",
            legend.title = element_blank())
  })

  # Fig2: boxplot of AF_X
  output$sim_multi_fig2 <- renderPlot({
    af_naive <- pull_vals(function(x) x$Naive_AF["AF_X_est"])
    af_lt    <- pull_vals(function(x) x$LT_AF["AF_X_est"])
    af_prop  <- pull_vals(function(x) x$Prop_AF["AF_X_est"])

    plot_df2 <- data.frame(
      AF = c(af_naive, af_lt, af_prop),
      Model = factor(rep(c("Naive AFT", "Standard LT AFT", "Proposed (Total Effect)"),
                         c(length(af_naive), length(af_lt), length(af_prop))),
                     levels = c("Naive AFT", "Standard LT AFT", "Proposed (Total Effect)"))
    )

    ggplot(plot_df2, aes(x = Model, y = AF, fill = Model)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 1) +
      geom_hline(yintercept = True_AF, color = "black", linetype = "dashed", linewidth = 1.2) +
      scale_fill_grey(start = 0.9, end = 0.4) +
      coord_cartesian(ylim = c(0, max(True_AF * 2.5, 3))) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 2: Distribution of Estimated AF",
           subtitle = "Dashed line = True AF",
           x = "",
           y = "Estimated Acceleration Factor (AF)") +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")
  })

  # Fig3: dependence scatter (T1 vs T2_true)
  output$sim_multi_fig3 <- renderPlot({
    t1_samp <- results[[1]]$sample_t1 / 365.25
    t2_samp <- results[[1]]$sample_t2 / 365.25
    plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)

    ggplot(plot_df3, aes(x = T1, y = T2)) +
      geom_point(color = "black", alpha = 0.4, size = 2) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 3: Dependent Left Truncation",
           subtitle = "Correlation between T1 and residual survival (true)",
           x = "Time to CGP Test (T1, Years)",
           y = "Residual Survival Time (T2, Years)") +
      theme(plot.title = element_text(face = "bold"))
  })
})

# -------------------------------------------------------------------------
# Robust Median SE and Ratio CI (Delta Method)
# -------------------------------------------------------------------------
median_se_uncens <- function(x) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  n <- length(x)
  if (n < 30) return(NA_real_)

  m <- median(x)
  if (m <= 0) return(NA_real_)

  # 対数スケール上でカーネル密度推定を行う（生存時間における最も安定したSE推定）
  lx <- log(x)
  den <- tryCatch(stats::density(lx, n = 512), error = function(e) NULL)

  if (is.null(den)) {
    # Fallback (Densityが失敗した場合)
    s <- mad(x, constant = 1.4826)
    if (!is.finite(s) || s <= 0) s <- sd(x)
    return(if (is.finite(s) && s > 0) 1.253314 * s / sqrt(n) else NA_real_)
  }

  f_lm <- approx(den$x, den$y, xout = log(m), rule = 2)$y
  f_m <- f_lm / m
  if (!is.finite(f_m) || f_m <= 0) return(NA_real_)

  # Medianの標準誤差
  sqrt(0.25 / (n * f_m^2))
}

ratio_ci_from_samples <- function(x0, x1, n_eff = NULL) {
  x0 <- x0[is.finite(x0) & !is.na(x0) & x0 > 0]
  x1 <- x1[is.finite(x1) & !is.na(x1) & x1 > 0]

  if (length(x0) < 30 || length(x1) < 30) return(c(est=NA_real_, L=NA_real_, U=NA_real_))

  m0 <- median(x0); m1 <- median(x1)
  if (!is.finite(m0) || !is.finite(m1) || m0 <= 0 || m1 <= 0) return(c(est=NA_real_, L=NA_real_, U=NA_real_))
  est <- m1 / m0

  se0 <- median_se_uncens(x0)
  se1 <- median_se_uncens(x1)
  if (!is.finite(se0) || !is.finite(se1) || se0 <= 0 || se1 <= 0) return(c(est=est, L=NA_real_, U=NA_real_))

  # Pseudo sample size から本来の ESS レベルへ信頼区間をペナルティ拡張
  if (!is.null(n_eff) && is.finite(n_eff) && n_eff > 0) {
    se0 <- se0 * sqrt(length(x0) / n_eff)
    se1 <- se1 * sqrt(length(x1) / n_eff)
  }

  se_log <- sqrt((se1/m1)^2 + (se0/m0)^2)
  L <- exp(log(est) - 1.95996 * se_log)
  U <- exp(log(est) + 1.95996 * se_log)
  c(est=est, L=L, U=U)
}

# -------------------------------------------------------------------------
# Proposed TOTAL EFFECT AFs by counterfactual OS distributions
# -------------------------------------------------------------------------
compute_prop_afs_counterfactual <- function(pseudo_base, fit_t1, fit_t2, fit_t1res, ns_obj, ns_colnames, Data_cgp_logT1_resid, n_eff) {

  # 内部関数: 与えられたデータフレームに対してOSをシミュレートする
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

    # ★極めて重要なCLAMP処理（外挿によるT2の暴走・0への潰れを防止）★
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

  # 1. Natural OS (Baseline)
  os_nat <- sim_os_do(pseudo_base)

  # 2. X (Mutated vs WT)
  df0 <- pseudo_base; df0$X <- factor("WT", levels = levels(pseudo_base$X))
  df1 <- pseudo_base; df1$X <- factor("Mutated", levels = levels(pseudo_base$X))
  af_x <- ratio_ci_from_samples(sim_os_do(df0), sim_os_do(df1), n_eff)

  # 3. Age (+10y shift)
  dfA0 <- pseudo_base
  dfA1 <- pseudo_base; dfA1$Age <- dfA1$Age + 10
  af_age <- ratio_ci_from_samples(sim_os_do(dfA0), sim_os_do(dfA1), n_eff)

  # 4. Sex (Female vs Male)
  dfS0 <- pseudo_base; dfS0$Sex <- factor("Male", levels = levels(pseudo_base$Sex))
  dfS1 <- pseudo_base; dfS1$Sex <- factor("Female", levels = levels(pseudo_base$Sex))
  af_sex <- ratio_ci_from_samples(sim_os_do(dfS0), sim_os_do(dfS1), n_eff)

  # 5. Histology (READ vs COAD)
  dfHR0 <- pseudo_base; dfHR0$Histology <- factor("COAD", levels = levels(pseudo_base$Histology))
  dfHR1 <- pseudo_base; dfHR1$Histology <- factor("READ", levels = levels(pseudo_base$Histology))
  af_read <- ratio_ci_from_samples(sim_os_do(dfHR0), sim_os_do(dfHR1), n_eff)

  # 6. Histology (COADREAD vs COAD)
  dfHC1 <- pseudo_base; dfHC1$Histology <- factor("COADREAD", levels = levels(pseudo_base$Histology))
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
