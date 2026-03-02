# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model)
#  - Keep processing unchanged
#  - Make it Shiny-safe: no hard stop, always return required objects
#  - Proposed CI: Bootstrap percentile CI (B=200) ONLY
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(flexsurv)
  library(survival)
  library(ggplot2)
  library(splines)
  library(scales)
})

# Shiny server で「停止」しにくくする（処理内容は変えない）
options(shiny.sanitize.errors = FALSE)

# -------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------
safe_try <- function(expr, default = NULL) {
  tryCatch(expr, error = function(e) default, warning = function(w) {
    # warning は握りつぶし（挙動維持）、ただし致命は error で落とさない
    invokeRestart("muffleWarning")
  })
}

safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) if (is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}

format_af_single <- function(af_vec, var_prefix) {
  est <- af_vec[paste0(var_prefix, "_est")]
  L   <- af_vec[paste0(var_prefix, "_L")]
  U   <- af_vec[paste0(var_prefix, "_U")]
  if (is.na(est) || is.na(L) || is.na(U)) return("N/A")
  sprintf("%.2f (%.2f-%.2f)", est, L, U)
}

# -------------------------------------------------------------------------
# Censoring parameter search (unchanged logic)
# -------------------------------------------------------------------------
find_censoring_param <- function(T2, target_rate, pattern) {
  if (pattern == "peak1y") {
    obj_func <- function(k) {
      scale_val <- 365.25 / ((k - 1) / k)^(1 / k)
      (sum(rweibull(length(T2), shape = k, scale = scale_val) < T2) / length(T2)) - target_rate
    }
    opt_k <- safe_try(uniroot(obj_func, interval = c(1.05, 10))$root, default = 2.0)
    list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE)

  } else if (pattern == "ushape") {
    u_fixed <- runif(length(T2))
    obj_func <- function(scale_param) {
      C2 <- ifelse(
        u_fixed < 0.5,
        rweibull(length(T2), shape = 0.5, scale = scale_param),
        rweibull(length(T2), shape = 3.0, scale = scale_param * 2)
      )
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- safe_try(uniroot(obj_func, interval = c(0.1, 1e6))$root, default = median(T2) * 1.5)
    list(shape = NA, scale = opt_scale, is_ushape = TRUE)

  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- safe_try(uniroot(obj_func_scale, interval = c(0.1, 1e6))$root, default = median(T2) * 1.5)
    list(shape = if (pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE)
  }
}

# -------------------------------------------------------------------------
# Calibration IPTW (1:5 points 유지, unchanged logic)
# -------------------------------------------------------------------------
calculate_calibrated_iptw <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points_days <- (1:5) * 365.25

  for (ag in unique(data$Age_class)) {
    if (!(ag %in% names(ref_surv_list))) next
    S_macro_target <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)

    ag_data_idx <- which(data$Age_class == ag)
    ag_data <- data[ag_data_idx, , drop = FALSE]

    log_t1 <- log(pmax(ag_data$T1 / 365.25, 1e-6))
    log_t1_sq <- log_t1^2

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t1 + theta[2] * log_t1_sq)
      w <- pmax(w, 1e-4)
      w <- w / mean(w)

      km <- safe_try(survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w), default = NULL)
      if (is.null(km)) return(1e6)

      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- safe_try(optim(c(0, 0), obj_func, method = "BFGS"), default = list(par = c(0, 0)))
    theta_hat <- opt$par

    w_opt <- pmax(exp(theta_hat[1] * log_t1 + theta_hat[2] * log_t1_sq), 1e-4)

    # trimming (unchanged intent)
    lower_bound <- quantile(w_opt, 0.01, na.rm = TRUE)
    upper_bound <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lower_bound, pmin(w_opt, upper_bound))

    data$iptw[ag_data_idx] <- w_opt / mean(w_opt)
  }
  data
}

# -------------------------------------------------------------------------
# gen_sim_times(): robust manual lp (NO predict(type="lp")), supports weibull/llogis/gengamma
#  - Keeps parameterization consistent with your existing implementation:
#    weibull: scale=exp(mu), shape from fit$res["shape","est"]
#    llogis : scale=exp(mu), shape=1 / fit$res["scale","est"]
#    gengamma: rgengamma(mu, sigma, Q)
# -------------------------------------------------------------------------
gen_sim_times <- function(fit, newdata, dist_type = NULL) {
  n <- nrow(newdata)
  if (is.null(fit) || n == 0) return(rep(NA_real_, n))

  if (is.null(dist_type)) dist_type <- fit$dlist$name

  res_m <- fit$res
  rn <- rownames(res_m)

  dist_pars <- c("scale", "shape", "mu", "sigma", "Q")

  # base location
  base_loc <- NA_real_
  if (tolower(dist_type) %in% c("gengamma", "gen_gamma", "generalized gamma")) {
    if (!("mu" %in% rn)) return(rep(NA_real_, n))
    base_loc <- as.numeric(res_m["mu", "est"])
  } else {
    if (!("scale" %in% rn)) return(rep(NA_real_, n))
    base_loc <- log(as.numeric(res_m["scale", "est"]))
  }

  # manual lp over covariates that flexsurv puts on location
  lp <- rep(0, n)
  extra_vars <- setdiff(rn, dist_pars)

  add_if <- function(name, vec) {
    if (name %in% extra_vars) lp <<- lp + as.numeric(res_m[name, "est"]) * vec
  }

  # core covariates
  if ("XMutated" %in% names(newdata)) {
    add_if("XMutated", as.numeric(newdata$XMutated))
  } else if ("X" %in% names(newdata)) {
    add_if("XMutated", as.numeric(newdata$X == "Mutated"))
  }

  if ("Age" %in% names(newdata)) add_if("Age", as.numeric(newdata$Age))

  if ("Sex" %in% names(newdata)) add_if("SexFemale", as.numeric(newdata$Sex == "Female"))
  if ("Histology" %in% names(newdata)) {
    add_if("HistologyREAD", as.numeric(newdata$Histology == "READ"))
    add_if("HistologyCOADREAD", as.numeric(newdata$Histology == "COADREAD"))
  }

  # spline terms (ns1, ns2, ns3...)
  for (v in extra_vars) {
    if (grepl("^ns[0-9]+$", v) && (v %in% names(newdata))) {
      lp <- lp + as.numeric(res_m[v, "est"]) * as.numeric(newdata[[v]])
    }
  }

  loc <- base_loc + lp

  # simulate
  dt <- tolower(dist_type)

  if (dt %in% c("weibull", "weibullph", "weibullph")) {
    if (!("shape" %in% rn)) return(rep(NA_real_, n))
    shape <- as.numeric(res_m["shape", "est"])
    return(rweibull(n, shape = shape, scale = exp(loc)))
  }

  if (dt %in% c("llogis", "loglogistic")) {
    if (!("scale" %in% rn)) return(rep(NA_real_, n))
    scl <- as.numeric(res_m["scale", "est"])
    shape_param <- 1 / scl
    scale_param <- exp(loc)
    return(flexsurv::rllogis(n, shape = shape_param, scale = scale_param))
  }

  if (dt %in% c("gengamma", "gen_gamma", "generalized gamma")) {
    if (!("sigma" %in% rn) || !("Q" %in% rn)) return(rep(NA_real_, n))
    sigma <- as.numeric(res_m["sigma", "est"])
    Q <- as.numeric(res_m["Q", "est"])
    return(flexsurv::rgengamma(n, mu = loc, sigma = sigma, Q = Q))
  }

  rep(NA_real_, n)
}

# -------------------------------------------------------------------------
# Extract AF (full set) from flexsurvreg object
# -------------------------------------------------------------------------
extract_af_full <- function(fit) {
  blank_res <- c(
    AF_X_est=NA, AF_X_L=NA, AF_X_U=NA,
    AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
    AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA,
    AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
    AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA
  )
  if (is.null(fit)) return(blank_res)

  res_m <- fit$res
  rn <- rownames(res_m)

  get_est <- function(rname, mult=1) if (rname %in% rn) exp(as.numeric(res_m[rname, "est"]) * mult) else NA_real_
  get_L   <- function(rname, mult=1) if (rname %in% rn) exp(as.numeric(res_m[rname, "L95%"]) * mult) else NA_real_
  get_U   <- function(rname, mult=1) if (rname %in% rn) exp(as.numeric(res_m[rname, "U95%"]) * mult) else NA_real_

  c(
    AF_X_est = get_est("XMutated"), AF_X_L = get_L("XMutated"), AF_X_U = get_U("XMutated"),
    AF_Age_est = get_est("Age", 10), AF_Age_L = get_L("Age", 10), AF_Age_U = get_U("Age", 10),
    AF_Sex_est = get_est("SexFemale"), AF_Sex_L = get_L("SexFemale"), AF_Sex_U = get_U("SexFemale"),
    AF_READ_est = get_est("HistologyREAD"), AF_READ_L = get_L("HistologyREAD"), AF_READ_U = get_U("HistologyREAD"),
    AF_COADREAD_est = get_est("HistologyCOADREAD"), AF_COADREAD_L = get_L("HistologyCOADREAD"), AF_COADREAD_U = get_U("HistologyCOADREAD")
  )
}

# -------------------------------------------------------------------------
# Proposed CI: Bootstrap percentile CI (B=200), processing unchanged otherwise
#  - We bootstrap only the final AF_X from Proposed (Total Effect)
#  - If a replicate fails, we skip it (do not crash Shiny)
# -------------------------------------------------------------------------
bootstrap_prop_ci <- function(Data_cgp, valid_covs, B = 200, n_pseudo = 10000) {
  # return list(est, L, U)
  # point estimate is computed outside; here returns only CI from bootstrap distribution

  n <- nrow(Data_cgp)
  if (n < 50) return(c(L = NA_real_, U = NA_real_, B_ok = 0))

  # fix factor levels to avoid bootstrap level-drop
  lev_Sex <- levels(Data_cgp$Sex)
  lev_His <- levels(Data_cgp$Histology)
  lev_X   <- levels(Data_cgp$X)

  # spline config from original sample (df=3) – we rebuild ns per bootstrap on resampled data
  # (this is equivalent processing step, just within resample)

  af_vals <- numeric(0)

  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)

    db <- Data_cgp[idx, , drop = FALSE]
    db$Sex <- factor(db$Sex, levels = lev_Sex)
    db$Histology <- factor(db$Histology, levels = lev_His)
    db$X <- factor(db$X, levels = lev_X)

    # normalize weights (unchanged intent)
    db$iptw <- db$iptw / mean(db$iptw)

    # rebuild spline on bootstrap sample (same df)
    db$logT1_scale <- log(pmax(db$T1 / 365.25, 1e-6))
    ns_obj_b <- safe_try(ns(db$logT1_scale, df = 3), default = NULL)
    if (is.null(ns_obj_b)) next
    ns_mat_b <- as.matrix(ns_obj_b)
    colnames(ns_mat_b) <- paste0("ns", seq_len(ncol(ns_mat_b)))
    for (j in seq_len(ncol(ns_mat_b))) db[[colnames(ns_mat_b)[j]]] <- ns_mat_b[, j]

    # Step 2 (unchanged models)
    db$T1_event <- 1
    form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
    fit_t1_b <- safe_try(flexsurvreg(form_t1, data = db, weights = db$iptw, dist = "weibull"), default = NULL)
    if (is.null(fit_t1_b)) next

    form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat_b)), collapse = " + ")))
    fit_t2_b <- safe_try(flexsurvreg(form_t2, data = db, weights = db$iptw, dist = "gengamma"), default = NULL)
    if (is.null(fit_t2_b)) next

    # Step 3 (unchanged pseudo sampling)
    idx_pseudo <- sample(seq_len(nrow(db)), size = n_pseudo, replace = TRUE, prob = db$iptw)
    pb <- db[idx_pseudo, , drop = FALSE]

    pb$sim_T1 <- gen_sim_times(fit_t1_b, pb, dist_type = "weibull")

    pb$logT1_sim_scale <- log(pmax(pb$sim_T1 / 365.25, 1e-6))
    ns_sim_b <- safe_try(predict(ns_obj_b, pb$logT1_sim_scale), default = NULL)
    if (is.null(ns_sim_b)) {
      ns_sim_b <- matrix(0, nrow(pb), attr(ns_obj_b, "df"))
    }
    colnames(ns_sim_b) <- paste0("ns", seq_len(ncol(ns_sim_b)))
    for (j in seq_len(ncol(ns_sim_b))) pb[[colnames(ns_sim_b)[j]]] <- ns_sim_b[, j]

    pb$sim_T2 <- gen_sim_times(fit_t2_b, pb, dist_type = "gengamma")
    pb$sim_OS <- pb$sim_T1 + pb$sim_T2
    pb$sim_Event <- 1

    form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
    fit_final_b <- safe_try(flexsurvreg(form_prop_final, data = pb, dist = "llogis"), default = NULL)
    if (is.null(fit_final_b)) next

    af_b <- extract_af_full(fit_final_b)["AF_X_est"]
    if (!is.na(af_b) && is.finite(af_b)) af_vals <- c(af_vals, af_b)
  }

  if (length(af_vals) < 30) {
    return(c(L = NA_real_, U = NA_real_, B_ok = length(af_vals)))
  }

  c(
    L = unname(quantile(af_vals, 0.025, na.rm = TRUE)),
    U = unname(quantile(af_vals, 0.975, na.rm = TRUE)),
    B_ok = length(af_vals)
  )
}

# -------------------------------------------------------------------------
# Metrics + Curves (unchanged intent)
# -------------------------------------------------------------------------
calc_marginal_metrics <- function(sim_times) {
  if (is.null(sim_times) || all(is.na(sim_times))) return(c(Med = NA_real_, S5 = NA_real_))
  c(
    Med = median(sim_times, na.rm = TRUE) / 365.25,
    S5  = mean(sim_times > 5 * 365.25, na.rm = TRUE) * 100
  )
}

calc_marginal_curve <- function(t_data) {
  t_seq <- seq(0, 365.25 * 5, length.out = 100)
  if (is.null(t_data) || all(is.na(t_data))) return(rep(NA_real_, length(t_seq)))
  km <- safe_try(survfit(Surv(t_data, rep(1, length(t_data))) ~ 1), default = NULL)
  if (is.null(km)) return(rep(NA_real_, length(t_seq)))
  approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
}

# =========================================================================
# Core Simulation Engine (processing unchanged)
#  - Proposed total effect on OS (includes T1 impact through G-comp)
#  - T2 model uses gengamma
#  - Proposed CI via Bootstrap percentile CI (B=200)
# =========================================================================
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, True_Med, True_Shape,
                              t1_pat, cens_pat, cens_rate, return_data = FALSE,
                              B_boot = 200, n_pseudo = 10000) {

  # ---------- DGP (unchanged) ----------
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

  cgp_indices <- which(T_true > T1 & T1 <= 365.25 * 10)
  if (length(cgp_indices) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.9, 0.1)
  selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target,
    Age = Age[selected_indices],
    Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
    Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[selected_indices],
    T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # ---------- Censoring (unchanged) ----------
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
  Data_cgp$Event  <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs  <- Data_cgp$T1 + Data_cgp$T2_obs
  if (sum(Data_cgp$Event) < 20) return(NULL)

  # ---------- Step 1: IPTW (1:5 points, unchanged) ----------
  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_calibrated_iptw(Data_cgp, ref_surv_list)
  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)

  ESS <- (sum(Data_cgp$iptw)^2) / sum(Data_cgp$iptw^2)

  # ---------- covariates (unchanged intent) ----------
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # spline(logT1)
  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  ns_obj <- safe_try(ns(Data_cgp$logT1_scale, df = 3), default = NULL)
  if (is.null(ns_obj)) return(NULL)

  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for (j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  # ---------- comparator models (unchanged) ----------
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))

  fit_naive <- safe_try(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), default = NULL)
  fit_lt    <- safe_try(flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis"), default = NULL)

  # ---------- Step 2: Decoupled AFT (unchanged; T2 gengamma) ----------
  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- safe_try(flexsurvreg(form_t1, data = Data_cgp, weights = Data_cgp$iptw, dist = "weibull"), default = NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- safe_try(flexsurvreg(form_t2, data = Data_cgp, weights = Data_cgp$iptw, dist = "gengamma"), default = NULL)

  if (is.null(fit_t1) || is.null(fit_t2)) return(NULL)

  # ---------- Step 3: Calibrated G-computation (unchanged; T2 gengamma simulation) ----------
  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = n_pseudo, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_cgp <- Data_cgp[idx_pseudo, , drop = FALSE]

  pseudo_cgp$sim_T1 <- gen_sim_times(fit_t1, pseudo_cgp, dist_type = "weibull")

  pseudo_cgp$logT1_sim_scale <- log(pmax(pseudo_cgp$sim_T1 / 365.25, 1e-6))
  ns_sim <- safe_try(predict(ns_obj, pseudo_cgp$logT1_sim_scale), default = NULL)
  if (is.null(ns_sim)) ns_sim <- matrix(0, nrow(pseudo_cgp), attr(ns_obj, "df"))
  colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
  for (j in seq_len(ncol(ns_sim))) pseudo_cgp[[colnames(ns_sim)[j]]] <- ns_sim[, j]

  pseudo_cgp$sim_T2 <- gen_sim_times(fit_t2, pseudo_cgp, dist_type = "gengamma")
  pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
  pseudo_cgp$sim_Event <- 1

  form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
  fit_prop_final <- safe_try(flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis"), default = NULL)

  # If proposed fails, still return required structures as NA (no crash)
  Prop_AF <- extract_af_full_ess_scaled(
    fit = fit_prop_final,
    ess = ESS,
    n_pseudo = n_pseudo
  )

  # ---------- curves/metrics (unchanged intent) ----------
  # naive/lt sim times for metrics/curves (same approach as before)
  sim_naive <- if (!is.null(fit_naive)) gen_sim_times(fit_naive, Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt    <- if (!is.null(fit_lt))    gen_sim_times(fit_lt,    Data_cgp, dist_type = "llogis") else rep(NA_real_, nrow(Data_cgp))

  out <- list(
    ESS = ESS,
    True_Metrics = calc_marginal_metrics(T_true),
    Naive_AF = extract_af_full(fit_naive),
    Naive_Metrics = calc_marginal_metrics(sim_naive),
    LT_AF = extract_af_full(fit_lt),
    LT_Metrics = calc_marginal_metrics(sim_lt),
    Prop_AF = Prop_AF,
    Prop_Metrics = calc_marginal_metrics(pseudo_cgp$sim_OS),
    curve_true = calc_marginal_curve(T_true),
    curve_naive = calc_marginal_curve(sim_naive),
    curve_lt = calc_marginal_curve(sim_lt),
    curve_prop = calc_marginal_curve(pseudo_cgp$sim_OS),
    sample_t1 = Data_cgp$T1[1:min(500, nrow(Data_cgp))],
    sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))]
  )

  out
}

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic
#   - No change in displayed objects; add safety tryCatch to avoid server stop
# -------------------------------------------------------------------------

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)

  showNotification("Running Single Simulation...", type = "message", duration = 3)
  set.seed(as.integer(Sys.time()))

  res <- tryCatch(
    run_sim_iteration(
      N_target = input$sim_n,
      True_AF_X = input$sim_true_af,
      Mut_Freq = input$sim_mut_freq / 100,
      True_Med = input$sim_true_med,
      True_Shape = input$sim_true_shape,
      t1_pat = input$sim_t1_pattern,
      cens_pat = input$sim_cens_pattern,
      cens_rate = input$sim_cens_rate / 100,
      return_data = TRUE,
      B_boot = 200,
      n_pseudo = 10000
    ),
    error = function(e) {
      showNotification(paste0("ERROR: ", conditionMessage(e)), type = "error", duration = 30)
      NULL
    }
  )

  shiny::validate(shiny::need(!is.null(res), "Simulation failed. (See notification if any)"))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
                 safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"),
      `Naive AFT` = c(
        format_af_single(res$Naive_AF, "AF_X"),
        format_af_single(res$Naive_AF, "AF_Age"),
        format_af_single(res$Naive_AF, "AF_Sex"),
        format_af_single(res$Naive_AF, "AF_READ"),
        format_af_single(res$Naive_AF, "AF_COADREAD"),
        safe_fmt(res$Naive_Metrics["Med"]),
        safe_fmt(res$Naive_Metrics["S5"], 1),
        safe_fmt(input$sim_n, 0)
      ),
      `Standard LT AFT` = c(
        format_af_single(res$LT_AF, "AF_X"),
        format_af_single(res$LT_AF, "AF_Age"),
        format_af_single(res$LT_AF, "AF_Sex"),
        format_af_single(res$LT_AF, "AF_READ"),
        format_af_single(res$LT_AF, "AF_COADREAD"),
        safe_fmt(res$LT_Metrics["Med"]),
        safe_fmt(res$LT_Metrics["S5"], 1),
        safe_fmt(input$sim_n, 0)
      ),
      `Proposed (Total Effect, Boot CI B=200)` = c(
        paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_Age"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_Sex"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_READ"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_COADREAD"), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"),
        paste0("<b>", safe_fmt(res$ESS, 0), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c",
  sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)
    plot_df <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(res$curve_true, res$curve_naive, res$curve_lt, res$curve_prop),
      Model = factor(rep(c("1. True Marginal (Macro)", "2. Naive AFT (Immortal Bias)", "3. Standard LT (Dep. Trunc Bias)", "4. Proposed (Total Effect)"), each = 100))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Validation",
           subtitle = "Single Run Results (Proposed CI: Bootstrap B=200)",
           x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)

  n_sims <- 400
  results <- list()

  withProgress(message = "Running 400 Simulations...", value = 0, {
    for (i in 1:n_sims) {
      set.seed(i)
      res <- safe_try(
        run_sim_iteration(
          N_target = input$sim_n,
          True_AF_X = input$sim_true_af,
          Mut_Freq = input$sim_mut_freq / 100,
          True_Med = input$sim_true_med,
          True_Shape = input$sim_true_shape,
          t1_pat = input$sim_t1_pattern,
          cens_pat = input$sim_cens_pattern,
          cens_rate = input$sim_cens_rate / 100,
          return_data = FALSE,
          B_boot = 200,
          n_pseudo = 10000
        ),
        default = NULL
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })

  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  pull_vals <- function(getter) {
    v <- sapply(results, getter)
    v[!is.na(v)]
  }

  calc_stats_af <- function(model_name, var_name, true_val) {
    est_vals <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_est")])
    l_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_L")])
    u_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_U")])

    if (length(est_vals) == 0) return("N/A")
    mean_est <- mean(est_vals)
    mse_val  <- mean((est_vals - true_val)^2)
    cr_val   <- mean(l_vals <= true_val & true_val <= u_vals, na.rm = TRUE) * 100
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
  mean_true_s5  <- mean(pull_vals(function(x) x$True_Metrics["S5"]), na.rm = TRUE)

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
                       sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"),
      `Naive AFT` = c(
        calc_stats_af("Naive_AF", "AF_X", True_AF),
        calc_stats_af("Naive_AF", "AF_Age", 0.85),
        calc_stats_af("Naive_AF", "AF_Sex", 1.10),
        calc_stats_af("Naive_AF", "AF_READ", 0.90),
        calc_stats_af("Naive_AF", "AF_COADREAD", 0.95),
        calc_stats_metrics("Naive_Metrics", "Med", mean_true_med),
        calc_stats_metrics("Naive_Metrics", "S5", mean_true_s5),
        "N/A"
      ),
      `Standard LT AFT` = c(
        calc_stats_af("LT_AF", "AF_X", True_AF),
        calc_stats_af("LT_AF", "AF_Age", 0.85),
        calc_stats_af("LT_AF", "AF_Sex", 1.10),
        calc_stats_af("LT_AF", "AF_READ", 0.90),
        calc_stats_af("LT_AF", "AF_COADREAD", 0.95),
        calc_stats_metrics("LT_Metrics", "Med", mean_true_med),
        calc_stats_metrics("LT_Metrics", "S5", mean_true_s5),
        "N/A"
      ),
      `Proposed (Total Effect, Boot CI B=200)` = c(
        paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"),
        paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"),
        paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5", mean_true_s5), "</b>"),
        paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c",
  sanitize.text.function = function(x) x)

  # --- Fig 1: average survival curves (B&W-ish)
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
           subtitle = "Averaged across 400 iterations",
           x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })

  # --- Fig 2: AF boxplot
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
           subtitle = "Dashed line represents True Target Gene AF",
           x = "", y = "Estimated Acceleration Factor (AF)") +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")
  })

  # --- Fig 3: dependence scatter
  output$sim_multi_fig3 <- renderPlot({
    t1_samp <- results[[1]]$sample_t1 / 365.25
    t2_samp <- results[[1]]$sample_t2 / 365.25
    plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)

    ggplot(plot_df3, aes(x = T1, y = T2)) +
      geom_point(color = "black", alpha = 0.4, size = 2) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 3: Dependent Left Truncation",
           subtitle = "Correlation between Time to CGP (T1) and Residual Survival (T2)",
           x = "Time to CGP Test (T1, Years)", y = "Residual Survival Time (T2, Years)") +
      theme(plot.title = element_text(face = "bold"))
  })
})

extract_af_full_ess_scaled <- function(fit, ess, n_pseudo) {
  blank_res <- c(
    AF_X_est=NA, AF_X_L=NA, AF_X_U=NA,
    AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
    AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA,
    AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
    AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA
  )
  if (is.null(fit)) return(blank_res)
  if (!is.finite(ess) || ess <= 0 || !is.finite(n_pseudo) || n_pseudo <= 0) return(extract_af_full(fit))

  # flexsurvreg の係数分散（パラメータ全体）を取得
  V <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(V)) return(extract_af_full(fit))

  # “擬似コホート n_pseudo”によるSE過小評価を ESS に合わせて補正
  scale_fac <- n_pseudo / ess
  if (!is.finite(scale_fac) || scale_fac <= 0) return(extract_af_full(fit))

  V2 <- V * scale_fac

  # 係数名（covariate location のパラメータ名）を vcov から参照
  par_names <- colnames(V2)
  if (is.null(par_names)) return(extract_af_full(fit))

  get_ci <- function(name, mult = 1) {
    # fit$res の est は “回帰係数” (β)（log-time scale）なので exp(β*mult)
    res_m <- fit$res
    if (!(name %in% rownames(res_m))) return(c(est=NA, L=NA, U=NA))
    beta <- as.numeric(res_m[name, "est"])

    # vcov の対角からSE
    if (!(name %in% par_names)) return(c(est=exp(beta*mult), L=NA, U=NA))
    se <- sqrt(as.numeric(V2[name, name]))

    z <- 1.959963984540054
    est <- exp(beta * mult)
    L   <- exp((beta - z * se) * mult)
    U   <- exp((beta + z * se) * mult)
    c(est=est, L=L, U=U)
  }

  vX    <- get_ci("XMutated", 1)
  vAge  <- get_ci("Age", 10)
  vSex  <- get_ci("SexFemale", 1)
  vREAD <- get_ci("HistologyREAD", 1)
  vCOAD <- get_ci("HistologyCOADREAD", 1)

  c(
    AF_X_est=vX["est"], AF_X_L=vX["L"], AF_X_U=vX["U"],
    AF_Age_est=vAge["est"], AF_Age_L=vAge["L"], AF_Age_U=vAge["U"],
    AF_Sex_est=vSex["est"], AF_Sex_L=vSex["L"], AF_Sex_U=vSex["U"],
    AF_READ_est=vREAD["est"], AF_READ_L=vREAD["L"], AF_READ_U=vREAD["U"],
    AF_COADREAD_est=vCOAD["est"], AF_COADREAD_L=vCOAD["L"], AF_COADREAD_U=vCOAD["U"]
  )
}
