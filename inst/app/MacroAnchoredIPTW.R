# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model)
#  - Proposed: Age-stratified external OS curve matching via weights on CGP cohort
#  - External info: 1–5y OS by Age only (5 points); internally smooth by llogis
#  - DGP: T1 depends on X/Histology as before (keep your current world)
#  - Assumption: X and Histology marginal frequencies match between CGP & population
#  - Goal: estimate conditional AF (AFT regression coefficients) for OS
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(flexsurv)
  library(survival)
  library(ggplot2)
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
# External survival: Fit llogis to 1–5y survival points (Age-specific)
#   Input: S_y (length 5) in %; Output: list(shape, scale) for rllogis/Surv
# -------------------------------------------------------------------------
fit_llogis_from_S15 <- function(S_y_percent) {
  # t in years: 1..5
  t_years <- 1:5
  S <- pmax(pmin(as.numeric(S_y_percent) / 100, 0.999), 0.001)

  # llogis survival: S(t) = 1 / (1 + (t/scale)^shape)
  # optimize on log-scale and log-shape for stability
  obj <- function(par) {
    log_shape <- par[1]
    log_scale <- par[2]
    shape <- exp(log_shape)
    scale <- exp(log_scale)
    S_hat <- 1 / (1 + (t_years / scale)^shape)
    sum((S_hat - S)^2)
  }

  # rough init: scale ~ median survival time around where S=0.5
  init <- c(log(1.5), log(3))
  opt <- tryCatch(
    optim(init, obj, method = "Nelder-Mead", control = list(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(shape = 1.5, scale = 3.0))
  }
  list(shape = exp(opt$par[1]), scale = exp(opt$par[2]))
}

llogis_surv <- function(t_years, shape, scale) {
  t_years <- pmax(t_years, 1e-8)
  1 / (1 + (t_years / scale)^shape)
}

# -------------------------------------------------------------------------
# Step 1: Direct T_obs Curve Matching (Age-stratified)
#  - Target curve is built ONLY from macro 1–5y survival points (Age-specific),
#    then smoothed to 0.5y grid by llogis.
# -------------------------------------------------------------------------
calculate_obs_matching_weights <- function(Data_cgp, macro_S15_by_ageclass) {
  stopifnot("Age_class" %in% names(Data_cgp))
  stopifnot("T_obs" %in% names(Data_cgp))
  stopifnot("Event" %in% names(Data_cgp))

  Data_cgp$iptw <- 1.0

  # 0.5-year grid up to 5 years (10 points)
  t_points_years <- seq(0.5, 5.0, by = 0.5)
  t_points_days  <- t_points_years * 365.25

  # per Age_class optimization
  for (ag in unique(Data_cgp$Age_class)) {
    if (!(ag %in% names(macro_S15_by_ageclass))) next

    idx <- which(Data_cgp$Age_class == ag)
    if (length(idx) < 50) next
    ag_data <- Data_cgp[idx, , drop = FALSE]

    # ---- build target curve from external 1–5y points (smoothed by llogis) ----
    pars <- fit_llogis_from_S15(macro_S15_by_ageclass[[ag]])
    S_target <- llogis_surv(t_points_years, pars$shape, pars$scale)
    S_target[!is.finite(S_target)] <- 0
    S_target <- pmax(pmin(S_target, 0.999), 0.001)

    # ---- smooth weight function of T_obs (polynomial of log time) ----
    log_t <- log(pmax(ag_data$T_obs / 365.25, 1e-6))
    log_t2 <- log_t^2
    log_t3 <- log_t^3

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t + theta[2] * log_t2 + theta[3] * log_t3)
      w <- pmax(w, 1e-6)
      w <- w / mean(w)

      km <- tryCatch(
        survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w),
        error = function(e) NULL
      )
      if (is.null(km)) return(1e9)

      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      S_est[!is.finite(S_est)] <- 0
      S_est <- pmax(pmin(S_est, 0.999), 0.001)

      # main alignment loss
      loss <- sum((S_est - S_target)^2)

      # regularization to avoid overfit / explosion
      loss + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(
      optim(c(0, 0, 0), obj_func, method = "Nelder-Mead", control = list(maxit = 2000)),
      error = function(e) list(par = c(0, 0, 0))
    )
    th <- opt$par

    w_opt <- exp(th[1] * log_t + th[2] * log_t2 + th[3] * log_t3)
    w_opt <- pmax(w_opt, 1e-6)

    # winsorize (avoid explosion)
    lb <- quantile(w_opt, 0.01, na.rm = TRUE)
    ub <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lb, pmin(w_opt, ub))

    Data_cgp$iptw[idx] <- w_opt / mean(w_opt)
  }

  # normalize weights
  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)

  # ESS for reporting (and optional CI interpretation)
  ESS <- (sum(Data_cgp$iptw)^2) / sum(Data_cgp$iptw^2)
  if (!is.finite(ESS) || is.na(ESS) || ESS <= 1) ESS <- nrow(Data_cgp)

  list(data = Data_cgp, ESS = ESS)
}

# -------------------------------------------------------------------------
# Extractors for AFT Models
# -------------------------------------------------------------------------
extract_af_full <- function(fit) {
  blank <- c(
    AF_X_est=NA, AF_X_L=NA, AF_X_U=NA,
    AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
    AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA,
    AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
    AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA
  )
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
# Safe sim times for naive/LT curve comparison (llogis only)
# -------------------------------------------------------------------------
gen_sim_times_llogis_safe <- function(fit, newdata) {
  if (is.null(fit) || is.null(newdata) || nrow(newdata) == 0) return(rep(NA_real_, nrow(newdata)))

  res_m <- fit$res
  n <- nrow(newdata)
  xb <- rep(0, n)
  ind <- function(x) { v <- as.numeric(x); v[is.na(v)] <- 0; v }

  if ("XMutated" %in% rownames(res_m) && "X" %in% names(newdata)) xb <- xb + res_m["XMutated","est"] * ind(newdata$X == "Mutated")
  if ("Age" %in% rownames(res_m) && "Age" %in% names(newdata)) xb <- xb + res_m["Age","est"] * ind(newdata$Age)
  if ("SexFemale" %in% rownames(res_m) && "Sex" %in% names(newdata)) xb <- xb + res_m["SexFemale","est"] * ind(newdata$Sex == "Female")
  if ("HistologyREAD" %in% rownames(res_m) && "Histology" %in% names(newdata)) xb <- xb + res_m["HistologyREAD","est"] * ind(newdata$Histology == "READ")
  if ("HistologyCOADREAD" %in% rownames(res_m) && "Histology" %in% names(newdata)) xb <- xb + res_m["HistologyCOADREAD","est"] * ind(newdata$Histology == "COADREAD")

  sim_t <- rep(NA_real_, n)
  if (all(c("shape","scale") %in% rownames(res_m))) {
    shape <- as.numeric(res_m["shape","est"])
    sc0   <- as.numeric(res_m["scale","est"])
    if (is.finite(shape) && shape > 0 && is.finite(sc0) && sc0 > 0) {
      mu <- log(sc0) + xb
      mu <- pmin(pmax(mu, log(0.01)), log(365.25 * 100))
      sim_t <- flexsurv::rllogis(n, shape = shape, scale = exp(mu))
    }
  }
  sim_t[!is.finite(sim_t)] <- NA_real_
  if (all(is.na(sim_t))) sim_t <- rep(365.25, n)
  sim_t[is.na(sim_t)] <- median(sim_t, na.rm = TRUE)
  pmin(sim_t, 365.25 * 100)
}

# -------------------------------------------------------------------------
# Metrics calculation
# -------------------------------------------------------------------------
calc_marginal_metrics <- function(sim_times_days) {
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(c(Med=NA_real_, S5=NA_real_))
  c(
    Med = median(sim_times_days, na.rm = TRUE) / 365.25,
    S5  = mean(sim_times_days > 5 * 365.25, na.rm = TRUE) * 100
  )
}

calc_curve_km <- function(sim_times_days) {
  t_seq <- seq(0, 365.25 * 5, length.out = 100)
  if (is.null(sim_times_days) || all(is.na(sim_times_days))) return(rep(NA_real_, length(t_seq)))
  km <- tryCatch(survfit(Surv(sim_times_days, rep(1, length(sim_times_days))) ~ 1), error = function(e) NULL)
  if (is.null(km)) return(rep(NA_real_, length(t_seq)))
  approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
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

  # --- DGP: generate TRUE OS time T (depends on covariates incl X) ---
  u <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true_base <- (True_Med * 365.25) * ((1 - u) / u)^(1 / True_Shape)

  # True OS time with covariate effects (AFT-style)
  T_true <- T_true_base * AF_total

  # --- Independent truncation world: generate T1 independent of T ---
  # IMPORTANT: T1 must NOT depend on T_true_base or AF_total
  if (is.null(t1_pat)) t1_pat <- "indep"

  min_days <- 30
  jitter_days <- 7

  # choose a fixed upper bound for entry-time distribution (independent of T)
  # (you can tune; must not reference T_true)
  max_entry_days <- 365.25 * 5   # e.g., allow CGP up to 5 years after diagnosis

  if (t1_pat == "early") {
    # earlier entry times (still independent of T)
    T1_base <- runif(N_macro, min_days, min_days + 0.3 * (max_entry_days - min_days))

  } else if (t1_pat %in% c("dep_1yr", "real")) {
    # later-biased (your "real" option) but independent of T
    a <- 0.5; b <- 1
    frac <- rbeta(N_macro, a, b)
    T1_base <- min_days + frac * (max_entry_days - min_days)

  } else if (t1_pat %in% c("dep_2yr", "rev")) {
    # later-biased (your "rev" option) but independent of T
    a <- 1; b <- 0.5
    frac <- rbeta(N_macro, a, b)
    T1_base <- min_days + frac * (max_entry_days - min_days)

  } else {
    # quasi-uniform entry (independent of T)
    T1_base <- runif(N_macro, min_days, max_entry_days)
  }

  # avoid spike at exactly min_days
  T1_base <- ifelse(T1_base <= min_days, min_days + runif(N_macro, 0, jitter_days), T1_base)

  # observed entry time
  T1 <- T1_base

  # residual survival after entry (deterministic decomposition; independence is ensured by construction above)
  T2_true <- pmax(0.1, T_true - T1)

  # Age class (external survival strata): <40, 40s, 50s, 60s, 70s, 80+
  Age_class_macro <- cut(
    Age,
    breaks = c(-Inf, 39, 49, 59, 69, 79, Inf),
    labels = c("40未満", "40代", "50代", "60代", "70代", "80以上"),
    right = TRUE
  )

  # External info (fixed): only 1–5y survival by age strata (5 points)
  # Use JSON: Data_age_survival_5_year[["大腸がん(大腸癌)"]][["40代"]] etc.
  macro_S15_by_ageclass <- Data_age_survival_5_year[["大腸がん(大腸癌)"]]

  # Safety: ensure required keys exist
  required_age_keys <- c("40未満","40代","50代","60代","70代","80以上")
  macro_S15_by_ageclass <- macro_S15_by_ageclass[required_age_keys]

  # --- Select CGP cohort (as before) ---
  cgp_idx <- which(T_true > T1 & T1 <= 365.25 * 10)
  if (length(cgp_idx) < N_target) return(NULL)

  age_cgp <- Age[cgp_idx]
  prob_select <- plogis((60 - age_cgp) / 7)   # 60歳を境に滑らかに
  prob_select <- 0.35 + 0.30 * prob_select    # -> おおむね [0.35, 0.65]
  sel <- sample(cgp_idx, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target,
    Age = Age[sel],
    Age_class = {
      age_tmp <- Age[sel]
      out <- cut(
        age_tmp,
        breaks = c(-Inf, 39, 49, 59, 69, 79, Inf),
        labels = c("40未満", "40代", "50代", "60代", "70代", "80以上"),
        right = TRUE
      )
      out[is.na(age_tmp)] <- "全年齢"
      factor(out, levels = c("40未満","40代","50代","60代","70代","80以上","全年齢"))
    },
    Sex = factor(Sex[sel], levels = c("Male", "Female")),
    Histology = factor(Histology[sel], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[sel] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[sel],
    T1 = T1[sel]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # --- Censoring on T2 (as before) ---
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

  # --- Proposed weights: Age-stratified matching to external curve (1–5y points) ---
  mw <- calculate_obs_matching_weights(Data_cgp, macro_S15_by_ageclass)
  Data_cgp <- mw$data
  ESS <- mw$ESS

  # --- Covariates for regression AF (conditional AF) ---
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")
  form_base <- paste(valid_covs, collapse = " + ")

  # --- Comparators ---
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", form_base))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", form_base))

  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt    <- tryCatch(flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  # --- Proposed: weighted AFT on observed (T_obs, Event) ---
  fit_prop <- tryCatch(
    flexsurvreg(form_naive, data = Data_cgp, weights = Data_cgp$iptw, dist = "llogis"),
    error = function(e) NULL
  )
  if (is.null(fit_prop)) return(NULL)

  # --- Curves / metrics ---
  # True curve is from full macro T_true (for simulation evaluation)
  curve_true <- calc_curve_km(T_true)
  True_Metrics <- calc_marginal_metrics(T_true)

  # Naive/LT: for curve comparison, use model-implied simulated times
  sim_naive <- if (!is.null(fit_naive)) gen_sim_times_llogis_safe(fit_naive, Data_cgp) else rep(NA_real_, nrow(Data_cgp))
  sim_lt    <- if (!is.null(fit_lt))    gen_sim_times_llogis_safe(fit_lt,    Data_cgp) else rep(NA_real_, nrow(Data_cgp))

  curve_naive <- calc_curve_km(sim_naive)
  curve_lt    <- calc_curve_km(sim_lt)

  # Proposed curve: weighted KM on observed data (this is what you are matching)
  t_seq_days <- seq(0, 365.25 * 5, length.out = 100)
  km_prop <- tryCatch(
    survfit(Surv(T_obs, Event) ~ 1, data = Data_cgp, weights = Data_cgp$iptw),
    error = function(e) NULL
  )
  if (!is.null(km_prop)) {
    curve_prop <- approx(km_prop$time, km_prop$surv, xout = t_seq_days, method = "constant", f = 0, rule = 2)$y
    curve_prop[!is.finite(curve_prop)] <- NA_real_
  } else {
    curve_prop <- rep(NA_real_, length(t_seq_days))
  }

  # Proposed "marginal" metrics from weighted KM (not model-implied)
  # median from KM object (if available)
  prop_med <- tryCatch({
    q <- quantile(km_prop, probs = 0.5)$quantile
    as.numeric(q) / 365.25
  }, error = function(e) NA_real_)

  prop_s5 <- tryCatch({
    s5 <- approx(km_prop$time, km_prop$surv, xout = 5 * 365.25, method = "constant", f = 0, rule = 2)$y
    as.numeric(s5) * 100
  }, error = function(e) NA_real_)

  Prop_Metrics <- c(Med = unname(prop_med), S5 = unname(prop_s5))

  list(
    ESS = ESS,
    True_Metrics = True_Metrics,
    Naive_AF = extract_af_full(fit_naive),
    Naive_Metrics = calc_marginal_metrics(sim_naive),
    LT_AF = extract_af_full(fit_lt),
    LT_Metrics = calc_marginal_metrics(sim_lt),
    Prop_AF = extract_af_full(fit_prop),
    Prop_Metrics = Prop_Metrics,
    curve_true = curve_true,
    curve_naive = curve_naive,
    curve_lt = curve_lt,
    curve_prop = curve_prop,
    sample_t1 = Data_cgp$T1[1:min(500, nrow(Data_cgp))],
    sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))]
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

  shiny::validate(shiny::need(!is.null(res), "Simulation failed (model fitting / generation). Try changing settings."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c(
        "Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
        "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"
      ),
      `True` = c(
        safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
        safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"
      ),
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
      `Proposed (Weighted AFT)` = c(
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
      Model = factor(rep(
        c("1. True Marginal (Macro)",
          "2. Naive AFT (No LT)",
          "3. Standard LT AFT (Non-inf LT)",
          "4. Proposed (Age-calibrated weights)"),
        each = 100
      ))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(
        title = "Tamura & Ikegami Simulation",
        subtitle = "Single Run: External Age-based Calibration + Weighted AFT",
        x = "Time from Diagnosis (Years)",
        y = "Overall Survival Probability"
      ) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank()
      )
  })
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)

  n_sims <- 100
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
      Metric = c(
        "Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
        "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"
      ),
      `True Value` = c(
        safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
        sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"
      ),
      `Naive AFT` = c(
        calc_stats_af("Naive_AF", "AF_X", True_AF),
        calc_stats_af("Naive_AF", "AF_Age", 0.85),
        calc_stats_af("Naive_AF", "AF_Sex", 1.10),
        calc_stats_af("Naive_AF", "AF_READ", 0.90),
        calc_stats_af("Naive_AF", "AF_COADREAD", 0.95),
        calc_stats_metrics("Naive_Metrics", "Med", mean_true_med),
        calc_stats_metrics("Naive_Metrics", "S5",  mean_true_s5),
        "N/A"
      ),
      `Standard LT AFT` = c(
        calc_stats_af("LT_AF", "AF_X", True_AF),
        calc_stats_af("LT_AF", "AF_Age", 0.85),
        calc_stats_af("LT_AF", "AF_Sex", 1.10),
        calc_stats_af("LT_AF", "AF_READ", 0.90),
        calc_stats_af("LT_AF", "AF_COADREAD", 0.95),
        calc_stats_metrics("LT_Metrics", "Med", mean_true_med),
        calc_stats_metrics("LT_Metrics", "S5",  mean_true_s5),
        "N/A"
      ),
      `Proposed (Weighted AFT)` = c(
        paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"),
        paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"),
        paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"),
        paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5",  mean_true_s5), "</b>"),
        paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>")
      )
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
      Model = factor(rep(
        c("1. True Marginal", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed (Age-calibrated)"),
        each = 100
      ))
    )

    title_txt <- paste(strwrap("Figure 1: Reconstructed Marginal Survival Curves", width = 38), collapse = "\n")
    sub_txt   <- paste(strwrap("Averaged across iterations", width = 42), collapse = "\n")

    ggplot(plot_df1, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(
        title = title_txt,
        subtitle = sub_txt,
        x = "Time from Diagnosis (Years)",
        y = "Overall Survival Probability"
      ) +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(face = "bold", size = 16, lineheight = 1.05),
        plot.subtitle = element_text(size = 12, lineheight = 1.05),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank()
      )
  })

  # Fig2: boxplot of AF_X
  output$sim_multi_fig2 <- renderPlot({
    af_naive <- pull_vals(function(x) x$Naive_AF["AF_X_est"])
    af_lt    <- pull_vals(function(x) x$LT_AF["AF_X_est"])
    af_prop  <- pull_vals(function(x) x$Prop_AF["AF_X_est"])

    plot_df2 <- data.frame(
      AF = c(af_naive, af_lt, af_prop),
      Model = factor(rep(
        c("Naive AFT", "Standard LT AFT", "Proposed (Weighted AFT)"),
        c(length(af_naive), length(af_lt), length(af_prop))
      ), levels = c("Naive AFT", "Standard LT AFT", "Proposed (Weighted AFT)"))
    )

    title_txt <- paste(strwrap("Figure 2: Distribution of Estimated AF (Dashed line = True AF)", width = 38), collapse = "\n")

    ggplot(plot_df2, aes(x = Model, y = AF, fill = Model)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 1) +
      geom_hline(yintercept = True_AF, color = "black", linetype = "dashed", linewidth = 1.2) +
      coord_cartesian(ylim = c(0, max(True_AF * 2.5, 3))) +
      theme_minimal(base_size = 15)+
      labs(title = title_txt, x = "", y = "Estimated AF") +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(face="bold", size=16, lineheight=1.05),
        plot.margin = margin(t=10, r=10, b=10, l=10),
        axis.text.x = element_text(angle = 15, hjust = 1)
      )
  })

  # Fig3: dependence scatter (T1 vs T2_true)
  output$sim_multi_fig3 <- renderPlot({
    t1_samp <- results[[1]]$sample_t1 / 365.25
    t2_samp <- results[[1]]$sample_t2 / 365.25
    plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)

    title_txt <- paste(strwrap("Figure 3: Dependent Truncation Structure (Correlation between T1 and residual survival)", width = 38), collapse = "\n")

    ggplot(plot_df3, aes(x = T1, y = T2)) +
      geom_point(color = "black", alpha = 0.4, size = 2) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE) +
      theme_minimal(base_size = 15) +
      labs(title = title_txt, x = "Time to CGP Test (T1, Years)", y = "Residual Survival Time (T2, Years)") +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(face="bold", size=16, lineheight=1.05),
        plot.margin = margin(t=10, r=10, b=10, l=10)
      )
  })
})
