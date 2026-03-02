# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.3.2)
#   + Proposed (Total Effect) uses Bootstrap CI (B=200, percentile)
#   + Refactored for Shiny Stability (Reactive Values & Error Handling)
# =========================================================================

library(shiny)
library(dplyr)
library(flexsurv)
library(survival)
library(ggplot2)
library(splines)
library(scales)

# =========================================================
# Helper Functions
# =========================================================

# Find censoring parameter
find_censoring_param <- function(T2, target_rate, pattern) {
  if (pattern == "peak1y") {
    obj_func <- function(k) {
      scale_val <- 365.25 / ((k - 1) / k)^(1 / k)
      (sum(rweibull(length(T2), shape = k, scale = scale_val) < T2) / length(T2)) - target_rate
    }
    opt_k <- tryCatch(uniroot(obj_func, interval = c(1.05, 10))$root, error = function(e) 2.0)
    list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE)

  } else if (pattern == "ushape") {
    u_fixed <- runif(length(T2))
    obj_func <- function(scale_param) {
      C2 <- ifelse(u_fixed < 0.5,
                   rweibull(length(T2), shape = 0.5, scale = scale_param),
                   rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func, interval = c(0.1, 1e6))$root, error = function(e) median(T2) * 1.5)
    list(shape = NA, scale = opt_scale, is_ushape = TRUE)

  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func_scale, interval = c(0.1, 1e6))$root, error = function(e) median(T2) * 1.5)
    list(shape = if (pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE)
  }
}

# Generate simulated times
gen_sim_times <- function(fit, newdata, dist_type = NULL) {
  if (is.null(fit)) stop("fit is NULL")
  if (is.null(dist_type)) dist_type <- fit$dlist$name
  n <- nrow(newdata)
  if (n == 0) return(numeric(0))

  # Try predict(type="lp")
  mu <- tryCatch(predict(fit, newdata = newdata, type = "lp"), error = function(e) NULL)

  # Fallback: manual LP (location)
  if (is.null(mu)) {
    res_m <- fit$res

    if (dist_type %in% c("gengamma", "gen_gamma", "generalized gamma")) {
      if (!("mu" %in% rownames(res_m))) stop("gengamma: baseline mu not found")
      base_mu <- as.numeric(res_m["mu", "est"])
    } else {
      if (!("scale" %in% rownames(res_m))) stop("baseline scale not found")
      base_mu <- log(as.numeric(res_m["scale", "est"]))
    }

    lp <- rep(0, n)

    if ("XMutated" %in% rownames(res_m) && "X" %in% names(newdata)) {
      lp <- lp + as.numeric(res_m["XMutated", "est"]) * as.numeric(newdata$X == "Mutated")
    }
    if ("Age" %in% rownames(res_m) && "Age" %in% names(newdata)) {
      lp <- lp + as.numeric(res_m["Age", "est"]) * as.numeric(newdata$Age)
    }
    if ("SexFemale" %in% rownames(res_m) && "Sex" %in% names(newdata)) {
      lp <- lp + as.numeric(res_m["SexFemale", "est"]) * as.numeric(newdata$Sex == "Female")
    }
    if ("HistologyREAD" %in% rownames(res_m) && "Histology" %in% names(newdata)) {
      lp <- lp + as.numeric(res_m["HistologyREAD", "est"]) * as.numeric(newdata$Histology == "READ")
    }
    if ("HistologyCOADREAD" %in% rownames(res_m) && "Histology" %in% names(newdata)) {
      lp <- lp + as.numeric(res_m["HistologyCOADREAD", "est"]) * as.numeric(newdata$Histology == "COADREAD")
    }

    ns_rows <- grep("^ns[0-9]+$", rownames(res_m), value = TRUE)
    for (v in ns_rows) {
      if (v %in% names(newdata)) lp <- lp + as.numeric(res_m[v, "est"]) * as.numeric(newdata[[v]])
    }

    mu <- base_mu + lp
  }

  # RNG by dist
  if (dist_type %in% c("weibull", "weibullPH", "weibullph")) {
    shape <- as.numeric(fit$res["shape", "est"])
    scale <- exp(mu)
    return(rweibull(n, shape = shape, scale = scale))
  }
  if (dist_type %in% c("llogis", "loglogistic")) {
    scl <- as.numeric(fit$res["scale", "est"])
    scale_param <- exp(mu)
    shape_param <- 1 / scl
    return(flexsurv::rllogis(n, shape = shape_param, scale = scale_param))
  }
  if (dist_type %in% c("gengamma", "gen_gamma", "generalized gamma")) {
    sigma <- as.numeric(fit$res["sigma", "est"])
    Q     <- as.numeric(fit$res["Q", "est"])
    return(flexsurv::rgengamma(n, mu = mu, sigma = sigma, Q = Q))
  }
  stop(paste0("Unsupported dist_type: ", dist_type))
}

# Calculate Calibrated IPTW
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
      km <- tryCatch(survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w),
                     error = function(e) NULL)
      if (is.null(km)) return(1e6)
      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(optim(c(0, 0), obj_func, method = "BFGS"),
                    error = function(e) list(par = c(0, 0)))
    theta_hat <- opt$par
    w_opt <- pmax(exp(theta_hat[1] * log_t1 + theta_hat[2] * log_t1_sq), 1e-4)

    lower_bound <- quantile(w_opt, 0.01, na.rm = TRUE)
    upper_bound <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lower_bound, pmin(w_opt, upper_bound))
    data$iptw[ag_data_idx] <- w_opt / mean(w_opt)
  }
  data
}

# Extract AF
extract_af_full <- function(fit) {
  blank_res <- c(AF_X_est=NA, AF_X_L=NA, AF_X_U=NA,
                 AF_Age_est=NA, AF_Age_L=NA, AF_Age_U=NA,
                 AF_Sex_est=NA, AF_Sex_L=NA, AF_Sex_U=NA,
                 AF_READ_est=NA, AF_READ_L=NA, AF_READ_U=NA,
                 AF_COADREAD_est=NA, AF_COADREAD_L=NA, AF_COADREAD_U=NA)
  if (is.null(fit)) return(blank_res)
  res_m <- fit$res

  get_est <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "est"] * mult) else NA_real_
  get_L <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "L95%"] * mult) else NA_real_
  get_U <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "U95%"] * mult) else NA_real_

  c(AF_X_est = get_est("XMutated"), AF_X_L = get_L("XMutated"), AF_X_U = get_U("XMutated"),
    AF_Age_est = get_est("Age", 10), AF_Age_L = get_L("Age", 10), AF_Age_U = get_U("Age", 10),
    AF_Sex_est = get_est("SexFemale"), AF_Sex_L = get_L("SexFemale"), AF_Sex_U = get_U("SexFemale"),
    AF_READ_est = get_est("HistologyREAD"), AF_READ_L = get_L("HistologyREAD"), AF_READ_U = get_U("HistologyREAD"),
    AF_COADREAD_est = get_est("HistologyCOADREAD"), AF_COADREAD_L = get_L("HistologyCOADREAD"), AF_COADREAD_U = get_U("HistologyCOADREAD"))
}

# Bootstrap CI for Proposed
bootstrap_prop_ci <- function(Data_cgp, ref_surv_list, valid_covs, B = 200, n_pseudo = 10000, ns_df = 3) {
  n <- nrow(Data_cgp)
  lev_Sex <- levels(Data_cgp$Sex)
  lev_His <- levels(Data_cgp$Histology)
  lev_X   <- levels(Data_cgp$X)

  store <- list(AF_X=numeric(0), AF_Age=numeric(0), AF_Sex=numeric(0), AF_READ=numeric(0), AF_COADREAD=numeric(0))

  for (b in 1:B) {
    idx <- sample.int(n, size = n, replace = TRUE)
    d_b <- Data_cgp[idx, , drop = FALSE]

    d_b$Sex <- factor(d_b$Sex, levels = lev_Sex)
    d_b$Histology <- factor(d_b$Histology, levels = lev_His)
    d_b$X <- factor(d_b$X, levels = lev_X)

    if (sum(d_b$Event) < 20) next

    d_b <- tryCatch(calculate_calibrated_iptw(d_b, ref_surv_list), error = function(e) NULL)
    if (is.null(d_b)) next
    d_b$iptw <- d_b$iptw / mean(d_b$iptw)

    d_b$logT1_scale <- log(pmax(d_b$T1 / 365.25, 1e-6))
    ns_obj_b <- tryCatch(ns(d_b$logT1_scale, df = ns_df), error = function(e) NULL)

    ns_terms_b <- character(0)
    if (!is.null(ns_obj_b)) {
      ns_mat_b <- as.matrix(ns_obj_b)
      colnames(ns_mat_b) <- paste0("ns", seq_len(ncol(ns_mat_b)))
      for (j in seq_len(ncol(ns_mat_b))) d_b[[colnames(ns_mat_b)[j]]] <- ns_mat_b[, j]
      ns_terms_b <- colnames(ns_mat_b)
    }

    d_b$T1_event <- 1
    form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
    form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, ns_terms_b), collapse = " + ")))

    fit_t1_b <- tryCatch(flexsurvreg(form_t1, data = d_b, weights = d_b$iptw, dist = "weibull"), error = function(e) NULL)
    fit_t2_b <- tryCatch(flexsurvreg(form_t2, data = d_b, weights = d_b$iptw, dist = "gengamma"), error = function(e) NULL)
    if (is.null(fit_t1_b) || is.null(fit_t2_b)) next

    idx_p <- sample(seq_len(nrow(d_b)), size = n_pseudo, replace = TRUE, prob = d_b$iptw)
    p_b <- d_b[idx_p, , drop = FALSE]

    p_b$sim_T1 <- tryCatch(gen_sim_times(fit_t1_b, p_b, dist_type = "weibull"), error = function(e) NA_real_)
    if (all(is.na(p_b$sim_T1))) next

    if (!is.null(ns_obj_b)) {
      p_b$logT1_sim_scale <- log(pmax(p_b$sim_T1 / 365.25, 1e-6))
      ns_sim <- tryCatch(predict(ns_obj_b, p_b$logT1_sim_scale), error = function(e) NULL)
      if (!is.null(ns_sim)) {
        colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
        for (j in seq_len(ncol(ns_sim))) p_b[[colnames(ns_sim)[j]]] <- ns_sim[, j]
      } else {
        for (term in ns_terms_b) p_b[[term]] <- 0
      }
    }

    p_b$sim_T2 <- tryCatch(gen_sim_times(fit_t2_b, p_b, dist_type = "gengamma"), error = function(e) NA_real_)
    if (all(is.na(p_b$sim_T2))) next

    p_b$sim_OS <- p_b$sim_T1 + p_b$sim_T2
    p_b$sim_Event <- 1

    form_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
    fit_final_b <- tryCatch(flexsurvreg(form_final, data = p_b, dist = "llogis"), error = function(e) NULL)
    if (is.null(fit_final_b)) next

    af <- extract_af_full(fit_final_b)

    if (!is.na(af["AF_X_est"])) store$AF_X <- c(store$AF_X, af["AF_X_est"])
    if (!is.na(af["AF_Age_est"])) store$AF_Age <- c(store$AF_Age, af["AF_Age_est"])
    if (!is.na(af["AF_Sex_est"])) store$AF_Sex <- c(store$AF_Sex, af["AF_Sex_est"])
    if (!is.na(af["AF_READ_est"])) store$AF_READ <- c(store$AF_READ, af["AF_READ_est"])
    if (!is.na(af["AF_COADREAD_est"])) store$AF_COADREAD <- c(store$AF_COADREAD, af["AF_COADREAD_est"])
  }

  q025_975 <- function(x) {
    if (length(x) < 30) return(c(L=NA_real_, U=NA_real_))
    as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, type = 6))
  }

  list(
    AF_X_L = q025_975(store$AF_X)[1], AF_X_U = q025_975(store$AF_X)[2],
    AF_Age_L = q025_975(store$AF_Age)[1], AF_Age_U = q025_975(store$AF_Age)[2],
    AF_Sex_L = q025_975(store$AF_Sex)[1], AF_Sex_U = q025_975(store$AF_Sex)[2],
    AF_READ_L = q025_975(store$AF_READ)[1], AF_READ_U = q025_975(store$AF_READ)[2],
    AF_COADREAD_L = q025_975(store$AF_COADREAD)[1], AF_COADREAD_U = q025_975(store$AF_COADREAD)[2],
    n_ok = min(length(store$AF_X), length(store$AF_Age))
  )
}

# Core Simulation Engine
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, True_Med, True_Shape,
                              t1_pat, cens_pat, cens_rate, return_data = FALSE) {

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
    ID = 1:N_target, Age = Age[selected_indices],
    Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
    Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[selected_indices],
    T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (isTRUE(params$is_ushape)) {
    u_c <- runif(N_target)
    C2 <- ifelse(u_c < 0.5,
                 rweibull(N_target, shape = 0.5, scale = params$scale),
                 rweibull(N_target, shape = 3.0, scale = params$scale * 2))
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

  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100)
  }
  Data_cgp <- tryCatch(calculate_calibrated_iptw(Data_cgp, ref_surv_list), error = function(e) Data_cgp)
  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)
  ESS <- (sum(Data_cgp$iptw)^2) / sum(Data_cgp$iptw^2)

  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # Safe Spline Generation
  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  ns_obj <- tryCatch(ns(Data_cgp$logT1_scale, df = 3), error = function(e) NULL)
  ns_terms <- character(0)

  if (!is.null(ns_obj)) {
    ns_mat <- as.matrix(ns_obj)
    colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
    for (j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]
    ns_terms <- colnames(ns_mat)
  }

  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt    <- tryCatch(flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurvreg(form_t1, data = Data_cgp, weights = Data_cgp$iptw, dist = "weibull"), error = function(e) NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, ns_terms), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurvreg(form_t2, data = Data_cgp, weights = Data_cgp$iptw, dist = "gengamma"), error = function(e) NULL)

  if (is.null(fit_t1) || is.null(fit_t2)) return(NULL)

  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = 10000, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_cgp <- Data_cgp[idx_pseudo, , drop = FALSE]

  pseudo_cgp$sim_T1 <- gen_sim_times(fit_t1, pseudo_cgp, dist_type = "weibull")

  if (!is.null(ns_obj)) {
    pseudo_cgp$logT1_sim_scale <- log(pmax(pseudo_cgp$sim_T1 / 365.25, 1e-6))
    ns_sim <- tryCatch(predict(ns_obj, pseudo_cgp$logT1_sim_scale), error = function(e) NULL)
    if (!is.null(ns_sim)) {
      colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
      for (j in seq_len(ncol(ns_sim))) pseudo_cgp[[colnames(ns_sim)[j]]] <- ns_sim[, j]
    } else {
      for (term in ns_terms) pseudo_cgp[[term]] <- 0
    }
  }

  pseudo_cgp$sim_T2 <- gen_sim_times(fit_t2, pseudo_cgp, dist_type = "gengamma")
  pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
  pseudo_cgp$sim_Event <- 1

  form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
  fit_prop_final <- tryCatch(flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis"), error = function(e) NULL)

  t_seq <- seq(0, 365.25 * 5, length.out = 100)
  calc_marginal_curve <- function(t_data) {
    if (is.null(t_data) || all(is.na(t_data))) return(rep(NA_real_, length(t_seq)))
    km <- survfit(Surv(t_data, rep(1, length(t_data))) ~ 1)
    approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
  }
  calc_marginal_metrics <- function(sim_times) {
    if (is.null(sim_times) || all(is.na(sim_times))) return(c(Med = NA_real_, S5 = NA_real_))
    c(Med = median(sim_times, na.rm = TRUE) / 365.25, S5 = mean(sim_times > 5 * 365.25, na.rm = TRUE) * 100)
  }

  sim_naive <- if(!is.null(fit_naive)) gen_sim_times(fit_naive, Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt    <- if(!is.null(fit_lt))    gen_sim_times(fit_lt,    Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))

  Prop_AF_point <- extract_af_full(fit_prop_final)

  boot_ci <- bootstrap_prop_ci(Data_cgp = Data_cgp, ref_surv_list = ref_surv_list, valid_covs = valid_covs, B = 200, n_pseudo = 10000, ns_df = 3)

  Prop_AF <- Prop_AF_point
  Prop_AF["AF_X_L"] <- boot_ci$AF_X_L; Prop_AF["AF_X_U"] <- boot_ci$AF_X_U
  Prop_AF["AF_Age_L"] <- boot_ci$AF_Age_L; Prop_AF["AF_Age_U"] <- boot_ci$AF_Age_U
  Prop_AF["AF_Sex_L"] <- boot_ci$AF_Sex_L; Prop_AF["AF_Sex_U"] <- boot_ci$AF_Sex_U
  Prop_AF["AF_READ_L"] <- boot_ci$AF_READ_L; Prop_AF["AF_READ_U"] <- boot_ci$AF_READ_U
  Prop_AF["AF_COADREAD_L"] <- boot_ci$AF_COADREAD_L; Prop_AF["AF_COADREAD_U"] <- boot_ci$AF_COADREAD_U

  list(
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
    sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))],
    Prop_boot_n_ok = boot_ci$n_ok
  )
}

# Formatting helpers
safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) if (is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}
format_af_single <- function(af_vec, var_name) {
  est <- af_vec[paste0(var_name, "_est")]
  L <- af_vec[paste0(var_name, "_L")]
  U <- af_vec[paste0(var_name, "_U")]
  if (is.na(est)) return("N/A")
  if (is.na(L) || is.na(U)) return(sprintf("%.2f (NA-NA)", est))
  sprintf("%.2f (%.2f-%.2f)", est, L, U)
}

# =========================================================================
# Server Logic (Reactive Architecture)
# =========================================================================

# Store simulation results safely outside observers
res_single <- reactiveVal(NULL)
res_multi <- reactiveVal(NULL)
stored_inputs <- reactiveVal(list())

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)

  tryCatch({
    showNotification("Running Single Simulation...", type = "message", duration = 3)
    set.seed(as.integer(Sys.time()))

    # Force casting to numeric to prevent fatal crash
    n_val <- as.numeric(input$sim_n)
    af_val <- as.numeric(input$sim_true_af)
    freq_val <- as.numeric(input$sim_mut_freq) / 100
    med_val <- as.numeric(input$sim_true_med)
    shape_val <- as.numeric(input$sim_true_shape)
    cens_val <- as.numeric(input$sim_cens_rate) / 100

    res <- run_sim_iteration(
      n_val, af_val, freq_val, med_val, shape_val,
      input$sim_t1_pattern, input$sim_cens_pattern, cens_val, return_data = TRUE
    )

    if (is.null(res)) {
      showNotification("Simulation generated too few events to model. Please adjust parameters.", type = "warning", duration = 5)
    } else {
      res_single(res)
      stored_inputs(list(n = n_val, af = af_val))
    }

  }, error = function(e) {
    showNotification(paste("Critical Error:", e$message), type = "error", duration = 10)
  })
})

output$sim_result_table <- renderTable({
  req(res_single(), stored_inputs()$af)
  res <- res_single()
  True_AF <- stored_inputs()$af

  data.frame(
    Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
               "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"),
    `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
               safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"),
    `Naive AFT` = c(format_af_single(res$Naive_AF, "AF_X"), format_af_single(res$Naive_AF, "AF_Age"),
                    format_af_single(res$Naive_AF, "AF_Sex"), format_af_single(res$Naive_AF, "AF_READ"),
                    format_af_single(res$Naive_AF, "AF_COADREAD"), safe_fmt(res$Naive_Metrics["Med"]),
                    safe_fmt(res$Naive_Metrics["S5"], 1), safe_fmt(stored_inputs()$n, 0)),
    `Standard LT AFT` = c(format_af_single(res$LT_AF, "AF_X"), format_af_single(res$LT_AF, "AF_Age"),
                          format_af_single(res$LT_AF, "AF_Sex"), format_af_single(res$LT_AF, "AF_READ"),
                          format_af_single(res$LT_AF, "AF_COADREAD"), safe_fmt(res$LT_Metrics["Med"]),
                          safe_fmt(res$LT_Metrics["S5"], 1), safe_fmt(stored_inputs()$n, 0)),
    `Proposed..Total.Effect.` = c(paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"),
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
  req(res_single())
  res <- res_single()
  t_seq <- seq(0, 5, length.out = 100)

  plot_df <- data.frame(
    Time = rep(t_seq, 4),
    Survival = c(res$curve_true, res$curve_naive, res$curve_lt, res$curve_prop),
    Model = factor(rep(c("1. True Marginal (Macro)", "2. Naive AFT (Immortal Bias)", "3. Standard LT (Dep. Trunc Bias)", "4. Proposed..Total.Effect."), each = 100))
  )

  ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    theme_minimal(base_size = 15) +
    labs(title = "Tamura & Ikegami Model Ver 2.3.2 Validation", subtitle = "Single Run Results",
         x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)

  tryCatch({
    n_sims <- 400
    results <- list()

    n_val <- as.numeric(input$sim_n)
    af_val <- as.numeric(input$sim_true_af)
    freq_val <- as.numeric(input$sim_mut_freq) / 100
    med_val <- as.numeric(input$sim_true_med)
    shape_val <- as.numeric(input$sim_true_shape)
    cens_val <- as.numeric(input$sim_cens_rate) / 100

    withProgress(message = "Running 400 Simulations...", value = 0, {
      for (i in 1:n_sims) {
        set.seed(i)
        res <- run_sim_iteration(
          n_val, af_val, freq_val, med_val, shape_val,
          input$sim_t1_pattern, input$sim_cens_pattern, cens_val, return_data = FALSE
        )
        if (!is.null(res)) results[[length(results) + 1]] <- res
        incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
      }
    })

    if (length(results) == 0) {
      showNotification("All multi-simulations failed.", type = "error")
    } else {
      res_multi(results)
      stored_inputs(list(n = n_val, af = af_val))
    }

  }, error = function(e) {
    showNotification(paste("Critical Error (Multi):", e$message), type = "error", duration = 10)
  })
})

output$sim_multi_result_table <- renderTable({
  req(res_multi(), stored_inputs()$af)
  results <- res_multi()
  True_AF <- stored_inputs()$af

  pull_vals <- function(getter) { v <- sapply(results, getter); v[!is.na(v)] }

  calc_stats_af <- function(model_name, var_name, true_val) {
    est_vals <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_est")])
    l_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_L")])
    u_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_U")])
    if (length(est_vals) == 0) return("N/A")
    mean_est <- mean(est_vals)
    mse_val <- mean((est_vals - true_val)^2)
    cr_val <- mean(l_vals <= true_val & true_val <= u_vals, na.rm=TRUE) * 100
    sprintf("%.2f (MSE: %.3f, CR: %.1f%%)", mean_est, mse_val, cr_val)
  }

  calc_stats_metrics <- function(model_name, metric_name, true_val = NULL) {
    vals <- pull_vals(function(x) x[[model_name]][metric_name])
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    if(!is.null(true_val)){
      mse_val <- mean((vals - true_val)^2)
      return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
    }
    return(sprintf("%.2f", mean_val))
  }

  mean_true_med <- mean(pull_vals(function(x) x$True_Metrics["Med"]), na.rm = TRUE)
  mean_true_s5  <- mean(pull_vals(function(x) x$True_Metrics["S5"]), na.rm = TRUE)

  data.frame(
    Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
               "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),
    `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"),
    `Naive AFT` = c(calc_stats_af("Naive_AF", "AF_X", True_AF), calc_stats_af("Naive_AF", "AF_Age", 0.85),
                    calc_stats_af("Naive_AF", "AF_Sex", 1.10), calc_stats_af("Naive_AF", "AF_READ", 0.90),
                    calc_stats_af("Naive_AF", "AF_COADREAD", 0.95), calc_stats_metrics("Naive_Metrics", "Med", mean_true_med),
                    calc_stats_metrics("Naive_Metrics", "S5", mean_true_s5), "N/A"),
    `Standard LT AFT` = c(calc_stats_af("LT_AF", "AF_X", True_AF), calc_stats_af("LT_AF", "AF_Age", 0.85),
                          calc_stats_af("LT_AF", "AF_Sex", 1.10), calc_stats_af("LT_AF", "AF_READ", 0.90),
                          calc_stats_af("LT_AF", "AF_COADREAD", 0.95), calc_stats_metrics("LT_Metrics", "Med", mean_true_med),
                          calc_stats_metrics("LT_Metrics", "S5", mean_true_s5), "N/A"),
    `Proposed..Total.Effect.` = c(paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"),
                                  paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"),
                                  paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"),
                                  paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"),
                                  paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"),
                                  paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"),
                                  paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5", mean_true_s5), "</b>"),
                                  paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>"))
  )
}, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

output$sim_multi_fig1 <- renderPlot({
  req(res_multi())
  results <- res_multi()
  t_seq <- seq(0, 5, length.out = 100)

  safe_curve_mean <- function(getter) {
    mat <- do.call(rbind, lapply(results, function(x) { v <- getter(x); if (is.null(v) || length(v) != 100) return(rep(NA_real_, 100)); return(v) }))
    colMeans(mat, na.rm = TRUE)
  }

  plot_df1 <- data.frame(
    Time = rep(t_seq, 4),
    Survival = c(safe_curve_mean(function(x) x$curve_true), safe_curve_mean(function(x) x$curve_naive),
                 safe_curve_mean(function(x) x$curve_lt), safe_curve_mean(function(x) x$curve_prop)),
    Model = factor(rep(c("1. True Marginal", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed Method"), each = 100))
  )

  ggplot(plot_df1, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("black", "gray50", "gray50", "black")) +
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "twodash")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    theme_minimal(base_size = 15) +
    labs(title = "Figure 1: Reconstructed Marginal Survival Curves", subtitle = "Averaged across 400 iterations",
         x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
})

output$sim_multi_fig2 <- renderPlot({
  req(res_multi(), stored_inputs()$af)
  results <- res_multi()
  True_AF <- stored_inputs()$af
  pull_vals <- function(getter) { v <- sapply(results, getter); v[!is.na(v)] }

  af_naive <- pull_vals(function(x) x$Naive_AF["AF_X_est"])
  af_lt    <- pull_vals(function(x) x$LT_AF["AF_X_est"])
  af_prop  <- pull_vals(function(x) x$Prop_AF["AF_X_est"])

  plot_df2 <- data.frame(
    AF = c(af_naive, af_lt, af_prop),
    Model = factor(rep(c("Naive AFT", "Standard LT AFT", "Proposed Method"), c(length(af_naive), length(af_lt), length(af_prop))), levels = c("Naive AFT", "Standard LT AFT", "Proposed Method"))
  )

  ggplot(plot_df2, aes(x = Model, y = AF, fill = Model)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 1) +
    geom_hline(yintercept = True_AF, color = "black", linetype = "dashed", linewidth = 1.2) +
    scale_fill_grey(start = 0.9, end = 0.4) +
    coord_cartesian(ylim = c(0, max(True_AF * 2.5, 3))) +
    theme_minimal(base_size = 15) +
    labs(title = "Figure 2: Distribution of Estimated AF", subtitle = "Dashed line represents True Target Gene AF",
         x = "", y = "Estimated Acceleration Factor (AF)") +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
})

output$sim_multi_fig3 <- renderPlot({
  req(res_multi())
  results <- res_multi()

  t1_samp <- results[[1]]$sample_t1 / 365.25
  t2_samp <- results[[1]]$sample_t2 / 365.25
  plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)

  ggplot(plot_df3, aes(x = T1, y = T2)) +
    geom_point(color = "black", alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE) +
    theme_minimal(base_size = 15) +
    labs(title = "Figure 3: Dependent Left Truncation", subtitle = "Correlation between Time to CGP (T1) and Residual Survival (T2)",
         x = "Time to CGP Test (T1, Years)", y = "Residual Survival Time (T2, Years)") +
    theme(plot.title = element_text(face = "bold"))
})
