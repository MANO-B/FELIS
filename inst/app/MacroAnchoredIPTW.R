# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.3.1)
# =========================================================================
library(dplyr)
library(flexsurv)
library(survival)
library(ggplot2)
library(splines)

# =========================================================================
# GLOBAL HELPERS
# =========================================================================

# Censoring parameter finder (randomness stabilized inside uniroot where needed)
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
    # NOTE: for indep/early we don't need u_fixed because we don't mix; still OK
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- tryCatch(uniroot(obj_func_scale, interval = c(0.1, 1e6))$root, error = function(e) median(T2) * 1.5)
    list(shape = if (pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE)
  }
}

# Parametric AFT sampler: Weibull or Log-logistic (flexsurvreg-style scale+shape)
gen_sim_times <- function(fit, newdata, dist_type = c("weibull", "llogis")) {
  dist_type <- match.arg(dist_type)
  if (is.null(fit) || nrow(newdata) == 0) return(rep(NA_real_, nrow(newdata)))

  res_m <- fit$res
  if (!all(c("scale", "shape") %in% rownames(res_m))) return(rep(NA_real_, nrow(newdata)))

  base_scale <- res_m["scale", "est"]
  shape <- res_m["shape", "est"]

  af_vec <- rep(1, nrow(newdata))

  # Standard covariates (match your factor coding)
  if ("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
  if ("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
  if ("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
  if ("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
  if ("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))

  # Extra vars (e.g., ns1..ns3)
  extra_vars <- setdiff(rownames(res_m), c("scale", "shape", "XMutated", "Age", "SexFemale", "HistologyREAD", "HistologyCOADREAD"))
  for (v in extra_vars) {
    if (v %in% names(newdata)) {
      af_vec <- af_vec * exp(res_m[v, "est"] * as.numeric(newdata[[v]]))
    }
  }

  scale_vec <- base_scale * af_vec

  u <- runif(nrow(newdata))
  u <- pmax(pmin(u, 0.9999), 0.0001)

  if (dist_type == "weibull") {
    # t = scale * (-log(u))^(1/shape)
    scale_vec * (-log(u))^(1 / shape)
  } else {
    # log-logistic inverse CDF: t = scale * ((1-u)/u)^(1/shape)
    scale_vec * ((1 - u) / u)^(1 / shape)
  }
}

# =========================================================================
# STEP 1: Calibration Weighting (KM at 1-5y matches macro, within Age_class)
# =========================================================================
calculate_calibrated_iptw <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points_days <- (1:5) * 365.25

  for (ag in unique(data$Age_class)) {
    if (!(ag %in% names(ref_surv_list))) next

    S_macro_target <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
    ag_data_idx <- which(data$Age_class == ag)
    ag_data <- data[ag_data_idx, , drop = FALSE]

    # weight model is function of log(T1) (flexible enough but stable)
    log_t1 <- log(pmax(ag_data$T1 / 365.25, 1e-6))
    log_t1_sq <- log_t1^2

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t1 + theta[2] * log_t1_sq)
      w <- pmax(w, 1e-4)
      w <- w / mean(w)

      km <- tryCatch(survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w), error = function(e) NULL)
      if (is.null(km)) return(1e6)

      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y

      # SSE + mild ridge
      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(optim(c(0, 0), obj_func, method = "BFGS"), error = function(e) list(par = c(0, 0)))
    theta_hat <- opt$par

    w_opt <- exp(theta_hat[1] * log_t1 + theta_hat[2] * log_t1_sq)
    w_opt <- pmax(w_opt, 1e-4)

    # winsorize extreme weights, then normalize mean=1
    lower_bound <- quantile(w_opt, 0.01, na.rm = TRUE)
    upper_bound <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lower_bound, pmin(w_opt, upper_bound))
    w_opt <- w_opt / mean(w_opt)

    data$iptw[ag_data_idx] <- w_opt
  }

  data
}

# =========================================================================
# Core Simulation Engine
# =========================================================================
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  # --- Macro population ---
  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  # AF for progression speed (your assumption: applies to T1 and T2 equally)
  AF_bg <- 1.0 *
    (0.85 ^ ((Age - 60) / 10)) *
    ifelse(Sex == "Female", 1.10, 1.0) *
    ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))

  AF_X_vec <- ifelse(X == 1, True_AF_X, 1.0)
  AF_total <- AF_bg * AF_X_vec

  # ---- FIX-1: Correct log-logistic baseline OS generation (scale=2y, shape=1.5) ----
  u <- runif(N_macro)
  u <- pmax(pmin(u, 0.999), 0.001)
  T_true_base <- (2.0 * 365.25) * ((1 - u) / u)^(1 / 1.5)

  # Baseline T1 pattern in AF=1 scale
  if (is.null(t1_pat)) t1_pat <- "indep"
  if (t1_pat == "early") {
    T1_base <- runif(N_macro, 30, pmax(31, T_true_base * 0.3))
  } else if (t1_pat == "dep_1yr" || t1_pat == "real") {
    T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat == "dep_2yr" || t1_pat == "rev") {
    T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else {
    T1_base <- runif(N_macro, 30, pmax(31, T_true_base))
  }

  T2_base <- pmax(0.1, T_true_base - T1_base)

  # Apply AF equally to T1 and T2 (your progression-speed assumption)
  T1 <- T1_base * AF_total
  T2_true <- T2_base * AF_total
  T_true <- T1 + T2_true

  # --- CGP inclusion: must survive to T1 (always true here since T_true=T1+T2_true) ---
  cgp_indices <- which(T_true > T1)
  if (length(cgp_indices) < N_target) return(NULL)

  # selection bias: younger oversampled
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

  # --- Censoring on T2 ---
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

  # --- Macro survival (age-stratified) targets (1-5y) ---
  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100)
  }

  # --- Calibration weights ---
  Data_cgp <- calculate_calibrated_iptw(Data_cgp, ref_surv_list)

  # ESS after calibration (mean-normalized weights)
  ESS <- (sum(Data_cgp$iptw, na.rm = TRUE)^2) / sum(Data_cgp$iptw^2, na.rm = TRUE)

  # --- covariates ---
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # spline basis from observed T1 (scale in years)
  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  ns_obj <- ns(Data_cgp$logT1_scale, df = 3)
  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for (j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  # --- Comparator models (observational world) ---
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))

  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt <- tryCatch(flexsurvreg(form_lt, data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  # --- Proposed: T1 model (weibull, always-event) + T2 model (llogis + spline) ---
  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurvreg(form_t1, data = Data_cgp, weights = iptw, dist = "weibull"), error = function(e) NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurvreg(form_t2, data = Data_cgp, weights = iptw, dist = "llogis"), error = function(e) NULL)

  # --- G-computation: weighted resampling ---
  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = 10000, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_cgp <- Data_cgp[idx_pseudo, , drop = FALSE]

  pseudo_cgp$sim_T1 <- gen_sim_times(fit_t1, pseudo_cgp, dist_type = "weibull")

  # regenerate spline terms from simulated T1
  pseudo_cgp$logT1_sim_scale <- log(pmax(pseudo_cgp$sim_T1 / 365.25, 1e-6))
  ns_sim <- tryCatch(predict(ns_obj, pseudo_cgp$logT1_sim_scale),
                     error = function(e) matrix(0, nrow(pseudo_cgp), attr(ns_obj, "df")))
  ns_sim <- as.matrix(ns_sim)
  colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
  for (j in seq_len(ncol(ns_sim))) pseudo_cgp[[colnames(ns_sim)[j]]] <- ns_sim[, j]

  pseudo_cgp$sim_T2 <- gen_sim_times(fit_t2, pseudo_cgp, dist_type = "llogis")
  pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
  pseudo_cgp$sim_Event <- 1

  # --- Final synced model on sim_OS (for AF output consistency on OS scale) ---
  form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
  fit_prop_final <- tryCatch(flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis"), error = function(e) NULL)

  # --- helpers for metrics ---
  calc_marginal_metrics <- function(sim_times) {
    if (is.null(sim_times) || all(is.na(sim_times))) return(c(Med = NA_real_, S5 = NA_real_))
    c(
      Med = median(sim_times, na.rm = TRUE) / 365.25,
      S5  = mean(sim_times > 5 * 365.25, na.rm = TRUE) * 100
    )
  }

  extract_af <- function(fit) {
    if (is.null(fit)) return(c(AF_X = NA, AF_Age = NA, AF_Sex = NA, AF_READ = NA, AF_COADREAD = NA))
    res_m <- fit$res
    get_val <- function(rname) if (rname %in% rownames(res_m)) res_m[rname, "est"] else NA_real_
    c(
      AF_X = exp(get_val("XMutated")),
      AF_Age = exp(get_val("Age") * 10),
      AF_Sex = exp(get_val("SexFemale")),
      AF_READ = exp(get_val("HistologyREAD")),
      AF_COADREAD = exp(get_val("HistologyCOADREAD"))
    )
  }

  # For Naive/LT "predicted curves", keep your existing interpretation:
  # - show what survival::analysis would imply if you (wrongly) treat those models as generative
  gen_naive_lt_sim <- function(fit) {
    if (is.null(fit)) return(rep(NA_real_, nrow(Data_cgp)))
    gen_sim_times(fit, Data_cgp, dist_type = "llogis")
  }
  sim_naive <- gen_naive_lt_sim(fit_naive)
  sim_lt <- gen_naive_lt_sim(fit_lt)

  out <- list(
    ESS = ESS,

    True_Metrics = calc_marginal_metrics(T_true),
    Naive_AF = extract_af(fit_naive),
    Naive_Metrics = calc_marginal_metrics(sim_naive),

    LT_AF = extract_af(fit_lt),
    LT_Metrics = calc_marginal_metrics(sim_lt),

    # FIX-5: Proposed AF from final synced OS model
    Prop_AF = extract_af(fit_prop_final),
    Prop_Metrics = calc_marginal_metrics(pseudo_cgp$sim_OS),

    # Plot sources
    t_true_plot = T_true,
    t_naive_plot = sim_naive,
    t_lt_plot = sim_lt,
    t_prop_plot = pseudo_cgp$sim_OS
  )

  if (return_data) out else out
}

safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) {
    if (is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v)
  })
}

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic (keep UI variables/output names)
# -------------------------------------------------------------------------

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation (Tamura & Ikegami Ver 2.3.1)...", type = "message", duration = 3)

  # FIX-2: Sys.time() -> integer seed (still "random each click")
  set.seed(as.integer(Sys.time()))

  res <- run_sim_iteration(
    input$sim_n,
    input$sim_true_af,
    input$sim_mut_freq / 100,
    input$sim_t1_pattern,
    input$sim_cens_pattern,
    input$sim_cens_rate / 100,
    return_data = TRUE
  )
  shiny::validate(shiny::need(!is.null(res), "Simulation failed. Not enough events generated."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c(
        "Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF",
        "Histology (READ) AF", "Histology (COADREAD) AF",
        "Marginal Med OS [Years]", "Marginal 5-Year OS [%]",
        "Effective Sample Size (ESS)"
      ),
      `True` = c(
        safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
        safe_fmt(res$True_Metrics["Med"]),
        safe_fmt(res$True_Metrics["S5"], 1),
        "N/A"   # FIX-3: True ESS is not defined here
      ),
      `Naive AFT` = c(
        safe_fmt(res$Naive_AF["AF_X"]),
        safe_fmt(res$Naive_AF["AF_Age"]),
        safe_fmt(res$Naive_AF["AF_Sex"]),
        safe_fmt(res$Naive_AF["AF_READ"]),
        safe_fmt(res$Naive_AF["AF_COADREAD"]),
        safe_fmt(res$Naive_Metrics["Med"]),
        safe_fmt(res$Naive_Metrics["S5"], 1),
        safe_fmt(input$sim_n, 0)  # "no weighting" reference
      ),
      `Standard LT` = c(
        safe_fmt(res$LT_AF["AF_X"]),
        safe_fmt(res$LT_AF["AF_Age"]),
        safe_fmt(res$LT_AF["AF_Sex"]),
        safe_fmt(res$LT_AF["AF_READ"]),
        safe_fmt(res$LT_AF["AF_COADREAD"]),
        safe_fmt(res$LT_Metrics["Med"]),
        safe_fmt(res$LT_Metrics["S5"], 1),
        safe_fmt(input$sim_n, 0)  # "no weighting" reference
      ),
      `Proposed (Ver 2.3.1)` = c(
        paste0("<b>", safe_fmt(res$Prop_AF["AF_X"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_AF["AF_Age"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_AF["AF_Sex"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_AF["AF_READ"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_AF["AF_COADREAD"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"),
        paste0("<b>", safe_fmt(res$ESS, 0), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    calc_marginal <- function(t_data) {
      if (is.null(t_data) || all(is.na(t_data))) return(rep(NA_real_, length(t_seq)))
      km <- survfit(Surv(t_data, rep(1, length(t_data))) ~ 1)
      approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
    }

    plot_df <- data.frame(
      Time = rep(t_seq / 365.25, 4),
      Survival = c(
        calc_marginal(res$t_true_plot),
        calc_marginal(res$t_naive_plot),
        calc_marginal(res$t_lt_plot),
        calc_marginal(res$t_prop_plot)
      ),
      Model = factor(rep(c(
        "1. True Marginal (Registry Macro)",
        "2. Naive AFT (Immortal Time Bias)",
        "3. Standard LT AFT (Dep. Truncation Bias)",
        "4. Proposed Method (Ver 2.3.1: Calibrated G-comp)"
      ), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(
        title = "Tamura & Ikegami Model Ver 2.3.1 Validation",
        subtitle = "DGP: AF applies to both T1 and T2 (progression-speed model).",
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
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)

  n_sims <- 400
  results <- list()

  withProgress(message = "Running 400 Simulations...", value = 0, {
    for (i in 1:n_sims) {
      set.seed(i)  # reproducible across runs, varies across iterations

      res <- run_sim_iteration(
        input$sim_n,
        input$sim_true_af,
        input$sim_mut_freq / 100,
        input$sim_t1_pattern,
        input$sim_cens_pattern,
        input$sim_cens_rate / 100,
        return_data = FALSE
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res

      incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })

  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  True_AF <- input$sim_true_af

  # FIX-4: proper aggregation for AF, Med, S5, ESS
  pull_vals <- function(getter) {
    v <- sapply(results, getter)
    v[!is.na(v)]
  }

  stat_mse <- function(vals, true_val) {
    if (length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    mse_val <- mean((vals - true_val)^2)
    sprintf("%.2f (MSE: %.3f)", mean_val, mse_val)
  }

  stat_mean <- function(vals, digits = 2) {
    if (length(vals) == 0) return("N/A")
    sprintf(paste0("%.", digits, "f"), mean(vals))
  }

  # True metrics estimated from macro each iteration (Monte Carlo truth)
  true_med_vals <- pull_vals(function(x) x$True_Metrics["Med"])
  true_s5_vals  <- pull_vals(function(x) x$True_Metrics["S5"])

  # Naive
  naive_af_x <- pull_vals(function(x) x$Naive_AF["AF_X"])
  naive_af_age <- pull_vals(function(x) x$Naive_AF["AF_Age"])
  naive_af_sex <- pull_vals(function(x) x$Naive_AF["AF_Sex"])
  naive_af_read <- pull_vals(function(x) x$Naive_AF["AF_READ"])
  naive_af_coadread <- pull_vals(function(x) x$Naive_AF["AF_COADREAD"])
  naive_med <- pull_vals(function(x) x$Naive_Metrics["Med"])
  naive_s5  <- pull_vals(function(x) x$Naive_Metrics["S5"])

  # LT
  lt_af_x <- pull_vals(function(x) x$LT_AF["AF_X"])
  lt_af_age <- pull_vals(function(x) x$LT_AF["AF_Age"])
  lt_af_sex <- pull_vals(function(x) x$LT_AF["AF_Sex"])
  lt_af_read <- pull_vals(function(x) x$LT_AF["AF_READ"])
  lt_af_coadread <- pull_vals(function(x) x$LT_AF["AF_COADREAD"])
  lt_med <- pull_vals(function(x) x$LT_Metrics["Med"])
  lt_s5  <- pull_vals(function(x) x$LT_Metrics["S5"])

  # Proposed
  prop_af_x <- pull_vals(function(x) x$Prop_AF["AF_X"])
  prop_af_age <- pull_vals(function(x) x$Prop_AF["AF_Age"])
  prop_af_sex <- pull_vals(function(x) x$Prop_AF["AF_Sex"])
  prop_af_read <- pull_vals(function(x) x$Prop_AF["AF_READ"])
  prop_af_coadread <- pull_vals(function(x) x$Prop_AF["AF_COADREAD"])
  prop_med <- pull_vals(function(x) x$Prop_Metrics["Med"])
  prop_s5  <- pull_vals(function(x) x$Prop_Metrics["S5"])
  ess_vals <- pull_vals(function(x) x$ESS)

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c(
        "Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF",
        "Histology (READ) AF", "Histology (COADREAD) AF",
        "Marginal Med OS [Years]", "Marginal 5-Year OS [%]",
        "Average ESS"
      ),
      `True Value` = c(
        safe_fmt(True_AF),
        "0.85", "1.10", "0.90", "0.95",
        stat_mean(true_med_vals, 2),
        stat_mean(true_s5_vals, 1),
        "N/A"
      ),
      `Naive AFT` = c(
        stat_mse(naive_af_x, True_AF),
        stat_mse(naive_af_age, 0.85),
        stat_mse(naive_af_sex, 1.10),
        stat_mse(naive_af_read, 0.90),
        stat_mse(naive_af_coadread, 0.95),
        stat_mean(naive_med, 2),
        stat_mean(naive_s5, 1),
        "N/A"
      ),
      `Standard LT AFT` = c(
        stat_mse(lt_af_x, True_AF),
        stat_mse(lt_af_age, 0.85),
        stat_mse(lt_af_sex, 1.10),
        stat_mse(lt_af_read, 0.90),
        stat_mse(lt_af_coadread, 0.95),
        stat_mean(lt_med, 2),
        stat_mean(lt_s5, 1),
        "N/A"
      ),
      `Proposed (Ver 2.3.1)` = c(
        paste0("<b>", stat_mse(prop_af_x, True_AF), "</b>"),
        paste0("<b>", stat_mse(prop_af_age, 0.85), "</b>"),
        paste0("<b>", stat_mse(prop_af_sex, 1.10), "</b>"),
        paste0("<b>", stat_mse(prop_af_read, 0.90), "</b>"),
        paste0("<b>", stat_mse(prop_af_coadread, 0.95), "</b>"),
        paste0("<b>", stat_mean(prop_med, 2), "</b>"),
        paste0("<b>", stat_mean(prop_s5, 1), "</b>"),
        paste0("<b>", stat_mean(ess_vals, 0), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
