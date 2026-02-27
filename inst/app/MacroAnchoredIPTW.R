# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.3.2)
# =========================================================================
library(dplyr)
library(flexsurv)
library(survival)
library(ggplot2)
library(splines)

find_censoring_param <- function(T2, target_rate, pattern) {
  if (pattern == "peak1y") {
    obj_func <- function(k) { scale_val <- 365.25 / ((k - 1) / k)^(1 / k); (sum(rweibull(length(T2), shape = k, scale = scale_val) < T2) / length(T2)) - target_rate }
    opt_k <- tryCatch(uniroot(obj_func, interval = c(1.05, 10))$root, error = function(e) 2.0)
    list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE)
  } else if (pattern == "ushape") {
    u_fixed <- runif(length(T2))
    obj_func <- function(scale_param) {
      C2 <- ifelse(u_fixed < 0.5, rweibull(length(T2), shape = 0.5, scale = scale_param), rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
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

gen_sim_times <- function(fit, newdata, dist_type = c("weibull", "llogis")) {
  dist_type <- match.arg(dist_type)
  if (is.null(fit) || nrow(newdata) == 0) return(rep(NA_real_, nrow(newdata)))
  res_m <- fit$res
  if (!all(c("scale", "shape") %in% rownames(res_m))) return(rep(NA_real_, nrow(newdata)))

  base_scale <- res_m["scale", "est"]
  shape <- res_m["shape", "est"]
  af_vec <- rep(1, nrow(newdata))

  if ("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
  if ("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
  if ("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
  if ("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
  if ("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))

  extra_vars <- setdiff(rownames(res_m), c("scale", "shape", "XMutated", "Age", "SexFemale", "HistologyREAD", "HistologyCOADREAD"))
  for (v in extra_vars) {
    if (v %in% names(newdata)) af_vec <- af_vec * exp(res_m[v, "est"] * as.numeric(newdata[[v]]))
  }

  scale_vec <- base_scale * af_vec
  u <- pmax(pmin(runif(nrow(newdata)), 0.9999), 0.0001)

  if (dist_type == "weibull") { scale_vec * (-log(u))^(1 / shape) }
  else { scale_vec * ((1 - u) / u)^(1 / shape) }
}

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
      km <- tryCatch(survfit(Surv(T_obs, Event) ~ 1, data = ag_data, weights = w), error = function(e) NULL)
      if (is.null(km)) return(1e6)
      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(optim(c(0, 0), obj_func, method = "BFGS"), error = function(e) list(par = c(0, 0)))
    theta_hat <- opt$par
    w_opt <- pmax(exp(theta_hat[1] * log_t1 + theta_hat[2] * log_t1_sq), 1e-4)

    lower_bound <- quantile(w_opt, 0.01, na.rm = TRUE)
    upper_bound <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lower_bound, pmin(w_opt, upper_bound))
    data$iptw[ag_data_idx] <- w_opt / mean(w_opt)
  }
  data
}

# =========================================================================
# Core Simulation Engine
# =========================================================================
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, True_Med, True_Shape, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", 1.10, 1.0) * ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))
  AF_total <- AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  # UIパラメータを使ってベースラインOSを生成
  u <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true_base <- (True_Med * 365.25) * ((1 - u) / u)^(1 / True_Shape)

  if (is.null(t1_pat)) t1_pat <- "indep"
  if (t1_pat == "early") { T1_base <- runif(N_macro, 30, pmax(31, T_true_base * 0.3))
  } else if (t1_pat == "dep_1yr" || t1_pat == "real") { T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat == "dep_2yr" || t1_pat == "rev") { T1_base <- pmax(30, T_true_base - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else { T1_base <- runif(N_macro, 30, pmax(31, T_true_base)) }

  T2_base <- pmax(0.1, T_true_base - T1_base)

  T1 <- T1_base * AF_total
  T2_true <- T2_base * AF_total
  T_true <- T1 + T2_true

  cgp_indices <- which(T_true > T1)
  if (length(cgp_indices) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.9, 0.1)
  selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target, Age = Age[selected_indices], Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")), Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")), T_true = T_true[selected_indices], T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

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

  ref_surv_list <- list()
  age_label_macro <- ifelse(Age < 60, "Young", "Old")
  for (ag in c("Young", "Old")) { ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[age_label_macro == ag] > y * 365.25) * 100) }

  Data_cgp <- calculate_calibrated_iptw(Data_cgp, ref_surv_list)
  ESS <- (sum(Data_cgp$iptw, na.rm = TRUE)^2) / sum(Data_cgp$iptw^2, na.rm = TRUE)

  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  Data_cgp$logT1_scale <- log(pmax(Data_cgp$T1 / 365.25, 1e-6))
  ns_obj <- ns(Data_cgp$logT1_scale, df = 3)
  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for (j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse = " + ")))
  form_lt <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse = " + ")))

  fit_naive <- tryCatch(flexsurvreg(form_naive, data = Data_cgp, dist = "llogis"), error = function(e) NULL)
  fit_lt <- tryCatch(flexsurvreg(form_lt, data = Data_cgp, dist = "llogis"), error = function(e) NULL)

  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurvreg(form_t1, data = Data_cgp, weights = iptw, dist = "weibull"), error = function(e) NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurvreg(form_t2, data = Data_cgp, weights = iptw, dist = "llogis"), error = function(e) NULL)

  idx_pseudo <- sample(seq_len(nrow(Data_cgp)), size = 10000, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_cgp <- Data_cgp[idx_pseudo, , drop = FALSE]

  pseudo_cgp$sim_T1 <- gen_sim_times(fit_t1, pseudo_cgp, dist_type = "weibull")
  pseudo_cgp$logT1_sim_scale <- log(pmax(pseudo_cgp$sim_T1 / 365.25, 1e-6))
  ns_sim <- tryCatch(predict(ns_obj, pseudo_cgp$logT1_sim_scale), error = function(e) matrix(0, nrow(pseudo_cgp), attr(ns_obj, "df")))
  colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
  for (j in seq_len(ncol(ns_sim))) pseudo_cgp[[colnames(ns_sim)[j]]] <- ns_sim[, j]

  pseudo_cgp$sim_T2 <- gen_sim_times(fit_t2, pseudo_cgp, dist_type = "llogis")
  pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
  pseudo_cgp$sim_Event <- 1

  form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse = " + ")))
  fit_prop_final <- tryCatch(flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis"), error = function(e) NULL)

  # --- Metrics extraction (95% CI を含むように拡張) ---
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
    get_est <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "est"] * mult) else NA_real_
    get_L <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "L95%"] * mult) else NA_real_
    get_U <- function(rname, mult=1) if (rname %in% rownames(res_m)) exp(res_m[rname, "U95%"] * mult) else NA_real_

    c(
      AF_X_est = get_est("XMutated"), AF_X_L = get_L("XMutated"), AF_X_U = get_U("XMutated"),
      AF_Age_est = get_est("Age", 10), AF_Age_L = get_L("Age", 10), AF_Age_U = get_U("Age", 10),
      AF_Sex_est = get_est("SexFemale"), AF_Sex_L = get_L("SexFemale"), AF_Sex_U = get_U("SexFemale"),
      AF_READ_est = get_est("HistologyREAD"), AF_READ_L = get_L("HistologyREAD"), AF_READ_U = get_U("HistologyREAD"),
      AF_COADREAD_est = get_est("HistologyCOADREAD"), AF_COADREAD_L = get_L("HistologyCOADREAD"), AF_COADREAD_U = get_U("HistologyCOADREAD")
    )
  }

  calc_marginal_metrics <- function(sim_times) {
    if (is.null(sim_times) || all(is.na(sim_times))) return(c(Med = NA_real_, S5 = NA_real_))
    c(Med = median(sim_times, na.rm = TRUE) / 365.25, S5 = mean(sim_times > 5 * 365.25, na.rm = TRUE) * 100)
  }

  sim_naive <- if(!is.null(fit_naive)) gen_sim_times(fit_naive, Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt <- if(!is.null(fit_lt)) gen_sim_times(fit_lt, Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))

  out <- list(
    ESS = ESS,
    True_Metrics = calc_marginal_metrics(T_true),
    Naive_AF = extract_af_full(fit_naive), Naive_Metrics = calc_marginal_metrics(sim_naive),
    LT_AF = extract_af_full(fit_lt), LT_Metrics = calc_marginal_metrics(sim_lt),
    Prop_AF = extract_af_full(fit_prop_final), Prop_Metrics = calc_marginal_metrics(pseudo_cgp$sim_OS),
    t_true_plot = T_true, t_naive_plot = sim_naive, t_lt_plot = sim_lt, t_prop_plot = pseudo_cgp$sim_OS
  )
  if (return_data) out else out
}

safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) if (is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}

# 単回実行用のフォーマット関数 (点推定と95%CI)
format_af_single <- function(af_vec, var_name) {
  est <- af_vec[paste0(var_name, "_est")]
  L <- af_vec[paste0(var_name, "_L")]
  U <- af_vec[paste0(var_name, "_U")]
  if (is.na(est)) return("N/A")
  sprintf("%.2f (%.2f-%.2f)", est, L, U)
}

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic
# -------------------------------------------------------------------------

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq, input$sim_true_med, input$sim_true_shape)
  showNotification("Running Single Simulation...", type = "message", duration = 3)

  set.seed(as.integer(Sys.time()))

  res <- run_sim_iteration(
    input$sim_n, input$sim_true_af, input$sim_mut_freq / 100,
    input$sim_true_med, input$sim_true_shape,
    input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE
  )
  shiny::validate(shiny::need(!is.null(res), "Simulation failed. Not enough events generated."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c(
        "Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF",
        "Histology (READ) AF", "Histology (COADREAD) AF",
        "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"
      ),
      `True` = c(
        safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95",
        safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"
      ),
      `Naive AFT` = c(
        format_af_single(res$Naive_AF, "AF_X"), format_af_single(res$Naive_AF, "AF_Age"), format_af_single(res$Naive_AF, "AF_Sex"),
        format_af_single(res$Naive_AF, "AF_READ"), format_af_single(res$Naive_AF, "AF_COADREAD"),
        safe_fmt(res$Naive_Metrics["Med"]), safe_fmt(res$Naive_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)
      ),
      `Standard LT AFT` = c(
        format_af_single(res$LT_AF, "AF_X"), format_af_single(res$LT_AF, "AF_Age"), format_af_single(res$LT_AF, "AF_Sex"),
        format_af_single(res$LT_AF, "AF_READ"), format_af_single(res$LT_AF, "AF_COADREAD"),
        safe_fmt(res$LT_Metrics["Med"]), safe_fmt(res$LT_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)
      ),
      `Proposed (Ver 2.3.2)` = c(
        paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_Age"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_Sex"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_READ"), "</b>"),
        paste0("<b>", format_af_single(res$Prop_AF, "AF_COADREAD"), "</b>"), paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"), paste0("<b>", safe_fmt(res$ESS, 0), "</b>")
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
      Survival = c(calc_marginal(res$t_true_plot), calc_marginal(res$t_naive_plot), calc_marginal(res$t_lt_plot), calc_marginal(res$t_prop_plot)),
      Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT (Immortal Time Bias)", "3. Standard LT AFT (Dep. Truncation Bias)", "4. Proposed Method (Ver 2.3.2: Calibrated G-comp)"), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Ver 2.3.2 Validation", subtitle = "DGP: AF applies to both T1 and T2. Single run includes 95% CIs.", x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
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
      res <- run_sim_iteration(
        input$sim_n, input$sim_true_af, input$sim_mut_freq / 100,
        input$sim_true_med, input$sim_true_shape,
        input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = FALSE
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1 / n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })

  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  pull_vals <- function(getter) { v <- sapply(results, getter); v[!is.na(v)] }

  # AFの集計関数: 点推定(Mean), MSE, そして Coverage Rate(CR) を計算
  calc_stats_af <- function(model_name, var_name, true_val) {
    est_vals <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_est")])
    l_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_L")])
    u_vals   <- pull_vals(function(x) x[[model_name]][paste0(var_name, "_U")])

    if (length(est_vals) == 0) return("N/A")
    mean_est <- mean(est_vals)
    mse_val <- mean((est_vals - true_val)^2)
    # L <= true <= U の割合を計算
    cr_val <- mean(l_vals <= true_val & true_val <= u_vals, na.rm=TRUE) * 100
    sprintf("%.2f (MSE: %.3f, CR: %.1f%%)", mean_est, mse_val, cr_val)
  }

  calc_stats_metrics <- function(model_name, metric_name, true_val = NULL) {
    if(!is.null(model_name)){ vals <- sapply(results, function(x) x[[model_name]][metric_name])
    } else { vals <- sapply(results, function(x) x[[metric_name]]) }
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    if(!is.null(true_val)){
      mse_val <- mean((vals - true_val)^2)
      return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
    }
    return(sprintf("%.2f", mean_val))
  }

  True_AF <- input$sim_true_af
  true_med_vals <- pull_vals(function(x) x$True_Metrics["Med"])
  true_s5_vals  <- pull_vals(function(x) x$True_Metrics["S5"])

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", calc_stats_metrics(NULL, "True_Metrics", 1), calc_stats_metrics(NULL, "True_Metrics", 2), "N/A"),

      `Naive AFT` = c(calc_stats_af("Naive_AF", "AF_X", True_AF), calc_stats_af("Naive_AF", "AF_Age", 0.85), calc_stats_af("Naive_AF", "AF_Sex", 1.10), calc_stats_af("Naive_AF", "AF_READ", 0.90), calc_stats_af("Naive_AF", "AF_COADREAD", 0.95), calc_stats_metrics("Naive_Metrics", "Med", mean(true_med_vals)), calc_stats_metrics("Naive_Metrics", "S5", mean(true_s5_vals)), "N/A"),

      `Standard LT AFT` = c(calc_stats_af("LT_AF", "AF_X", True_AF), calc_stats_af("LT_AF", "AF_Age", 0.85), calc_stats_af("LT_AF", "AF_Sex", 1.10), calc_stats_af("LT_AF", "AF_READ", 0.90), calc_stats_af("LT_AF", "AF_COADREAD", 0.95), calc_stats_metrics("LT_Metrics", "Med", mean(true_med_vals)), calc_stats_metrics("LT_Metrics", "S5", mean(true_s5_vals)), "N/A"),

      `Proposed (Ver 2.3.2)` = c(paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean(true_med_vals)), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5", mean(true_s5_vals)), "</b>"), paste0("<b>", calc_stats_metrics(NULL, "ESS"), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
