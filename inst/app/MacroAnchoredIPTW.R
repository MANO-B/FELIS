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

  sim_naive <- if(!is.null(fit_naive)) gen_sim_times(fit_naive, Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))
  sim_lt <- if(!is.null(fit_lt)) gen_sim_times(fit_lt, Data_cgp, "llogis") else rep(NA_real_, nrow(Data_cgp))

  # 省メモリ化: 100点の曲線データのみを抽出して返す
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

  out <- list(
    ESS = ESS,
    True_Metrics = calc_marginal_metrics(T_true),
    Naive_AF = extract_af_full(fit_naive), Naive_Metrics = calc_marginal_metrics(sim_naive),
    LT_AF = extract_af_full(fit_lt), LT_Metrics = calc_marginal_metrics(sim_lt),
    Prop_AF = extract_af_full(fit_prop_final), Prop_Metrics = calc_marginal_metrics(pseudo_cgp$sim_OS),

    # 図描画用のデータ抽出
    curve_true = calc_marginal_curve(T_true),
    curve_naive = calc_marginal_curve(sim_naive),
    curve_lt = calc_marginal_curve(sim_lt),
    curve_prop = calc_marginal_curve(pseudo_cgp$sim_OS),

    # 散布図用に最初の500例だけT1とT2を保存
    sample_t1 = Data_cgp$T1[1:min(500, nrow(Data_cgp))],
    sample_t2 = Data_cgp$T2_true[1:min(500, nrow(Data_cgp))]
  )
  out
}

safe_fmt <- function(x, digits = 2) {
  sapply(x, function(v) if (is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}

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
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Effective Sample Size (ESS)"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(res$True_Metrics["Med"]), safe_fmt(res$True_Metrics["S5"], 1), "N/A"),
      `Naive AFT` = c(format_af_single(res$Naive_AF, "AF_X"), format_af_single(res$Naive_AF, "AF_Age"), format_af_single(res$Naive_AF, "AF_Sex"), format_af_single(res$Naive_AF, "AF_READ"), format_af_single(res$Naive_AF, "AF_COADREAD"), safe_fmt(res$Naive_Metrics["Med"]), safe_fmt(res$Naive_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)),
      `Standard LT AFT` = c(format_af_single(res$LT_AF, "AF_X"), format_af_single(res$LT_AF, "AF_Age"), format_af_single(res$LT_AF, "AF_Sex"), format_af_single(res$LT_AF, "AF_READ"), format_af_single(res$LT_AF, "AF_COADREAD"), safe_fmt(res$LT_Metrics["Med"]), safe_fmt(res$LT_Metrics["S5"], 1), safe_fmt(input$sim_n, 0)),
      `Proposed (Ver 2.3.2)` = c(paste0("<b>", format_af_single(res$Prop_AF, "AF_X"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_Age"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_Sex"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_READ"), "</b>"), paste0("<b>", format_af_single(res$Prop_AF, "AF_COADREAD"), "</b>"), paste0("<b>", safe_fmt(res$Prop_Metrics["Med"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_Metrics["S5"], 1), "</b>"), paste0("<b>", safe_fmt(res$ESS, 0), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)
    plot_df <- data.frame(
      Time = rep(t_seq, 4),
      Survival = c(res$curve_true, res$curve_naive, res$curve_lt, res$curve_prop),
      Model = factor(rep(c("1. True Marginal (Macro)", "2. Naive AFT (Immortal Bias)", "3. Standard LT (Dep. Trunc Bias)", "4. Proposed Method (Ver 2.3.2)"), each = 100))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Ver 2.3.2 Validation", subtitle = "Single Run Results", x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
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

  # AFの集計関数: 点推定(Mean), MSE, Coverage Rate(CR)
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

  # OS指標の集計関数 (True_Metricsからの抽出を修正)
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

  True_AF <- input$sim_true_af

  # 【修正】True Valueの指標を安全に抽出して平均化
  mean_true_med <- mean(pull_vals(function(x) x$True_Metrics["Med"]), na.rm = TRUE)
  mean_true_s5  <- mean(pull_vals(function(x) x$True_Metrics["S5"]), na.rm = TRUE)

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]", "Average ESS"),

      # True Value には MSE を計算させない (ただの平均値を表示)
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", sprintf("%.2f", mean_true_med), sprintf("%.2f", mean_true_s5), "N/A"),

      `Naive AFT` = c(calc_stats_af("Naive_AF", "AF_X", True_AF), calc_stats_af("Naive_AF", "AF_Age", 0.85), calc_stats_af("Naive_AF", "AF_Sex", 1.10), calc_stats_af("Naive_AF", "AF_READ", 0.90), calc_stats_af("Naive_AF", "AF_COADREAD", 0.95), calc_stats_metrics("Naive_Metrics", "Med", mean_true_med), calc_stats_metrics("Naive_Metrics", "S5", mean_true_s5), "N/A"),

      `Standard LT AFT` = c(calc_stats_af("LT_AF", "AF_X", True_AF), calc_stats_af("LT_AF", "AF_Age", 0.85), calc_stats_af("LT_AF", "AF_Sex", 1.10), calc_stats_af("LT_AF", "AF_READ", 0.90), calc_stats_af("LT_AF", "AF_COADREAD", 0.95), calc_stats_metrics("LT_Metrics", "Med", mean_true_med), calc_stats_metrics("LT_Metrics", "S5", mean_true_s5), "N/A"),

      `Proposed (Ver 2.3.2)` = c(paste0("<b>", calc_stats_af("Prop_AF", "AF_X", True_AF), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats_af("Prop_AF", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "Med", mean_true_med), "</b>"), paste0("<b>", calc_stats_metrics("Prop_Metrics", "S5", mean_true_s5), "</b>"), paste0("<b>", sprintf("%.2f", mean(pull_vals(function(x) x$ESS))), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  # --- Fig 1: 400回の平均生存曲線 (安全な抽出処理) ---
  output$sim_multi_fig1 <- renderPlot({
    t_seq <- seq(0, 5, length.out = 100)

    safe_curve_mean <- function(getter) {
      mat <- do.call(rbind, lapply(results, function(x) {
        v <- getter(x)
        if (is.null(v) || length(v) != 100) return(rep(NA_real_, 100))
        return(v)
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
      Model = factor(rep(c("1. True Marginal", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed Method"), each = 100))
    )

    ggplot(plot_df1, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 1: Reconstructed Marginal Survival Curves", subtitle = "Averaged across 400 iterations", x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })

  # --- Fig 2: AF推定値の箱ひげ図 ---
  output$sim_multi_fig2 <- renderPlot({
    af_naive <- pull_vals(function(x) x$Naive_AF["AF_X_est"])
    af_lt    <- pull_vals(function(x) x$LT_AF["AF_X_est"])
    af_prop  <- pull_vals(function(x) x$Prop_AF["AF_X_est"])

    plot_df2 <- data.frame(
      AF = c(af_naive, af_lt, af_prop),
      Model = factor(rep(c("Naive AFT", "Standard LT AFT", "Proposed Method"),
                         c(length(af_naive), length(af_lt), length(af_prop))),
                     levels = c("Naive AFT", "Standard LT AFT", "Proposed Method"))
    )

    ggplot(plot_df2, aes(x = Model, y = AF, fill = Model)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 1) +
      geom_hline(yintercept = True_AF, color = "black", linetype = "dashed", linewidth = 1.2) +
      scale_fill_manual(values = c("#e74c3c", "#f39c12", "#27ae60")) +
      coord_cartesian(ylim = c(0, max(True_AF * 2.5, 3))) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 2: Distribution of Estimated AF", subtitle = "Dashed line represents True Target Gene AF", x = "", y = "Estimated Acceleration Factor (AF)") +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")
  })

  # --- Fig 3: 依存性切断の散布図 (T1 vs T2) ---
  output$sim_multi_fig3 <- renderPlot({
    t1_samp <- results[[1]]$sample_t1 / 365.25
    t2_samp <- results[[1]]$sample_t2 / 365.25
    plot_df3 <- data.frame(T1 = t1_samp, T2 = t2_samp)

    ggplot(plot_df3, aes(x = T1, y = T2)) +
      geom_point(color = "#34495e", alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", color = "red", linetype = "dashed", linewidth = 1) +
      theme_minimal(base_size = 15) +
      labs(title = "Figure 3: Dependent Left Truncation", subtitle = "Correlation between Time to CGP (T1) and Residual Survival (T2)", x = "Time to CGP Test (T1, Years)", y = "Residual Survival Time (T2, Years)") +
      theme(plot.title = element_text(face = "bold"))
  })
})
