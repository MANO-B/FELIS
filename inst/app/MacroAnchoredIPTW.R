# =========================================================================
# Server Logic for Simulation Study (OS-Matching IPTW: True Tamura & Ikegami Method)
# =========================================================================
library(dplyr)
library(flexsurv)
library(survival)
library(ggplot2)

find_censoring_param <- function(T2, target_rate, pattern) {
  if (pattern == "peak1y") {
    obj_func <- function(k) {
      scale_val <- 365.25 / ((k - 1) / k)^(1 / k)
      C2 <- rweibull(length(T2), shape = k, scale = scale_val)
      return((sum(C2 < T2) / length(T2)) - target_rate)
    }
    opt_k <- tryCatch({ uniroot(obj_func, interval = c(1.05, 10))$root }, error = function(e) { 2.0 })
    opt_scale <- 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k)
    return(list(shape = opt_k, scale = opt_scale))
  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- rexp(length(T2), rate = 1 / scale_param)
      return((sum(C2 < T2) / length(T2)) - target_rate)
    }
    opt_scale <- tryCatch({ uniroot(obj_func_scale, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
    return(list(shape = 1.0, scale = opt_scale))
  }
}

run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  # =========================================================================
  # 1. Macro Population (Registry Anchor: Median OS ~2 years)
  # =========================================================================
  # ※全例がCGPを受けるわけではありません。10万人のレジストリから、後でN_target人だけが抽出されます。
  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  True_AF_Age_10yr <- 0.85
  True_AF_Sex_Female <- 1.10
  True_AF_Hist_READ <- 0.90
  True_AF_Hist_COADREAD <- 0.95

  AF_bg <- 1.0 * (True_AF_Age_10yr ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", True_AF_Sex_Female, 1.0) *
    ifelse(Histology == "READ", True_AF_Hist_READ, ifelse(Histology == "COADREAD", True_AF_Hist_COADREAD, 1.0))

  shape_base <- 1.5
  scale_base <- 2.0 * 365.25

  u <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true <- scale_base * ((1 - u) / u)^(1 / shape_base) * AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  # =========================================================================
  # 2. Timing of CGP (T1)
  # =========================================================================
  if (t1_pat == "indep") {
    T1 <- runif(N_macro, 30, pmax(31, T_true))
  } else if (t1_pat == "dep_1yr") {
    # 最悪の依存性切断: 死亡の約1年前にCGP
    T2_exp <- rlnorm(N_macro, meanlog = log(365.25 * 1.0), sdlog = 0.3)
    T1 <- pmax(30, T_true - T2_exp)
  } else if (t1_pat == "dep_2yr") {
    T2_exp <- rlnorm(N_macro, meanlog = log(365.25 * 2.0), sdlog = 0.3)
    T1 <- pmax(30, T_true - T2_exp)
  }

  cgp_indices <- which(T_true > T1)
  if(length(cgp_indices) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.8, 0.4)
  selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target,
    Age = Age[selected_indices],
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
    Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[selected_indices],
    T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # =========================================================================
  # 3. Censoring (C2)
  # =========================================================================
  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (cens_pat == "indep") { C2 <- rexp(N_target, rate = 1 / params$scale)
  } else { C2 <- rweibull(N_target, shape = params$shape, scale = params$scale) }

  Data_cgp$C2 <- C2
  Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)
  Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs # OS (全生存期間)
  if(sum(Data_cgp$Event) < 20) return(NULL)

  # =========================================================================
  # 4. 真の Tamura & Ikegami メソッド (OS-Matching IPTW)
  # 院内がん登録(Macro)のOSに完全一致させる重み付け
  # =========================================================================
  # OSのビン(半年刻み)を設定
  breaks <- c(seq(0, 10 * 365.25, by = 0.5 * 365.25), Inf)

  # Macro集団が各ビンで死亡する確率 (院内がん登録の真のOS分布)
  macro_counts <- hist(T_true, breaks = breaks, plot = FALSE)$counts
  p_macro <- macro_counts / sum(macro_counts)

  # CGP集団が各ビンで死亡する観測確率 (カプランマイヤーで打ち切りを考慮)
  km_cgp <- survfit(Surv(T_obs, Event) ~ 1, data = Data_cgp)
  surv_probs_cgp <- approx(c(0, km_cgp$time), c(1, km_cgp$surv), xout = breaks, method = "constant", f = 0, rule = 2)$y
  p_cgp <- -diff(surv_probs_cgp)
  p_cgp <- pmax(p_cgp, 1e-5) # ゼロ割り防止

  # 全生存期間(OS)の分布が一致するように重みを割り当て
  Data_cgp$os_bin <- cut(Data_cgp$T_obs, breaks = breaks, labels = FALSE)
  Data_cgp$raw_weight <- p_macro[Data_cgp$os_bin] / p_cgp[Data_cgp$os_bin]

  # 爆発を防ぐキャッピング処理 (必須)
  lower_cap <- quantile(Data_cgp$raw_weight, 0.05, na.rm = TRUE)
  upper_cap <- quantile(Data_cgp$raw_weight, 0.95, na.rm = TRUE)
  Data_cgp$raw_weight <- pmax(lower_cap, pmin(Data_cgp$raw_weight, upper_cap))
  Data_cgp$iptw <- Data_cgp$raw_weight / mean(Data_cgp$raw_weight)

  # =========================================================================
  # 5. Model Fitting
  # =========================================================================
  form <- as.formula("~ X + Age + Sex + Histology")

  # ① Naive (切断無視, 重みなし)
  fit_naive <- tryCatch({ flexsurvreg(update(Surv(T_obs, Event) ~ ., form), data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  # ② Standard LT (標準の左側切断モデル. 依存性切断の前では崩壊する)
  fit_lt <- tryCatch({ flexsurvreg(update(Surv(T1, T_obs, Event) ~ ., form), data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  # ③ Proposed Method (OS-Matching IPTW)
  # ※OSを一致させたことでT=0からの仮想コホートになっているため、LT(T1)は不要。通常のSurvを使用。
  fit_prop <- tryCatch({ flexsurvreg(update(Surv(T_obs, Event) ~ ., form), data = Data_cgp, weights = iptw, dist = "llogis") }, error = function(e) NULL)

  extract_metrics <- function(fit) {
    if (is.null(fit)) return(rep(NA, 11))
    res <- fit$res
    AF_X <- exp(res["XMutated", "est"]); AF_X_L <- exp(res["XMutated", "L95%"]); AF_X_U <- exp(res["XMutated", "U95%"])
    AF_Age <- exp(res["Age", "est"] * 10); AF_Sex <- exp(res["SexFemale", "est"])
    AF_READ <- exp(res["HistologyREAD", "est"]); AF_COADREAD <- exp(res["HistologyCOADREAD", "est"])

    scale_base <- res["scale", "est"]; shape <- res["shape", "est"]
    af_age60 <- exp(res["Age", "est"] * 60)
    med_wt_days <- scale_base * af_age60
    med_wt <- med_wt_days / 365.25; med_mut <- med_wt * AF_X
    s5_wt <- 1 / (1 + (5 * 365.25 / med_wt_days)^shape) * 100
    s5_mut <- 1 / (1 + (5 * 365.25 / (med_wt_days * AF_X))^shape) * 100

    return(c(AF_X=AF_X, AF_X_L=AF_X_L, AF_X_U=AF_X_U, AF_Age=AF_Age, AF_Sex=AF_Sex, AF_READ=AF_READ, AF_COADREAD=AF_COADREAD, Med_WT=med_wt, Med_Mut=med_mut, S5_WT=s5_wt, S5_Mut=s5_mut))
  }

  out <- list(Naive = extract_metrics(fit_naive), LT = extract_metrics(fit_lt), Prop = extract_metrics(fit_prop),
              fit_naive = fit_naive, fit_lt = fit_lt, fit_prop = fit_prop, Data_cgp = Data_cgp, T_true_macro = T_true)
  if(return_data) return(out) else return(out[1:3])
}

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic
# -------------------------------------------------------------------------
observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation with OS-Matching IPTW...", type = "message", duration = 3)

  res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE)
  shiny::validate(shiny::need(!is.null(res), "Simulation failed. Please try again."))

  True_AF <- input$sim_true_af; True_Med_WT <- 2.0; True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Median OS (WT) [Years]", "5-Year OS (WT) [%]"),
      `True` = c(sprintf("%.2f", True_AF), "0.85", "1.10", "0.90", "0.95", sprintf("%.2f", True_Med_WT), sprintf("%.1f", True_S5_WT)),
      `Naive AFT` = c(sprintf("%.2f", res$Naive["AF_X"]), sprintf("%.2f", res$Naive["AF_Age"]), sprintf("%.2f", res$Naive["AF_Sex"]), sprintf("%.2f", res$Naive["AF_READ"]), sprintf("%.2f", res$Naive["AF_COADREAD"]), sprintf("%.2f", res$Naive["Med_WT"]), sprintf("%.1f", res$Naive["S5_WT"])),
      `Standard LT` = c(sprintf("%.2f", res$LT["AF_X"]), sprintf("%.2f", res$LT["AF_Age"]), sprintf("%.2f", res$LT["AF_Sex"]), sprintf("%.2f", res$LT["AF_READ"]), sprintf("%.2f", res$LT["AF_COADREAD"]), sprintf("%.2f", res$LT["Med_WT"]), sprintf("%.1f", res$LT["S5_WT"])),
      `Proposed Method` = c(sprintf("<b>%.2f</b>", res$Prop["AF_X"]), sprintf("<b>%.2f</b>", res$Prop["AF_Age"]), sprintf("<b>%.2f</b>", res$Prop["AF_Sex"]), sprintf("<b>%.2f</b>", res$Prop["AF_READ"]), sprintf("<b>%.2f</b>", res$Prop["AF_COADREAD"]), sprintf("<b>%.2f</b>", res$Prop["Med_WT"]), sprintf("<b>%.1f</b>", res$Prop["S5_WT"]))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    library(ggplot2)
    library(survival)
    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    km_true <- survfit(Surv(res$T_true_macro, rep(1, length(res$T_true_macro))) ~ 1)
    surv_true <- approx(km_true$time, km_true$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y

    calc_marginal <- function(fit) {
      if(is.null(fit)) return(rep(NA, length(t_seq)))
      samp_data <- res$Data_cgp[sample(1:nrow(res$Data_cgp), min(300, nrow(res$Data_cgp))), ]
      summ <- summary(fit, newdata = samp_data, t = t_seq, type = "survival", tidy = TRUE)
      agg <- aggregate(est ~ time, data = summ, FUN = mean)
      return(agg$est)
    }

    plot_df <- data.frame(
      Time = rep(t_seq / 365.25, 4),
      Survival = c(surv_true, calc_marginal(res$fit_naive), calc_marginal(res$fit_lt), calc_marginal(res$fit_prop)),
      Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT", "3. Standard LT AFT", "4. Proposed Method (OS-Anchored IPTW)"), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Marginal Survival Recovery (OS Matching Validation)",
           subtitle = "Proposed Method forces baseline alignment to registry, completely neutralizing dependent truncation.",
           x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
  })
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  n_sims <- 400; results <- list()
  withProgress(message = 'Running 400 Simulations...', value = 0, {
    for (i in 1:n_sims) {
      res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100)
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1/n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })
  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  True_AF <- input$sim_true_af; True_Med_WT <- 2.0; True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100
  calc_stats <- function(model_name, metric_name, true_val) {
    vals <- sapply(results, function(x) x[[model_name]][metric_name])
    mean_val <- mean(vals, na.rm = TRUE); mse_val <- mean((vals - true_val)^2, na.rm = TRUE)
    return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
  }

  af_l <- sapply(results, function(x) x$Prop["AF_X_L"]); af_u <- sapply(results, function(x) x$Prop["AF_X_U"])
  cp_val <- mean(af_l <= True_AF & True_AF <= af_u, na.rm = TRUE) * 100

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Median OS (WT)", "5-Year OS (WT) [%]"),
      `True Value` = c(sprintf("%.2f", True_AF), "0.85", "1.10", "0.90", "0.95", sprintf("%.2f", True_Med_WT), sprintf("%.1f", True_S5_WT)),
      `Naive AFT` = c(calc_stats("Naive", "AF_X", True_AF), calc_stats("Naive", "AF_Age", 0.85), calc_stats("Naive", "AF_Sex", 1.10), calc_stats("Naive", "AF_READ", 0.90), calc_stats("Naive", "AF_COADREAD", 0.95), calc_stats("Naive", "Med_WT", True_Med_WT), calc_stats("Naive", "S5_WT", True_S5_WT)),
      `Standard LT AFT` = c(calc_stats("LT", "AF_X", True_AF), calc_stats("LT", "AF_Age", 0.85), calc_stats("LT", "AF_Sex", 1.10), calc_stats("LT", "AF_READ", 0.90), calc_stats("LT", "AF_COADREAD", 0.95), calc_stats("LT", "Med_WT", True_Med_WT), calc_stats("LT", "S5_WT", True_S5_WT)),
      `Proposed Method` = c(paste0("<b>", calc_stats("Prop", "AF_X", True_AF), " [CP: ", sprintf("%.1f", cp_val), "%]</b>"), sprintf("<b>%s</b>", calc_stats("Prop", "AF_Age", 0.85)), sprintf("<b>%s</b>", calc_stats("Prop", "AF_Sex", 1.10)), sprintf("<b>%s</b>", calc_stats("Prop", "AF_READ", 0.90)), sprintf("<b>%s</b>", calc_stats("Prop", "AF_COADREAD", 0.95)), sprintf("<b>%s</b>", calc_stats("Prop", "Med_WT", True_Med_WT)), sprintf("<b>%s</b>", calc_stats("Prop", "S5_WT", True_S5_WT)))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
