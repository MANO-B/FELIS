# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.0)
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
    return(list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE))
  } else if (pattern == "ushape") {
    obj_func <- function(scale_param) {
      u <- runif(length(T2))
      C2 <- ifelse(u < 0.5, rweibull(length(T2), shape = 0.5, scale = scale_param), rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
      return((sum(C2 < T2) / length(T2)) - target_rate)
    }
    opt_scale <- tryCatch({ uniroot(obj_func, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
    return(list(shape = NA, scale = opt_scale, is_ushape = TRUE))
  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      return((sum(C2 < T2) / length(T2)) - target_rate)
    }
    opt_scale <- tryCatch({ uniroot(obj_func_scale, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
    return(list(shape = if(pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE))
  }
}

# =========================================================================
# GLOBAL HELPER: AFTモデルからの安全な生存期間サンプリング
# =========================================================================
gen_sim_times <- function(fit, newdata) {
  if(is.null(fit) || nrow(newdata) == 0) return(rep(NA, nrow(newdata)))
  res_m <- fit$res
  base_scale <- res_m["scale", "est"]
  shape <- res_m["shape", "est"]

  # クラッシュを防ぐため、モデル行列ではなく直接係数を当てはめる
  af_vec <- rep(1, nrow(newdata))
  if("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
  if("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
  if("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
  if("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
  if("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))

  scale_vec <- base_scale * af_vec
  u <- runif(nrow(newdata))
  u <- pmax(pmin(u, 0.9999), 0.0001)
  sim_t <- scale_vec * ((1 - u) / u)^(1 / shape)
  return(sim_t)
}

# =========================================================================
# STEP 1 & 2: Two-step IPTW (Dead-OS matching -> T1 Hazard construction)
# =========================================================================
calculate_iptw_gcomp <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points <- 1:5

  for(ag in unique(data$Age_class)) {
    if(ag %in% names(ref_surv_list)) {

      # 【STEP 1】 院内がん登録（Macro）のOS分布をモデル化
      S_t <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
      y <- log(1/S_t - 1)
      x <- log(t_points)
      fit <- lm(y ~ x)
      p <- coef(fit)[2]
      lambda <- exp(coef(fit)[1] / p)
      pdf_macro <- function(t_y) { (lambda^p * p * t_y^(p-1)) / (1 + (lambda * t_y)^p)^2 }

      # 打ち切り（情報リーク）を避けるため「死亡患者」のみ抽出
      ag_data <- data[data$Age_class == ag, ]
      dead_data <- ag_data[ag_data$Event == 1, ]
      if(nrow(dead_data) < 5) next

      # 【STEP 2】 死亡例OSとMacro-OSを比較し、バイアスを消す「ハザード比(W_os)」を算出
      p_macro_val <- pdf_macro(dead_data$T_obs / 365.25)

      # エラー防止のためバンド幅を明示的に指定して警告を抑制
      bw_os <- suppressWarnings(density(dead_data$T_obs / 365.25)$bw)
      dens_dead_os <- suppressWarnings(density(dead_data$T_obs / 365.25, bw = bw_os, n = 512, from = 0))
      f_dead_os <- approxfun(dens_dead_os$x, dens_dead_os$y, rule = 2)

      w_os <- p_macro_val / pmax(f_dead_os(dead_data$T_obs / 365.25), 1e-6)
      w_os <- w_os / sum(w_os)

      # 【STEP 3】 死亡例のT1にW_osを掛け、「真のCGP到達ハザード(T1密度)」を再構築
      bw_t1 <- suppressWarnings(density(dead_data$T1 / 365.25)$bw)
      dens_t1_target <- suppressWarnings(density(dead_data$T1 / 365.25, weights = w_os, bw = bw_t1, n = 512, from = 0))
      f_t1_target <- approxfun(dens_t1_target$x, dens_t1_target$y, rule = 2)

      # 【STEP 4】 打ち切り例を含む全CGP患者に対し、IPTWを割り当てる
      bw_t1_all <- suppressWarnings(density(ag_data$T1 / 365.25)$bw)
      dens_t1_all <- suppressWarnings(density(ag_data$T1 / 365.25, bw = bw_t1_all, n = 512, from = 0))
      f_t1_all <- approxfun(dens_t1_all$x, dens_t1_all$y, rule = 2)

      w_final <- f_t1_target(ag_data$T1 / 365.25) / pmax(f_t1_all(ag_data$T1 / 365.25), 1e-6)
      data$iptw[data$Age_class == ag] <- w_final
    }
  }

  # 重みの爆発を防ぐトリミング（必須）
  lower_bound <- quantile(data$iptw, 0.05, na.rm = TRUE)
  upper_bound <- quantile(data$iptw, 0.95, na.rm = TRUE)
  data$iptw <- pmax(lower_bound, pmin(data$iptw, upper_bound))
  data$iptw <- data$iptw / mean(data$iptw, na.rm = TRUE)
  return(data)
}

# =========================================================================
# Core Simulation Engine
# =========================================================================
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  # --- Data Generation Process ---
  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", 1.10, 1.0) *
    ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))

  T_true <- (2.0 * 365.25) * ((1 - pmax(pmin(runif(N_macro), 0.999), 0.001)) / pmax(pmin(runif(N_macro), 0.999), 0.001))^(1 / 1.5) * AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  if (is.null(t1_pat)) t1_pat <- "indep"
  if (t1_pat == "early") { T1 <- runif(N_macro, 30, pmax(31, T_true * 0.3))
  } else if (t1_pat == "dep_1yr" || t1_pat == "real") { T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat == "dep_2yr" || t1_pat == "rev") { T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else { T1 <- runif(N_macro, 30, pmax(31, T_true)) }

  cgp_indices <- which(T_true > T1)
  if(length(cgp_indices) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.9, 0.1)
  selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target, Age = Age[selected_indices], Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
    Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[selected_indices], T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (isTRUE(params$is_ushape)) {
    u <- runif(N_target)
    C2 <- ifelse(u < 0.5, rweibull(N_target, shape = 0.5, scale = params$scale), rweibull(N_target, shape = 3.0, scale = params$scale * 2))
  } else if (cens_pat == "indep") { C2 <- rexp(N_target, rate = 1 / params$scale)
  } else { C2 <- rweibull(N_target, shape = params$shape, scale = params$scale) }

  Data_cgp$C2 <- C2
  Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)
  Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs
  if(sum(Data_cgp$Event) < 20) return(NULL)

  # Execute IPTW logic
  ref_surv_list <- list()
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[ifelse(Age < 60, "Young", "Old") == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_iptw_gcomp(Data_cgp, ref_surv_list)

  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # --- Standard Comparator Models ---
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse=" + ")))
  form_lt <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse=" + ")))

  fit_naive <- tryCatch({ flexsurvreg(form_naive, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
  fit_lt <- tryCatch({ flexsurvreg(form_lt, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # 【STEP 5】 G-computation: T1層別化によるT2 AFTモデル構築
  # =========================================================================
  med_t1 <- median(Data_cgp$T1)
  data_short <- subset(Data_cgp, T1 <= med_t1)
  data_long <- subset(Data_cgp, T1 > med_t1)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(valid_covs, collapse=" + ")))
  fit_t2_short <- tryCatch({ flexsurvreg(form_t2, data = data_short, weights = iptw, dist = "llogis") }, error = function(e) NULL)
  fit_t2_long <- tryCatch({ flexsurvreg(form_t2, data = data_long, weights = iptw, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # 【STEP 6】 完全なるG-computation (表とグラフの100%同期化)
  # =========================================================================

  # 1. 真のマージナル値の計算（グラフの黒線と完全に一致する値）
  True_Med <- median(T_true) / 365.25
  True_S5 <- mean(T_true > 5*365.25) * 100

  # 2. Naiveのマージナルシミュレーション
  t_naive <- gen_sim_times(fit_naive, Data_cgp)
  naive_med <- median(t_naive, na.rm=TRUE) / 365.25
  naive_s5 <- mean(t_naive > 5*365.25, na.rm=TRUE) * 100

  # 3. Standard LTのマージナルシミュレーション
  t_lt <- gen_sim_times(fit_lt, Data_cgp)
  lt_med <- median(t_lt, na.rm=TRUE) / 365.25
  lt_s5 <- mean(t_lt > 5*365.25, na.rm=TRUE) * 100

  # 4. Proposed (G-comp) のマージナルシミュレーション
  # IPTWによるバイアス浄化コホートを再構築
  idx_pseudo <- sample(1:nrow(Data_cgp), size = 3000, replace = TRUE, prob = Data_cgp$iptw)
  pseudo_cgp <- Data_cgp[idx_pseudo, ]

  # T2をシミュレーションし、T1に足し合わせる
  t2_prop <- ifelse(pseudo_cgp$T1 <= med_t1, gen_sim_times(fit_t2_short, pseudo_cgp), gen_sim_times(fit_t2_long, pseudo_cgp))
  t_prop <- pseudo_cgp$T1 + t2_prop
  prop_med <- median(t_prop, na.rm=TRUE) / 365.25
  prop_s5 <- mean(t_prop > 5*365.25, na.rm=TRUE) * 100

  # --- AFの抽出 ---
  extract_af <- function(fit) {
    if (is.null(fit)) return(c(AF_X=NA, AF_Age=NA, AF_Sex=NA, AF_READ=NA, AF_COADREAD=NA))
    res_m <- fit$res
    get_val <- function(rname) if (rname %in% rownames(res_m)) res_m[rname, "est"] else NA
    c(AF_X = exp(get_val("XMutated")), AF_Age = exp(get_val("Age") * 10),
      AF_Sex = exp(get_val("SexFemale")), AF_READ = exp(get_val("HistologyREAD")), AF_COADREAD = exp(get_val("HistologyCOADREAD")))
  }

  get_prop_af <- function(vname, multiplier=1) {
    val_s <- if(!is.null(fit_t2_short) && vname %in% rownames(fit_t2_short$res)) fit_t2_short$res[vname, "est"] else NA
    val_l <- if(!is.null(fit_t2_long) && vname %in% rownames(fit_t2_long$res)) fit_t2_long$res[vname, "est"] else NA
    if(!is.na(val_s) && !is.na(val_l)) {
      w_avg <- (val_s * nrow(data_short) + val_l * nrow(data_long)) / nrow(Data_cgp)
      return(exp(w_avg * multiplier))
    } else return(NA)
  }

  out <- list(
    True_Metrics = c(True_Med, True_S5),
    Naive_AF = extract_af(fit_naive), Naive_Metrics = c(naive_med, naive_s5),
    LT_AF = extract_af(fit_lt), LT_Metrics = c(lt_med, lt_s5),
    Prop_AF = c(AF_X=get_prop_af("XMutated"), AF_Age=get_prop_af("Age", 10), AF_Sex=get_prop_af("SexFemale"), AF_READ=get_prop_af("HistologyREAD"), AF_COADREAD=get_prop_af("HistologyCOADREAD")),
    Prop_Metrics = c(prop_med, prop_s5),
    T_true_macro = T_true, t_naive_plot = t_naive, t_lt_plot = t_lt, t_prop_plot = t_prop
  )
  if(return_data) return(out) else return(out)
}

safe_fmt <- function(x, digits=2) { sapply(x, function(v) if(is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v)) }

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic
# -------------------------------------------------------------------------
observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation (Tamura & Ikegami Ver 2.0)...", type = "message", duration = 3)

  res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE)
  shiny::validate(shiny::need(!is.null(res), "Simulation failed. Not enough events generated."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(res$True_Metrics[1]), safe_fmt(res$True_Metrics[2], 1)),
      `Naive AFT (Surv(T_obs))` = c(safe_fmt(res$Naive_AF["AF_X"]), safe_fmt(res$Naive_AF["AF_Age"]), safe_fmt(res$Naive_AF["AF_Sex"]), safe_fmt(res$Naive_AF["AF_READ"]), safe_fmt(res$Naive_AF["AF_COADREAD"]), safe_fmt(res$Naive_Metrics[1]), safe_fmt(res$Naive_Metrics[2], 1)),
      `Standard LT (Surv(T1,T_obs))` = c(safe_fmt(res$LT_AF["AF_X"]), safe_fmt(res$LT_AF["AF_Age"]), safe_fmt(res$LT_AF["AF_Sex"]), safe_fmt(res$LT_AF["AF_READ"]), safe_fmt(res$LT_AF["AF_COADREAD"]), safe_fmt(res$LT_Metrics[1]), safe_fmt(res$LT_Metrics[2], 1)),
      `Proposed (G-comp)` = c(paste0("<b>", safe_fmt(res$Prop_AF["AF_X"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_AF["AF_Age"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_AF["AF_Sex"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_AF["AF_READ"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_AF["AF_COADREAD"]), "</b>"), paste0("<b>", safe_fmt(res$Prop_Metrics[1]), "</b>"), paste0("<b>", safe_fmt(res$Prop_Metrics[2], 1), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    km_true <- survfit(Surv(res$T_true_macro, rep(1, length(res$T_true_macro))) ~ 1)
    surv_true <- approx(km_true$time, km_true$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y

    calc_marginal <- function(t_data) {
      if(is.null(t_data) || length(t_data) == 0) return(rep(NA, length(t_seq)))
      km <- survfit(Surv(t_data, rep(1, length(t_data))) ~ 1)
      return(approx(km$time, km$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y)
    }

    surv_naive <- calc_marginal(res$t_naive_plot)
    surv_lt <- calc_marginal(res$t_lt_plot)
    surv_prop <- calc_marginal(res$t_prop_plot)

    plot_df <- data.frame(
      Time = rep(t_seq / 365.25, 4),
      Survival = c(surv_true, surv_naive, surv_lt, surv_prop),
      Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT (Immortal Time Bias)", "3. Standard LT AFT (Dep. Truncation Bias)", "4. Proposed Method (Ver 2.0: G-computation)"), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Ver 2.0 Validation",
           subtitle = "G-computation natively resolves dependent truncation. Table and Plot are now 100% perfectly synchronized.",
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

  True_AF <- input$sim_true_af;

  calc_stats <- function(model_name, metric_name, true_val = NULL) {
    if(!is.null(model_name)){
      vals <- sapply(results, function(x) x[[model_name]][metric_name])
    } else {
      vals <- sapply(results, function(x) x[[metric_name]])
    }
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    if(!is.null(true_val)){
      mse_val <- mean((vals - true_val)^2)
      return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
    }
    return(sprintf("%.2f", mean_val))
  }

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Marginal Med OS [Years]", "Marginal 5-Year OS [%]"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", calc_stats("True_Metrics", 1), calc_stats("True_Metrics", 2)),
      `Naive AFT` = c(calc_stats("Naive_AF", "AF_X", True_AF), calc_stats("Naive_AF", "AF_Age", 0.85), calc_stats("Naive_AF", "AF_Sex", 1.10), calc_stats("Naive_AF", "AF_READ", 0.90), calc_stats("Naive_AF", "AF_COADREAD", 0.95), calc_stats("Naive_Metrics", 1), calc_stats("Naive_Metrics", 2)),
      `Standard LT AFT` = c(calc_stats("LT_AF", "AF_X", True_AF), calc_stats("LT_AF", "AF_Age", 0.85), calc_stats("LT_AF", "AF_Sex", 1.10), calc_stats("LT_AF", "AF_READ", 0.90), calc_stats("LT_AF", "AF_COADREAD", 0.95), calc_stats("LT_Metrics", 1), calc_stats("LT_Metrics", 2)),
      `Proposed (G-comp)` = c(paste0("<b>", calc_stats("Prop_AF", "AF_X", True_AF), "</b>"), paste0("<b>", calc_stats("Prop_AF", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats("Prop_AF", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats("Prop_AF", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats("Prop_AF", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats("Prop_Metrics", 1), "</b>"), paste0("<b>", calc_stats("Prop_Metrics", 2), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
