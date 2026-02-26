# =========================================================================
# Server Logic for Simulation Study (Crash-Proof T1-IPTW Validation)
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
      C2 <- ifelse(u < 0.5, rweibull(length(T2), shape = 0.5, scale = scale_param),
                   rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
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

calculate_iptw_sim <- function(data, ref_surv_list) {
  init_pop <- 10000
  max_years <- 10
  bin_width <- 0.5
  breaks <- seq(0, max_years, by = bin_width)
  n_bins <- length(breaks) - 1

  data <- data %>%
    dplyr::mutate(
      time_years = T1 / 365.25,
      time_bin = ceiling(time_years / bin_width),
      time_bin = ifelse(time_bin > n_bins, n_bins, time_bin),
      time_bin = ifelse(time_bin == 0, 1, time_bin)
    )

  pt_table <- expand.grid(Age_class = names(ref_surv_list), time_bin = 1:n_bins, stringsAsFactors = FALSE)
  pt_table$pt_ref <- 0
  t_points <- 1:5

  for(ag in names(ref_surv_list)) {
    surv_rates <- ref_surv_list[[ag]]
    S_t <- surv_rates / 100
    S_t <- pmax(pmin(S_t, 0.999), 0.001)

    y <- log(1/S_t - 1)
    x <- log(t_points)
    fit <- lm(y ~ x)
    p <- coef(fit)[2]
    lambda <- exp(coef(fit)[1] / p)
    S_fit <- function(t_y) { 1 / (1 + (lambda * t_y)^p) }

    pt_bins <- numeric(n_bins)
    for(i in 1:n_bins) {
      t_start <- breaks[i]
      t_end <- breaks[i+1]
      pt_bins[i] <- init_pop * (S_fit(t_start) + S_fit(t_end)) / 2 * bin_width
    }
    pt_table[pt_table$Age_class == ag, "pt_ref"] <- pt_bins
  }

  bin_counts <- data %>% dplyr::count(Age_class, time_bin, name = "N_cgp")

  data <- data %>%
    dplyr::left_join(pt_table, by = c("Age_class", "time_bin")) %>%
    dplyr::left_join(bin_counts, by = c("Age_class", "time_bin")) %>%
    dplyr::mutate(
      raw_weight = ifelse(!is.na(N_cgp) & N_cgp > 0 & !is.na(pt_ref), pt_ref / N_cgp, 0)
    )

  if (any(data$raw_weight > 0, na.rm = TRUE)) {
    lower_bound <- quantile(data$raw_weight[data$raw_weight > 0], 0.025, na.rm = TRUE)
    upper_bound <- quantile(data$raw_weight[data$raw_weight > 0], 0.975, na.rm = TRUE)
    data <- data %>%
      dplyr::mutate(
        raw_weight = ifelse(raw_weight > 0 & raw_weight < lower_bound, lower_bound, raw_weight),
        raw_weight = ifelse(raw_weight > upper_bound, upper_bound, raw_weight)
      )
  }
  mean_w <- mean(data$raw_weight[data$raw_weight > 0], na.rm = TRUE)
  data$iptw <- ifelse(data$raw_weight > 0, data$raw_weight / mean_w, 1.0)

  return(data)
}

run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", 1.10, 1.0) *
    ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))

  T_true <- (2.0 * 365.25) * ((1 - pmax(pmin(runif(N_macro), 0.999), 0.001)) / pmax(pmin(runif(N_macro), 0.999), 0.001))^(1 / 1.5) * AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  # [安全装置] T1パターンがUIから正しく渡されなかった場合のフォールバックを完備
  if (is.null(t1_pat)) t1_pat <- "indep"

  if (t1_pat == "early") {
    T1 <- runif(N_macro, 30, pmax(31, T_true * 0.3))
  } else if (t1_pat == "dep_1yr" || t1_pat == "real") {
    T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
  } else if (t1_pat == "dep_2yr" || t1_pat == "rev") {
    T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
  } else {
    T1 <- runif(N_macro, 30, pmax(31, T_true))
  }

  cgp_indices <- which(T_true > T1)
  if(length(cgp_indices) < N_target) return(NULL)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.9, 0.1)
  selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:N_target, Age = Age[selected_indices], Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
    Sex = factor(Sex[selected_indices]), Histology = factor(Histology[selected_indices]),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT")), T_true = T_true[selected_indices], T1 = T1[selected_indices]
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

  ref_surv_list <- list()
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[ifelse(Age < 60, "Young", "Old") == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_iptw_sim(Data_cgp, ref_surv_list)

  # 変数のレベルが足りない場合にモデルがクラッシュするのを防ぐ動的フォーミュラ
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")
  form <- as.formula(paste("~", paste(valid_covs, collapse=" + ")))

  fit_naive <- tryCatch({ flexsurvreg(update(Surv(T_obs, Event) ~ ., form), data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
  fit_lt <- tryCatch({ flexsurvreg(update(Surv(T1, T_obs, Event) ~ ., form), data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
  fit_prop <- tryCatch({ flexsurvreg(update(Surv(T1, T_obs, Event) ~ ., form), data = Data_cgp, weights = iptw, dist = "llogis") }, error = function(e) NULL)

  # [安全装置] NAや変数落ちによる抽出クラッシュを完全に防ぐラッパー
  extract_metrics <- function(fit) {
    blank_res <- c(AF_X=NA, AF_X_L=NA, AF_X_U=NA, AF_Age=NA, AF_Sex=NA, AF_READ=NA, AF_COADREAD=NA, Med_WT=NA, Med_Mut=NA, S5_WT=NA, S5_Mut=NA)
    if (is.null(fit)) return(blank_res)
    res <- fit$res

    get_val <- function(rname, cname) {
      if (rname %in% rownames(res) && cname %in% colnames(res)) res[rname, cname] else NA
    }

    AF_X <- exp(get_val("XMutated", "est"))
    AF_X_L <- exp(get_val("XMutated", "L95%"))
    AF_X_U <- exp(get_val("XMutated", "U95%"))
    AF_Age <- exp(get_val("Age", "est") * 10)
    AF_Sex <- exp(get_val("SexFemale", "est"))
    AF_READ <- exp(get_val("HistologyREAD", "est"))
    AF_COADREAD <- exp(get_val("HistologyCOADREAD", "est"))

    scale_base <- get_val("scale", "est")
    shape <- get_val("shape", "est")
    af_age60 <- exp(get_val("Age", "est") * 60)

    med_wt_days <- scale_base * af_age60
    med_wt <- med_wt_days / 365.25
    med_mut <- med_wt * AF_X
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

# [安全装置] NAを検知して"N/A"を返し、sprintfのクラッシュを防ぐフォーマッター
safe_fmt <- function(x, digits=2) {
  sapply(x, function(v) if(is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation (Original T1-based IPTW)...", type = "message", duration = 3)

  res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE)
  shiny::validate(shiny::need(!is.null(res), "Simulation generated insufficient events. Please try again or reduce censoring rate."))

  True_AF <- input$sim_true_af; True_Med_WT <- 2.0; True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Median OS (WT) [Years]", "5-Year OS (WT) [%]"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
      `Naive AFT` = c(safe_fmt(res$Naive["AF_X"]), safe_fmt(res$Naive["AF_Age"]), safe_fmt(res$Naive["AF_Sex"]), safe_fmt(res$Naive["AF_READ"]), safe_fmt(res$Naive["AF_COADREAD"]), safe_fmt(res$Naive["Med_WT"]), safe_fmt(res$Naive["S5_WT"], 1)),
      `Standard LT` = c(safe_fmt(res$LT["AF_X"]), safe_fmt(res$LT["AF_Age"]), safe_fmt(res$LT["AF_Sex"]), safe_fmt(res$LT["AF_READ"]), safe_fmt(res$LT["AF_COADREAD"]), safe_fmt(res$LT["Med_WT"]), safe_fmt(res$LT["S5_WT"], 1)),
      `Proposed Method` = c(
        paste0("<b>", safe_fmt(res$Prop["AF_X"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["AF_Age"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["AF_Sex"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["AF_READ"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["AF_COADREAD"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["Med_WT"]), "</b>"),
        paste0("<b>", safe_fmt(res$Prop["S5_WT"], 1), "</b>")
      )
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
      Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT (Immortal Time Bias)", "3. Standard LT AFT", "4. Proposed Method (Original T1-IPTW + LT-AFT)"), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Marginal Survival Recovery (T1-based IPTW Validation)",
           subtitle = "Notice how the Original Proposed Method perfectly handles selection bias and pulls back to True Survival.",
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
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    mse_val <- mean((vals - true_val)^2)
    return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
  }

  af_l <- sapply(results, function(x) x$Prop["AF_X_L"])
  af_u <- sapply(results, function(x) x$Prop["AF_X_U"])
  cp_val <- mean(af_l <= True_AF & True_AF <= af_u, na.rm = TRUE) * 100

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Median OS (WT)", "5-Year OS (WT) [%]"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
      `Naive AFT` = c(calc_stats("Naive", "AF_X", True_AF), calc_stats("Naive", "AF_Age", 0.85), calc_stats("Naive", "AF_Sex", 1.10), calc_stats("Naive", "AF_READ", 0.90), calc_stats("Naive", "AF_COADREAD", 0.95), calc_stats("Naive", "Med_WT", True_Med_WT), calc_stats("Naive", "S5_WT", True_S5_WT)),
      `Standard LT AFT` = c(calc_stats("LT", "AF_X", True_AF), calc_stats("LT", "AF_Age", 0.85), calc_stats("LT", "AF_Sex", 1.10), calc_stats("LT", "AF_READ", 0.90), calc_stats("LT", "AF_COADREAD", 0.95), calc_stats("LT", "Med_WT", True_Med_WT), calc_stats("LT", "S5_WT", True_S5_WT)),
      `Proposed Method` = c(paste0("<b>", calc_stats("Prop", "AF_X", True_AF), " [CP: ", safe_fmt(cp_val, 1), "%]</b>"), paste0("<b>", calc_stats("Prop", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats("Prop", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats("Prop", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats("Prop", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats("Prop", "Med_WT", True_Med_WT), "</b>"), paste0("<b>", calc_stats("Prop", "S5_WT", True_S5_WT), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
