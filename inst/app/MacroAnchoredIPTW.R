# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.0 - Bug Fix)
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

# 先生ご提案の「死亡例OSマッチング IPTW」
calculate_iptw_v2 <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points <- 1:5

  for(ag in unique(data$Age_class)) {
    if(ag %in% names(ref_surv_list)) {
      S_t <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
      y <- log(1/S_t - 1)
      x <- log(t_points)
      fit <- lm(y ~ x)
      p <- coef(fit)[2]
      lambda <- exp(coef(fit)[1] / p)

      ag_data <- data[data$Age_class == ag, ]
      dead_data <- ag_data[ag_data$Event == 1, ]
      if(nrow(dead_data) < 5) next

      pdf_macro <- function(t_y) { (lambda^p * p * t_y^(p-1)) / (1 + (lambda * t_y)^p)^2 }
      p_macro_val <- pdf_macro(dead_data$T_obs / 365.25)

      dens_dead_os <- density(dead_data$T_obs / 365.25, n = 512, from = 0)
      f_dead_os <- approxfun(dens_dead_os$x, dens_dead_os$y, rule = 2)

      w_os <- p_macro_val / pmax(f_dead_os(dead_data$T_obs / 365.25), 1e-4)
      w_os <- w_os / sum(w_os)

      # suppressWarningsでdensityの帯域幅計算警告を消去
      dens_t1_target <- suppressWarnings(density(dead_data$T1 / 365.25, weights = w_os, n = 512, from = 0))
      f_t1_target <- approxfun(dens_t1_target$x, dens_t1_target$y, rule = 2)

      dens_t1_all <- density(ag_data$T1 / 365.25, n = 512, from = 0)
      f_t1_all <- approxfun(dens_t1_all$x, dens_t1_all$y, rule = 2)

      w_final <- f_t1_target(ag_data$T1 / 365.25) / pmax(f_t1_all(ag_data$T1 / 365.25), 1e-4)
      data$iptw[data$Age_class == ag] <- w_final
    }
  }

  lower_bound <- quantile(data$iptw, 0.05, na.rm = TRUE)
  upper_bound <- quantile(data$iptw, 0.95, na.rm = TRUE)
  data$iptw <- pmax(lower_bound, pmin(data$iptw, upper_bound))
  data$iptw <- data$iptw / mean(data$iptw, na.rm = TRUE)

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

  ref_surv_list <- list()
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[ifelse(Age < 60, "Young", "Old") == ag] > y * 365.25) * 100)
  }
  Data_cgp <- calculate_iptw_v2(Data_cgp, ref_surv_list)

  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")
  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(valid_covs, collapse=" + ")))

  med_t1 <- median(Data_cgp$T1)
  data_short <- subset(Data_cgp, T1 <= med_t1)
  data_long <- subset(Data_cgp, T1 > med_t1)

  fit_naive <- tryCatch({ flexsurvreg(update(Surv(T_obs, Event) ~ ., as.formula(paste("~", paste(valid_covs, collapse=" + ")))), data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  fit_t2_short <- tryCatch({ flexsurvreg(form_t2, data = data_short, weights = iptw, dist = "llogis") }, error = function(e) NULL)
  fit_t2_long <- tryCatch({ flexsurvreg(form_t2, data = data_long, weights = iptw, dist = "llogis") }, error = function(e) NULL)

  extract_af <- function(rname) {
    af_s <- if(!is.null(fit_t2_short) && rname %in% rownames(fit_t2_short$res)) fit_t2_short$res[rname, "est"] else NA
    af_l <- if(!is.null(fit_t2_long) && rname %in% rownames(fit_t2_long$res)) fit_t2_long$res[rname, "est"] else NA

    if(!is.na(af_s) && !is.na(af_l)) { return(exp((af_s * nrow(data_short) + af_l * nrow(data_long)) / nrow(Data_cgp)))
    } else if (!is.na(af_s)) { return(exp(af_s))
    } else if (!is.na(af_l)) { return(exp(af_l))
    } else { return(NA) }
  }

  af_res <- list(
    AF_X = extract_af("XMutated"),
    AF_Age = extract_af("Age")^10,
    AF_Sex = extract_af("SexFemale"),
    AF_READ = extract_af("HistologyREAD"),
    AF_COADREAD = extract_af("HistologyCOADREAD")
  )

  out <- list(Naive = fit_naive, Prop_Short = fit_t2_short, Prop_Long = fit_t2_long,
              Metrics = af_res, Data_cgp = Data_cgp, med_t1 = med_t1, T_true_macro = T_true)
  if(return_data) return(out) else return(out$Metrics)
}

safe_fmt <- function(x, digits=2) { sapply(x, function(v) if(is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v)) }

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation (Tamura & Ikegami Ver 2.0)...", type = "message", duration = 3)

  res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE)
  shiny::validate(shiny::need(!is.null(res), "Simulation failed."))

  True_AF <- input$sim_true_af

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90"),
      `Proposed Method` = c(
        paste0("<b>", safe_fmt(res$Metrics$AF_X), "</b>"), paste0("<b>", safe_fmt(res$Metrics$AF_Age), "</b>"),
        paste0("<b>", safe_fmt(res$Metrics$AF_Sex), "</b>"), paste0("<b>", safe_fmt(res$Metrics$AF_READ), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    km_true <- survfit(Surv(res$T_true_macro, rep(1, length(res$T_true_macro))) ~ 1)
    surv_true <- approx(km_true$time, km_true$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y

    surv_naive <- rep(NA, length(t_seq))
    if(!is.null(res$Naive)) {
      summ_n <- summary(res$Naive, newdata = res$Data_cgp[sample(1:nrow(res$Data_cgp), min(300, nrow(res$Data_cgp))), ], t = t_seq, tidy = TRUE)
      surv_naive <- aggregate(est ~ time, data = summ_n, FUN = mean)$est
    }

    idx_pseudo <- sample(1:nrow(res$Data_cgp), size = 3000, replace = TRUE, prob = res$Data_cgp$iptw)
    pseudo_cgp <- res$Data_cgp[idx_pseudo, ]

    gen_t2 <- function(fit, newdata) {
      if(is.null(fit)) return(rep(NA, nrow(newdata)))
      res_m <- fit$res
      base_scale <- res_m["scale", "est"]; shape <- res_m["shape", "est"]

      af_vec <- rep(1, nrow(newdata))
      if("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * (newdata$X == "Mutated"))
      if("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
      if("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * (newdata$Sex == "Female"))
      if("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * (newdata$Histology == "READ"))
      if("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * (newdata$Histology == "COADREAD"))

      return(flexsurv::rllogis(nrow(newdata), shape = shape, scale = base_scale * af_vec))
    }

    pseudo_cgp$sim_T2 <- ifelse(pseudo_cgp$T1 <= res$med_t1, gen_t2(res$Prop_Short, pseudo_cgp), gen_t2(res$Prop_Long, pseudo_cgp))
    pseudo_cgp$sim_OS <- pseudo_cgp$T1 + pseudo_cgp$sim_T2

    km_prop <- survfit(Surv(sim_OS, rep(1, nrow(pseudo_cgp))) ~ 1, data = pseudo_cgp)
    surv_prop <- approx(km_prop$time, km_prop$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y

    plot_df <- data.frame(
      Time = rep(t_seq / 365.25, 3),
      Survival = c(surv_true, surv_naive, surv_prop),
      Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT (Immortal Time Bias)", "3. Proposed Method (Ver 2.0: IPTW + T1-Stratified G-comp)"), each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Ver 2.0 Validation",
           subtitle = "Stratified T2 simulation completely neutralizing dependent truncation & early censoring.",
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

  True_AF <- input$sim_true_af
  calc_stats <- function(metric_name, true_val) {
    vals <- sapply(results, function(x) x[[metric_name]])
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    mse_val <- mean((vals - true_val)^2)
    return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
  }

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95"),
      `Proposed Method` = c(
        paste0("<b>", calc_stats("AF_X", True_AF), "</b>"), paste0("<b>", calc_stats("AF_Age", 0.85), "</b>"),
        paste0("<b>", calc_stats("AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats("AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats("AF_COADREAD", 0.95), "</b>")
      )
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
