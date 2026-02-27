# # =========================================================================
# # Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.0)
# # =========================================================================
# library(dplyr)
# library(flexsurv)
# library(survival)
# library(ggplot2)
#
# # Helper function to generate complex censoring patterns
# find_censoring_param <- function(T2, target_rate, pattern) {
#   if (pattern == "peak1y") {
#     obj_func <- function(k) { scale_val <- 365.25 / ((k - 1) / k)^(1 / k); return((sum(rweibull(length(T2), shape = k, scale = scale_val) < T2) / length(T2)) - target_rate) }
#     opt_k <- tryCatch({ uniroot(obj_func, interval = c(1.05, 10))$root }, error = function(e) { 2.0 })
#     return(list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE))
#   } else if (pattern == "ushape") {
#     obj_func <- function(scale_param) {
#       u <- runif(length(T2))
#       C2 <- ifelse(u < 0.5, rweibull(length(T2), shape = 0.5, scale = scale_param), rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
#       return((sum(C2 < T2) / length(T2)) - target_rate)
#     }
#     opt_scale <- tryCatch({ uniroot(obj_func, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
#     return(list(shape = NA, scale = opt_scale, is_ushape = TRUE))
#   } else {
#     obj_func_scale <- function(scale_param) {
#       C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
#       return((sum(C2 < T2) / length(T2)) - target_rate)
#     }
#     opt_scale <- tryCatch({ uniroot(obj_func_scale, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
#     return(list(shape = if(pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE))
#   }
# }
#
# # =========================================================================
# # GLOBAL HELPERS: AFTモデルからの安全な生存期間サンプリング関数
# # =========================================================================
#
# # 1. Weibull分布用サンプリング関数 (T1生成用)
# gen_sim_times_weibull <- function(fit, newdata) {
#   if(is.null(fit) || nrow(newdata) == 0) return(rep(NA, nrow(newdata)))
#   res_m <- fit$res
#   base_scale <- res_m["scale", "est"]
#   shape <- res_m["shape", "est"]
#
#   af_vec <- rep(1, nrow(newdata))
#   if("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
#   if("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
#   if("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
#   if("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
#   if("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))
#
#   scale_vec <- base_scale * af_vec
#   u <- runif(nrow(newdata))
#   u <- pmax(pmin(u, 0.9999), 0.0001)
#   # Weibullの逆関数によるサンプリング: t = scale * (-log(u))^(1/shape)
#   sim_t <- scale_vec * (-log(u))^(1 / shape)
#   return(sim_t)
# }
#
# # 2. llogis分布用サンプリング関数 (T2生成・最終モデル用)
# gen_sim_times_llogis <- function(fit, newdata) {
#   if(is.null(fit) || nrow(newdata) == 0) return(rep(NA, nrow(newdata)))
#   res_m <- fit$res
#   base_scale <- res_m["scale", "est"]
#   shape <- res_m["shape", "est"]
#
#   af_vec <- rep(1, nrow(newdata))
#   if("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
#   if("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
#   if("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
#   if("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
#   if("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))
#
#   scale_vec <- base_scale * af_vec
#   u <- runif(nrow(newdata))
#   u <- pmax(pmin(u, 0.9999), 0.0001)
#   # llogisの逆関数によるサンプリング: t = scale * ((1-u)/u)^(1/shape)
#   sim_t <- scale_vec * ((1 - u) / u)^(1 / shape)
#   return(sim_t)
# }
#
# # =========================================================================
# # 【STEP 1〜3】 Two-step IPTW (Dead-OS matching -> True T1 Hazard construction)
# # =========================================================================
# calculate_iptw_gcomp <- function(data, ref_surv_list) {
#   data$iptw <- 1.0
#   t_points <- 1:5
#
#   for(ag in unique(data$Age_class)) {
#     if(ag %in% names(ref_surv_list)) {
#
#       # 【STEP 1】 院内がん登録（Macro）のOS分布をモデル化 (llogis)
#       S_t <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
#       y <- log(1/S_t - 1)
#       x <- log(t_points)
#       fit <- lm(y ~ x)
#       p <- coef(fit)[2]
#       lambda <- exp(coef(fit)[1] / p)
#       pdf_macro <- function(t_y) { (lambda^p * p * t_y^(p-1)) / (1 + (lambda * t_y)^p)^2 }
#
#       # 推奨治療非実施かつ「死亡例(Event=1)」のみ抽出
#       ag_data <- data[data$Age_class == ag, ]
#       dead_data <- ag_data[ag_data$Event == 1, ]
#       if(nrow(dead_data) < 5) next
#
#       # 【STEP 2】 死亡例のOS分布を院内がん登録に一致させる重み(W_os)を算出
#       p_macro_val <- pdf_macro(dead_data$T_obs / 365.25)
#       bw_os <- suppressWarnings(density(dead_data$T_obs / 365.25)$bw)
#       dens_dead_os <- suppressWarnings(density(dead_data$T_obs / 365.25, bw = bw_os, n = 512, from = 0))
#       f_dead_os <- approxfun(dens_dead_os$x, dens_dead_os$y, rule = 2)
#
#       w_os <- p_macro_val / pmax(f_dead_os(dead_data$T_obs / 365.25), 1e-6)
#       w_os <- w_os / sum(w_os)
#
#       # 【STEP 3】 W_osを用いて「真のCGP到達ハザード(T1密度)」を再構築し、全患者に適用
#       bw_t1 <- suppressWarnings(density(dead_data$T1 / 365.25)$bw)
#       dens_t1_target <- suppressWarnings(density(dead_data$T1 / 365.25, weights = w_os, bw = bw_t1, n = 512, from = 0))
#       f_t1_target <- approxfun(dens_t1_target$x, dens_t1_target$y, rule = 2)
#
#       bw_t1_all <- suppressWarnings(density(ag_data$T1 / 365.25)$bw)
#       dens_t1_all <- suppressWarnings(density(ag_data$T1 / 365.25, bw = bw_t1_all, n = 512, from = 0))
#       f_t1_all <- approxfun(dens_t1_all$x, dens_t1_all$y, rule = 2)
#
#       w_final <- f_t1_target(ag_data$T1 / 365.25) / pmax(f_t1_all(ag_data$T1 / 365.25), 1e-6)
#       data$iptw[data$Age_class == ag] <- w_final
#     }
#   }
#
#   lower_bound <- quantile(data$iptw, 0.05, na.rm = TRUE)
#   upper_bound <- quantile(data$iptw, 0.95, na.rm = TRUE)
#   data$iptw <- pmax(lower_bound, pmin(data$iptw, upper_bound))
#   data$iptw <- data$iptw / mean(data$iptw, na.rm = TRUE)
#   return(data)
# }
#
# # =========================================================================
# # Core Simulation Engine
# # =========================================================================
# run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {
#
#   # --- Data Generation Process ---
#   N_macro <- 100000
#   Age <- round(runif(N_macro, 40, 80))
#   Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
#   Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
#   X <- rbinom(N_macro, 1, Mut_Freq)
#
#   AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) * ifelse(Sex == "Female", 1.10, 1.0) *
#     ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))
#
#   T_true <- (2.0 * 365.25) * ((1 - pmax(pmin(runif(N_macro), 0.999), 0.001)) / pmax(pmin(runif(N_macro), 0.999), 0.001))^(1 / 1.5) * AF_bg * ifelse(X == 1, True_AF_X, 1.0)
#
#   if (is.null(t1_pat)) t1_pat <- "indep"
#   if (t1_pat == "early") { T1 <- runif(N_macro, 30, pmax(31, T_true * 0.3))
#   } else if (t1_pat == "dep_1yr" || t1_pat == "real") { T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 1.0), 0.4))
#   } else if (t1_pat == "dep_2yr" || t1_pat == "rev") { T1 <- pmax(30, T_true - rlnorm(N_macro, log(365.25 * 2.0), 0.4))
#   } else { T1 <- runif(N_macro, 30, pmax(31, T_true)) }
#
#   cgp_indices <- which(T_true > T1)
#   if(length(cgp_indices) < N_target) return(NULL)
#
#   prob_select <- ifelse(Age[cgp_indices] < 60, 0.9, 0.1)
#   selected_indices <- sample(cgp_indices, size = N_target, prob = prob_select, replace = FALSE)
#
#   Data_cgp <- data.frame(
#     ID = 1:N_target, Age = Age[selected_indices], Age_class = ifelse(Age[selected_indices] < 60, "Young", "Old"),
#     Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
#     Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
#     X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
#     T_true = T_true[selected_indices], T1 = T1[selected_indices]
#   )
#   Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1
#
#   params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
#   if (isTRUE(params$is_ushape)) {
#     u <- runif(N_target)
#     C2 <- ifelse(u < 0.5, rweibull(N_target, shape = 0.5, scale = params$scale), rweibull(N_target, shape = 3.0, scale = params$scale * 2))
#   } else if (cens_pat == "indep") { C2 <- rexp(N_target, rate = 1 / params$scale)
#   } else { C2 <- rweibull(N_target, shape = params$shape, scale = params$scale) }
#
#   Data_cgp$C2 <- C2
#   Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)
#   Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
#   Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs
#   if(sum(Data_cgp$Event) < 20) return(NULL)
#
#   # IPTW実行
#   ref_surv_list <- list()
#   for (ag in c("Young", "Old")) {
#     ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[ifelse(Age < 60, "Young", "Old") == ag] > y * 365.25) * 100)
#   }
#   Data_cgp <- calculate_iptw_gcomp(Data_cgp, ref_surv_list)
#
#   valid_covs <- c("X", "Age")
#   if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
#   if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")
#
#   # --- Standard Comparator Models ---
#   form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse=" + ")))
#   form_lt <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse=" + ")))
#
#   fit_naive <- tryCatch({ flexsurvreg(form_naive, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
#   fit_lt <- tryCatch({ flexsurvreg(form_lt, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
#
#   # =========================================================================
#   # 【STEP 4】 T1のモデルを Weibull で構築
#   # =========================================================================
#   Data_cgp$T1_event <- 1
#   form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse=" + ")))
#   fit_t1 <- tryCatch({ flexsurvreg(form_t1, data = Data_cgp, weights = iptw, dist = "weibull") }, error = function(e) NULL)
#
#   # =========================================================================
#   # 【STEP 5】 T1中央値で層別化した T2 モデルを llogis で構築
#   # 先生の指定通り、長期生存の裾野(Long-tail)を表現するため llogis を採用
#   # =========================================================================
#   med_t1 <- median(Data_cgp$T1)
#   data_short <- subset(Data_cgp, T1 <= med_t1)
#   data_long <- subset(Data_cgp, T1 > med_t1)
#
#   form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(valid_covs, collapse=" + ")))
#   fit_t2_short <- tryCatch({ flexsurvreg(form_t2, data = data_short, weights = iptw, dist = "llogis") }, error = function(e) NULL)
#   fit_t2_long <- tryCatch({ flexsurvreg(form_t2, data = data_long, weights = iptw, dist = "llogis") }, error = function(e) NULL)
#
#   # =========================================================================
#   # 【STEP 6】 G-computation: T1(weibull) + T2(llogis) のサンプリングと合算
#   # =========================================================================
#   idx_pseudo <- sample(1:nrow(Data_cgp), size = 10000, replace = TRUE, prob = Data_cgp$iptw)
#   pseudo_cgp <- Data_cgp[idx_pseudo, ]
#
#   # T1はWeibullでサンプリング
#   pseudo_cgp$sim_T1 <- gen_sim_times_weibull(fit_t1, pseudo_cgp)
#
#   # T2はllogisでサンプリング
#   pseudo_cgp$sim_T2 <- ifelse(pseudo_cgp$sim_T1 <= med_t1,
#                               gen_sim_times_llogis(fit_t2_short, pseudo_cgp),
#                               gen_sim_times_llogis(fit_t2_long, pseudo_cgp))
#
#   pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
#   pseudo_cgp$sim_Event <- 1
#
#   # 生成された完全OSに対して最終モデルを当てはめ、表とグラフを完全に同期させる
#   form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse=" + ")))
#   fit_prop_final <- tryCatch({ flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis") }, error = function(e) NULL)
#
#
#   # --- 指標抽出関数 ---
#   extract_metrics <- function(fit) {
#     blank_res <- c(AF_X=NA, AF_Age=NA, AF_Sex=NA, AF_READ=NA, AF_COADREAD=NA, Med_WT=NA, S5_WT=NA)
#     if (is.null(fit)) return(blank_res)
#     res_m <- fit$res
#     get_val <- function(rname) if (rname %in% rownames(res_m)) res_m[rname, "est"] else NA
#
#     AF_X <- exp(get_val("XMutated")); AF_Age <- exp(get_val("Age") * 10)
#     AF_Sex <- exp(get_val("SexFemale")); AF_READ <- exp(get_val("HistologyREAD"))
#     AF_COADREAD <- exp(get_val("HistologyCOADREAD"))
#
#     scale_base <- get_val("scale"); shape <- get_val("shape")
#     scale_wt_baseline <- scale_base * exp(get_val("Age") * 60)
#
#     med_wt <- scale_wt_baseline / 365.25
#     s5_wt <- 1 / (1 + (5 * 365.25 / scale_wt_baseline)^shape) * 100
#
#     return(c(AF_X=AF_X, AF_Age=AF_Age, AF_Sex=AF_Sex, AF_READ=AF_READ, AF_COADREAD=AF_COADREAD, Med_WT=med_wt, S5_WT=s5_wt))
#   }
#
#   out <- list(Naive = extract_metrics(fit_naive),
#               LT = extract_metrics(fit_lt),
#               Prop = extract_metrics(fit_prop_final),
#               fit_naive = fit_naive, fit_lt = fit_lt, fit_prop_final = fit_prop_final,
#               pseudo_cgp = pseudo_cgp, T_true_macro = T_true)
#
#   if(return_data) return(out) else return(list(Naive=out$Naive, LT=out$LT, Prop=out$Prop))
# }
#
# safe_fmt <- function(x, digits=2) { sapply(x, function(v) if(is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v)) }
#
# # -------------------------------------------------------------------------
# # UI Observers and Rendering Logic
# # -------------------------------------------------------------------------
# observeEvent(input$run_sim, {
#   req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
#   showNotification("Running simulation (Tamura & Ikegami Ver 2.0)...", type = "message", duration = 3)
#
#   res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100, return_data = TRUE)
#   shiny::validate(shiny::need(!is.null(res), "Simulation failed. Not enough events generated."))
#
#   True_AF <- input$sim_true_af; True_Med_WT <- 2.0; True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100
#
#   output$sim_result_table <- renderTable({
#     data.frame(
#       Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Baseline Med OS [Years]", "Baseline 5-Year OS [%]"),
#       `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
#       `Naive AFT (Surv(T_obs))` = c(safe_fmt(res$Naive["AF_X"]), safe_fmt(res$Naive["AF_Age"]), safe_fmt(res$Naive["AF_Sex"]), safe_fmt(res$Naive["AF_READ"]), safe_fmt(res$Naive["AF_COADREAD"]), safe_fmt(res$Naive["Med_WT"]), safe_fmt(res$Naive["S5_WT"], 1)),
#       `Standard LT (Surv(T1,T_obs))` = c(safe_fmt(res$LT["AF_X"]), safe_fmt(res$LT["AF_Age"]), safe_fmt(res$LT["AF_Sex"]), safe_fmt(res$LT["AF_READ"]), safe_fmt(res$LT["AF_COADREAD"]), safe_fmt(res$LT["Med_WT"]), safe_fmt(res$LT["S5_WT"], 1)),
#       `Proposed (G-comp)` = c(paste0("<b>", safe_fmt(res$Prop["AF_X"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["AF_Age"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["AF_Sex"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["AF_READ"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["AF_COADREAD"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["Med_WT"]), "</b>"), paste0("<b>", safe_fmt(res$Prop["S5_WT"], 1), "</b>"))
#     )
#   }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
#
#   output$sim_survival_plot <- renderPlot({
#     t_seq <- seq(0, 365.25 * 5, length.out = 100)
#
#     km_true <- survfit(Surv(res$T_true_macro, rep(1, length(res$T_true_macro))) ~ 1)
#     surv_true <- approx(km_true$time, km_true$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y
#
#     calc_marginal <- function(fit, dataset) {
#       if(is.null(fit)) return(rep(NA, length(t_seq)))
#       samp <- dataset[sample(1:nrow(dataset), min(1000, nrow(dataset))), ]
#       summ <- summary(fit, newdata = samp, t = t_seq, tidy = TRUE)
#       return(aggregate(est ~ time, data = summ, FUN = mean)$est)
#     }
#
#     surv_naive <- calc_marginal(res$fit_naive, res$pseudo_cgp)
#     surv_lt <- calc_marginal(res$fit_lt, res$pseudo_cgp)
#     surv_prop <- calc_marginal(res$fit_prop_final, res$pseudo_cgp)
#
#     plot_df <- data.frame(
#       Time = rep(t_seq / 365.25, 4),
#       Survival = c(surv_true, surv_naive, surv_lt, surv_prop),
#       Model = factor(rep(c("1. True Marginal (Registry Macro)", "2. Naive AFT (Immortal Time Bias)", "3. Standard LT AFT (Dep. Truncation Bias)", "4. Proposed Method (Ver 2.0: G-computation)"), each = length(t_seq)))
#     )
#
#     ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
#       geom_line(size = 1.2) +
#       scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
#       scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
#       scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
#       theme_minimal(base_size = 15) +
#       labs(title = "Tamura & Ikegami Model Ver 2.0 Validation",
#            subtitle = "G-computation resolving dependent truncation (T1=Weibull, T2=Llogis).",
#            x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
#       theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())
#   })
# })
#
# observeEvent(input$run_sim_multi, {
#   req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
#   n_sims <- 400; results <- list()
#   withProgress(message = 'Running 400 Simulations...', value = 0, {
#     for (i in 1:n_sims) {
#       res <- run_sim_iteration(input$sim_n, input$sim_true_af, input$sim_mut_freq / 100, input$sim_t1_pattern, input$sim_cens_pattern, input$sim_cens_rate / 100)
#       if (!is.null(res)) results[[length(results) + 1]] <- res
#       incProgress(1/n_sims, detail = paste("Iteration", i, "of", n_sims))
#     }
#   })
#   shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))
#
#   True_AF <- input$sim_true_af; True_Med_WT <- 2.0; True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100
#
#   calc_stats <- function(model_name, metric_name, true_val) {
#     vals <- sapply(results, function(x) x[[model_name]][metric_name])
#     vals <- vals[!is.na(vals)]
#     if(length(vals) == 0) return("N/A")
#     mean_val <- mean(vals)
#     mse_val <- mean((vals - true_val)^2)
#     return(sprintf("%.2f (MSE: %.3f)", mean_val, mse_val))
#   }
#
#   output$sim_multi_result_table <- renderTable({
#     data.frame(
#       Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF", "Baseline Med OS", "Baseline 5-Year OS [%]"),
#       `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
#       `Naive AFT` = c(calc_stats("Naive", "AF_X", True_AF), calc_stats("Naive", "AF_Age", 0.85), calc_stats("Naive", "AF_Sex", 1.10), calc_stats("Naive", "AF_READ", 0.90), calc_stats("Naive", "AF_COADREAD", 0.95), calc_stats("Naive", "Med_WT", True_Med_WT), calc_stats("Naive", "S5_WT", True_S5_WT)),
#       `Standard LT AFT` = c(calc_stats("LT", "AF_X", True_AF), calc_stats("LT", "AF_Age", 0.85), calc_stats("LT", "AF_Sex", 1.10), calc_stats("LT", "AF_READ", 0.90), calc_stats("LT", "AF_COADREAD", 0.95), calc_stats("LT", "Med_WT", True_Med_WT), calc_stats("LT", "S5_WT", True_S5_WT)),
#       `Proposed (G-comp)` = c(paste0("<b>", calc_stats("Prop", "AF_X", True_AF), "</b>"), paste0("<b>", calc_stats("Prop", "AF_Age", 0.85), "</b>"), paste0("<b>", calc_stats("Prop", "AF_Sex", 1.10), "</b>"), paste0("<b>", calc_stats("Prop", "AF_READ", 0.90), "</b>"), paste0("<b>", calc_stats("Prop", "AF_COADREAD", 0.95), "</b>"), paste0("<b>", calc_stats("Prop", "Med_WT", True_Med_WT), "</b>"), paste0("<b>", calc_stats("Prop", "S5_WT", True_S5_WT), "</b>"))
#     )
#   }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
# })




# =========================================================================
# Server Logic for Simulation Study (Tamura & Ikegami Model Ver 2.0.1)
#  - UI/outputs/inputs are kept the same
#  - Main changes:
#     (1) IPTW is now "survival-function calibration" to Macro S(t) (age-stratified 1–5y),
#         not density-matching among deaths.
#     (2) T2 is modeled as a continuous function of T1 via spline basis columns (ns1..ns3),
#         avoiding median stratification while keeping your objects/outputs unchanged.
#     (3) Random numbers are fixed/reused inside the calibration objective for stability.
#     (4) Median/5y survival are computed from model survival directly (summary + uniroot).
# =========================================================================

library(dplyr)
library(flexsurv)
library(survival)
library(ggplot2)
library(splines)

# =========================================================================
# Helper: generate complex censoring patterns (unchanged)
# =========================================================================
find_censoring_param <- function(T2, target_rate, pattern) {
  if (pattern == "peak1y") {
    obj_func <- function(k) {
      scale_val <- 365.25 / ((k - 1) / k)^(1 / k)
      (sum(rweibull(length(T2), shape = k, scale = scale_val) < T2) / length(T2)) - target_rate
    }
    opt_k <- tryCatch({ uniroot(obj_func, interval = c(1.05, 10))$root }, error = function(e) { 2.0 })
    return(list(shape = opt_k, scale = 365.25 / ((opt_k - 1) / opt_k)^(1 / opt_k), is_ushape = FALSE))
  } else if (pattern == "ushape") {
    obj_func <- function(scale_param) {
      u <- runif(length(T2))
      C2 <- ifelse(u < 0.5,
                   rweibull(length(T2), shape = 0.5, scale = scale_param),
                   rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- tryCatch({ uniroot(obj_func, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
    return(list(shape = NA, scale = opt_scale, is_ushape = TRUE))
  } else {
    obj_func_scale <- function(scale_param) {
      C2 <- if (pattern == "indep") rexp(length(T2), rate = 1 / scale_param) else rweibull(length(T2), shape = 0.5, scale = scale_param)
      (sum(C2 < T2) / length(T2)) - target_rate
    }
    opt_scale <- tryCatch({ uniroot(obj_func_scale, interval = c(0.1, 1000000))$root }, error = function(e) { median(T2) * 1.5 })
    return(list(shape = if(pattern == "early") 0.5 else 1.0, scale = opt_scale, is_ushape = FALSE))
  }
}

# =========================================================================
# GLOBAL HELPERS: Deterministic inverse-CDF sampling with supplied uniforms
#   - Works with flexsurvreg parameterization (AFT: scale multiplied by exp(eta))
#   - We keep your earlier coefficient handling but allow extra numeric columns (ns1..ns3)
# =========================================================================

.build_af_vec <- function(res_m, newdata) {
  af_vec <- rep(1, nrow(newdata))

  # Known factors
  if("XMutated" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["XMutated", "est"] * as.numeric(newdata$X == "Mutated"))
  if("Age" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["Age", "est"] * newdata$Age)
  if("SexFemale" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["SexFemale", "est"] * as.numeric(newdata$Sex == "Female"))
  if("HistologyREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyREAD", "est"] * as.numeric(newdata$Histology == "READ"))
  if("HistologyCOADREAD" %in% rownames(res_m)) af_vec <- af_vec * exp(res_m["HistologyCOADREAD", "est"] * as.numeric(newdata$Histology == "COADREAD"))

  # Extra numeric columns (e.g., ns1..ns3). Must exist in newdata with same names as coefficients.
  extra_names <- setdiff(rownames(res_m), c("scale", "shape", "XMutated", "Age", "SexFemale", "HistologyREAD", "HistologyCOADREAD"))
  for(nm in extra_names) {
    if(nm %in% names(newdata)) {
      af_vec <- af_vec * exp(res_m[nm, "est"] * as.numeric(newdata[[nm]]))
    }
  }
  af_vec
}

gen_sim_times_weibull_u <- function(fit, newdata, u) {
  if(is.null(fit) || nrow(newdata) == 0) return(rep(NA, nrow(newdata)))
  res_m <- fit$res
  base_scale <- res_m["scale", "est"]
  shape <- res_m["shape", "est"]

  af_vec <- .build_af_vec(res_m, newdata)
  scale_vec <- base_scale * af_vec

  u <- pmax(pmin(u, 0.9999), 0.0001)
  scale_vec * (-log(u))^(1 / shape)
}

gen_sim_times_llogis_u <- function(fit, newdata, u) {
  if(is.null(fit) || nrow(newdata) == 0) return(rep(NA, nrow(newdata)))
  res_m <- fit$res
  base_scale <- res_m["scale", "est"]
  shape <- res_m["shape", "est"]

  af_vec <- .build_af_vec(res_m, newdata)
  scale_vec <- base_scale * af_vec

  u <- pmax(pmin(u, 0.9999), 0.0001)
  scale_vec * ((1 - u) / u)^(1 / shape)
}

# =========================================================================
# NEW: Survival-function calibration weights to Macro S(t) at 1..5 years (age-stratified)
#   - We calibrate weights as a smooth function of T1 within each age class
#   - Objective uses deterministic pseudo-sampling with fixed uniforms
#   - Then refit T1/T2 models using calibrated weights
# =========================================================================
.calib_weights_from_theta <- function(data, theta) {
  # theta = c(Young0, Young1, Young2, Old0, Old1, Old2)
  logt1 <- log(pmax(data$T1 / 365.25, 1e-6))
  logt1_sq <- logt1^2

  is_y <- as.numeric(data$Age_class == "Young")
  is_o <- 1 - is_y

  eta <- is_y * (theta[1] + theta[2] * logt1 + theta[3] * logt1_sq) +
    is_o * (theta[4] + theta[5] * logt1 + theta[6] * logt1_sq)

  w <- exp(eta)

  # Trim to stabilize
  lo <- quantile(w, 0.01, na.rm = TRUE)
  hi <- quantile(w, 0.99, na.rm = TRUE)
  w <- pmax(lo, pmin(w, hi))

  # Normalize mean to 1
  w / mean(w, na.rm = TRUE)
}

# Deterministic resampling given probabilities and fixed uniforms
.resample_idx_from_probs <- function(probs, u_vec) {
  p <- probs / sum(probs)
  cdf <- cumsum(p)
  findInterval(u_vec, cdf) + 1L
}

# Compute age-stratified S_hat(t) from simulated OS (in years) in a pseudo sample
.calc_Shat_from_simOS <- function(sim_OS_days, age_class, t_years = 1:5) {
  out <- list()
  for(ag in c("Young", "Old")) {
    idx <- which(age_class == ag)
    if(length(idx) < 10) {
      out[[ag]] <- rep(NA_real_, length(t_years))
    } else {
      os_y <- sim_OS_days[idx]
      out[[ag]] <- sapply(t_years, function(ty) mean(os_y > ty * 365.25))
    }
  }
  out
}

# Calibrate weights so that g-comp generated S_hat(g,t) matches Macro S_macro(g,t)
calibrate_weights_surv <- function(Data_cgp, fit_t1_init, fit_t2_init, ns_obj,
                                   S_macro_list, # list(Young=vector length5 in [0,1], Old=...)
                                   M_pseudo = 8000,
                                   seed = 12345) {

  if(is.null(fit_t1_init) || is.null(fit_t2_init)) {
    Data_cgp$iptw <- 1.0
    return(Data_cgp)
  }

  set.seed(seed)
  u_idx <- runif(M_pseudo)      # for deterministic resampling
  u_t1  <- runif(M_pseudo)      # for T1 inverse sampling
  u_t2  <- runif(M_pseudo)      # for T2 inverse sampling

  # Objective: squared error of S_hat vs S_macro at t=1..5, both age groups, plus weight-variance penalty
  obj <- function(theta) {
    w <- .calib_weights_from_theta(Data_cgp, theta)
    idx <- .resample_idx_from_probs(w, u_idx)
    pseudo <- Data_cgp[idx, , drop = FALSE]

    # simulate T1
    pseudo$sim_T1 <- gen_sim_times_weibull_u(fit_t1_init, pseudo, u_t1)

    # build spline basis for T2 using simulated T1 (continuous dependence)
    pseudo$logT1_sim <- log(pmax(pseudo$sim_T1, 1e-6))
    ns_mat <- tryCatch({ predict(ns_obj, newx = pseudo$logT1_sim) },
                       error = function(e) matrix(0, nrow(pseudo), attr(ns_obj, "df")))
    colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
    for(j in seq_len(ncol(ns_mat))) pseudo[[colnames(ns_mat)[j]]] <- ns_mat[, j]

    # simulate T2
    pseudo$sim_T2 <- gen_sim_times_llogis_u(fit_t2_init, pseudo, u_t2)
    pseudo$sim_OS <- pseudo$sim_T1 + pseudo$sim_T2

    Shat <- .calc_Shat_from_simOS(pseudo$sim_OS, pseudo$Age_class, t_years = 1:5)

    # SSE over available values
    sse <- 0
    for(ag in c("Young", "Old")) {
      target <- S_macro_list[[ag]]
      est <- Shat[[ag]]
      ok <- is.finite(target) & is.finite(est)
      if(any(ok)) sse <- sse + sum((est[ok] - target[ok])^2)
    }

    # penalty to avoid extreme weights
    pen <- 0.01 * var(log(pmax(w, 1e-12)))
    sse + pen
  }

  theta0 <- rep(0, 6)
  opt <- tryCatch({
    optim(theta0, obj, method = "BFGS", control = list(maxit = 80))
  }, error = function(e) NULL)

  theta_hat <- if(is.null(opt)) theta0 else opt$par
  Data_cgp$iptw <- .calib_weights_from_theta(Data_cgp, theta_hat)
  Data_cgp
}

# =========================================================================
# Core Simulation Engine (UPDATED)
#  - Keeps return objects / names used by UI
#  - Naive and LT comparators unchanged
#  - Proposed:
#     * Initial T1/T2 fit (unweighted)
#     * Calibrate weights to Macro S(t) via survival-function matching
#     * Refit T1/T2 with weights
#     * G-comp to generate pseudo, then fit final OS model to pseudo for synced table/plot
# =========================================================================
run_sim_iteration <- function(N_target, True_AF_X, Mut_Freq, t1_pat, cens_pat, cens_rate, return_data = FALSE) {

  # --- Data Generation Process ---
  N_macro <- 100000
  Age <- round(runif(N_macro, 40, 80))
  Sex <- sample(c("Male", "Female"), N_macro, replace = TRUE, prob = c(0.55, 0.45))
  Histology <- sample(c("COAD", "READ", "COADREAD"), N_macro, replace = TRUE, prob = c(0.60, 0.30, 0.10))
  X <- rbinom(N_macro, 1, Mut_Freq)

  AF_bg <- 1.0 * (0.85 ^ ((Age - 60) / 10)) *
    ifelse(Sex == "Female", 1.10, 1.0) *
    ifelse(Histology == "READ", 0.90, ifelse(Histology == "COADREAD", 0.95, 1.0))

  # True OS (days): log-logistic-like via inverse transform (median ~ 2y, shape=1.5)
  u1 <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  u2 <- pmax(pmin(runif(N_macro), 0.999), 0.001)
  T_true <- (2.0 * 365.25) * ((1 - u1) / u2)^(1 / 1.5) * AF_bg * ifelse(X == 1, True_AF_X, 1.0)

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

  params <- find_censoring_param(Data_cgp$T2_true, cens_rate, cens_pat)
  if (isTRUE(params$is_ushape)) {
    u <- runif(N_target)
    C2 <- ifelse(u < 0.5,
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
  if(sum(Data_cgp$Event) < 20) return(NULL)

  # Macro reference S(t): age-stratified 1..5 years (truth in simulation; real use would be registry)
  ref_surv_list <- list()
  for (ag in c("Young", "Old")) {
    ref_surv_list[[ag]] <- sapply(1:5, function(y) mean(T_true[ifelse(Age < 60, "Young", "Old") == ag] > y * 365.25) * 100)
  }
  S_macro_list <- list(
    Young = pmax(pmin(ref_surv_list$Young[1:5] / 100, 0.999), 0.001),
    Old   = pmax(pmin(ref_surv_list$Old[1:5]   / 100, 0.999), 0.001)
  )

  # --- Covariates used in models (kept similar to your original logic) ---
  valid_covs <- c("X", "Age")
  if (length(unique(Data_cgp$Sex)) > 1) valid_covs <- c(valid_covs, "Sex")
  if (length(unique(Data_cgp$Histology)) > 1) valid_covs <- c(valid_covs, "Histology")

  # --- Standard Comparator Models (unchanged) ---
  form_naive <- as.formula(paste("Surv(T_obs, Event) ~", paste(valid_covs, collapse=" + ")))
  form_lt    <- as.formula(paste("Surv(T1, T_obs, Event) ~", paste(valid_covs, collapse=" + ")))

  fit_naive <- tryCatch({ flexsurvreg(form_naive, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)
  fit_lt    <- tryCatch({ flexsurvreg(form_lt,    data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # Proposed method (UPDATED):
  #   Step P1: Build spline basis columns for log(T1) in observed CGP data
  # =========================================================================
  Data_cgp$logT1 <- log(pmax(Data_cgp$T1, 1e-6))
  ns_obj <- ns(Data_cgp$logT1, df = 3)
  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for(j in seq_len(ncol(ns_mat))) Data_cgp[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  # =========================================================================
  # Step P2: Initial (unweighted) T1 model (Weibull) and T2 model (llogis + spline(T1))
  # =========================================================================
  Data_cgp$T1_event <- 1
  form_t1 <- as.formula(paste("Surv(T1, T1_event) ~", paste(valid_covs, collapse=" + ")))
  fit_t1_init <- tryCatch({ flexsurvreg(form_t1, data = Data_cgp, dist = "weibull") }, error = function(e) NULL)

  form_t2 <- as.formula(paste("Surv(T2_obs, Event) ~", paste(c(valid_covs, colnames(ns_mat)), collapse=" + ")))
  fit_t2_init <- tryCatch({ flexsurvreg(form_t2, data = Data_cgp, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # Step P3: Calibrate weights by matching Macro S(t) (age-stratified 1..5y)
  # =========================================================================
  Data_cgp <- calibrate_weights_surv(
    Data_cgp = Data_cgp,
    fit_t1_init = fit_t1_init,
    fit_t2_init = fit_t2_init,
    ns_obj = ns_obj,
    S_macro_list = S_macro_list,
    M_pseudo = 8000,
    seed = 20240227
  )

  # =========================================================================
  # Step P4: Refit T1/T2 using calibrated weights (iptw)
  # =========================================================================
  fit_t1 <- tryCatch({ flexsurvreg(form_t1, data = Data_cgp, weights = iptw, dist = "weibull") }, error = function(e) NULL)
  fit_t2 <- tryCatch({ flexsurvreg(form_t2, data = Data_cgp, weights = iptw, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # Step P5: G-computation to generate complete OS and fit final OS model
  # =========================================================================
  M_final <- 10000
  set.seed(20240228)
  u_idx_final <- runif(M_final)
  idx_pseudo <- .resample_idx_from_probs(Data_cgp$iptw, u_idx_final)
  pseudo_cgp <- Data_cgp[idx_pseudo, , drop = FALSE]

  # deterministic uniforms for times
  u_t1_final <- runif(M_final)
  u_t2_final <- runif(M_final)

  pseudo_cgp$sim_T1 <- gen_sim_times_weibull_u(fit_t1, pseudo_cgp, u_t1_final)

  pseudo_cgp$logT1_sim <- log(pmax(pseudo_cgp$sim_T1, 1e-6))
  ns_sim <- tryCatch({ predict(ns_obj, newx = pseudo_cgp$logT1_sim) },
                     error = function(e) matrix(0, nrow(pseudo_cgp), attr(ns_obj, "df")))
  colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
  for(j in seq_len(ncol(ns_sim))) pseudo_cgp[[colnames(ns_sim)[j]]] <- ns_sim[, j]

  pseudo_cgp$sim_T2 <- gen_sim_times_llogis_u(fit_t2, pseudo_cgp, u_t2_final)
  pseudo_cgp$sim_OS <- pseudo_cgp$sim_T1 + pseudo_cgp$sim_T2
  pseudo_cgp$sim_Event <- 1

  form_prop_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(valid_covs, collapse=" + ")))
  fit_prop_final <- tryCatch({ flexsurvreg(form_prop_final, data = pseudo_cgp, dist = "llogis") }, error = function(e) NULL)

  # =========================================================================
  # Metrics extraction (UPDATED median/5y via model survival directly)
  # =========================================================================
  .get_coef <- function(res_m, nm) if (nm %in% rownames(res_m)) res_m[nm, "est"] else NA_real_

  .surv_pred <- function(fit, newdata, t_days) {
    if(is.null(fit)) return(NA_real_)
    s <- tryCatch({
      ss <- summary(fit, newdata = newdata, t = t_days, tidy = TRUE)
      # tidy TRUE returns data.frame with est, time
      # when multiple rows, take mean (should be one row here)
      as.numeric(ss$est[1])
    }, error = function(e) NA_real_)
    s
  }

  .median_pred <- function(fit, newdata, tmax_days = 365.25 * 20) {
    if(is.null(fit)) return(NA_real_)
    f <- function(t) .surv_pred(fit, newdata, t) - 0.5
    # bracket: start near 1 day to 20 years
    out <- tryCatch({
      lo <- 1
      hi <- tmax_days
      # ensure sign change; if not, return NA
      flo <- f(lo); fhi <- f(hi)
      if(!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
      uniroot(f, interval = c(lo, hi))$root
    }, error = function(e) NA_real_)
    out
  }

  extract_metrics <- function(fit) {
    blank_res <- c(AF_X=NA, AF_Age=NA, AF_Sex=NA, AF_READ=NA, AF_COADREAD=NA, Med_WT=NA, S5_WT=NA)
    if (is.null(fit)) return(blank_res)
    res_m <- fit$res

    AF_X <- exp(.get_coef(res_m, "XMutated"))
    AF_Age <- exp(.get_coef(res_m, "Age") * 10)
    AF_Sex <- exp(.get_coef(res_m, "SexFemale"))
    AF_READ <- exp(.get_coef(res_m, "HistologyREAD"))
    AF_COADREAD <- exp(.get_coef(res_m, "HistologyCOADREAD"))

    # Baseline: WT, Male, COAD, Age=60
    nd <- data.frame(
      Age = 60,
      X = factor("WT", levels = c("WT", "Mutated")),
      Sex = factor("Male", levels = c("Male", "Female")),
      Histology = factor("COAD", levels = c("COAD", "READ", "COADREAD"))
    )

    med_days <- .median_pred(fit, nd)
    s5 <- .surv_pred(fit, nd, 5 * 365.25)

    med_years <- ifelse(is.finite(med_days), med_days / 365.25, NA_real_)
    s5_pct <- ifelse(is.finite(s5), s5 * 100, NA_real_)

    c(AF_X=AF_X, AF_Age=AF_Age, AF_Sex=AF_Sex, AF_READ=AF_READ, AF_COADREAD=AF_COADREAD,
      Med_WT=med_years, S5_WT=s5_pct)
  }

  out <- list(
    Naive = extract_metrics(fit_naive),
    LT    = extract_metrics(fit_lt),
    Prop  = extract_metrics(fit_prop_final),
    fit_naive = fit_naive,
    fit_lt = fit_lt,
    fit_prop_final = fit_prop_final,
    pseudo_cgp = pseudo_cgp,
    T_true_macro = T_true
  )

  if(return_data) return(out) else return(list(Naive=out$Naive, LT=out$LT, Prop=out$Prop))
}

safe_fmt <- function(x, digits=2) {
  sapply(x, function(v) if(is.na(v) || !is.numeric(v)) "N/A" else sprintf(paste0("%.", digits, "f"), v))
}

# -------------------------------------------------------------------------
# UI Observers and Rendering Logic (UNCHANGED VARIABLE NAMES)
# -------------------------------------------------------------------------
observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  showNotification("Running simulation (Tamura & Ikegami Ver 2.0.1)...", type = "message", duration = 3)

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
  True_Med_WT <- 2.0
  True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100

  output$sim_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Baseline Med OS [Years]", "Baseline 5-Year OS [%]"),
      `True` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
      `Naive AFT (Surv(T_obs))` = c(safe_fmt(res$Naive["AF_X"]), safe_fmt(res$Naive["AF_Age"]), safe_fmt(res$Naive["AF_Sex"]),
                                    safe_fmt(res$Naive["AF_READ"]), safe_fmt(res$Naive["AF_COADREAD"]),
                                    safe_fmt(res$Naive["Med_WT"]), safe_fmt(res$Naive["S5_WT"], 1)),
      `Standard LT (Surv(T1,T_obs))` = c(safe_fmt(res$LT["AF_X"]), safe_fmt(res$LT["AF_Age"]), safe_fmt(res$LT["AF_Sex"]),
                                         safe_fmt(res$LT["AF_READ"]), safe_fmt(res$LT["AF_COADREAD"]),
                                         safe_fmt(res$LT["Med_WT"]), safe_fmt(res$LT["S5_WT"], 1)),
      `Proposed (G-comp)` = c(paste0("<b>", safe_fmt(res$Prop["AF_X"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["AF_Age"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["AF_Sex"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["AF_READ"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["AF_COADREAD"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["Med_WT"]), "</b>"),
                              paste0("<b>", safe_fmt(res$Prop["S5_WT"], 1), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)

  output$sim_survival_plot <- renderPlot({
    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    km_true <- survfit(Surv(res$T_true_macro, rep(1, length(res$T_true_macro))) ~ 1)
    surv_true <- approx(km_true$time, km_true$surv, xout = t_seq, method = "constant", f = 0, rule = 2)$y

    calc_marginal <- function(fit, dataset) {
      if(is.null(fit)) return(rep(NA, length(t_seq)))
      samp <- dataset[sample(1:nrow(dataset), min(1000, nrow(dataset))), ]
      summ <- summary(fit, newdata = samp, t = t_seq, tidy = TRUE)
      aggregate(est ~ time, data = summ, FUN = mean)$est
    }

    surv_naive <- calc_marginal(res$fit_naive, res$pseudo_cgp)
    surv_lt    <- calc_marginal(res$fit_lt,    res$pseudo_cgp)
    surv_prop  <- calc_marginal(res$fit_prop_final, res$pseudo_cgp)

    plot_df <- data.frame(
      Time = rep(t_seq / 365.25, 4),
      Survival = c(surv_true, surv_naive, surv_lt, surv_prop),
      Model = factor(rep(c("1. True Marginal (Registry Macro)",
                           "2. Naive AFT (Immortal Time Bias)",
                           "3. Standard LT AFT (Dep. Truncation Bias)",
                           "4. Proposed Method (Ver 2.0.1: Calibrated G-computation)"),
                         each = length(t_seq)))
    )

    ggplot(plot_df, aes(x = Time, y = Survival, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "#e74c3c", "#f39c12", "#27ae60")) +
      scale_linetype_manual(values = c("solid", "dotted", "dashed", "solid")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(title = "Tamura & Ikegami Model Ver 2.0.1 Validation",
           subtitle = "Survival-function calibration to Macro S(t) + continuous T1->T2 dependence (spline).",
           x = "Time from Diagnosis (Years)", y = "Overall Survival Probability") +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.direction = "vertical",
            legend.title = element_blank())
  })
})

observeEvent(input$run_sim_multi, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)
  n_sims <- 400
  results <- list()

  withProgress(message = 'Running 400 Simulations...', value = 0, {
    for (i in 1:n_sims) {
      res <- run_sim_iteration(
        input$sim_n,
        input$sim_true_af,
        input$sim_mut_freq / 100,
        input$sim_t1_pattern,
        input$sim_cens_pattern,
        input$sim_cens_rate / 100
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
      incProgress(1/n_sims, detail = paste("Iteration", i, "of", n_sims))
    }
  })
  shiny::validate(shiny::need(length(results) > 0, "All simulations failed."))

  True_AF <- input$sim_true_af
  True_Med_WT <- 2.0
  True_S5_WT <- 1 / (1 + (5 / True_Med_WT)^1.5) * 100

  calc_stats <- function(model_name, metric_name, true_val) {
    vals <- sapply(results, function(x) x[[model_name]][metric_name])
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) return("N/A")
    mean_val <- mean(vals)
    mse_val <- mean((vals - true_val)^2)
    sprintf("%.2f (MSE: %.3f)", mean_val, mse_val)
  }

  output$sim_multi_result_table <- renderTable({
    data.frame(
      Metric = c("Target Gene AF", "Age (+10 yrs) AF", "Sex (Female) AF", "Histology (READ) AF", "Histology (COADREAD) AF",
                 "Baseline Med OS", "Baseline 5-Year OS [%]"),
      `True Value` = c(safe_fmt(True_AF), "0.85", "1.10", "0.90", "0.95", safe_fmt(True_Med_WT), safe_fmt(True_S5_WT, 1)),
      `Naive AFT` = c(calc_stats("Naive", "AF_X", True_AF),
                      calc_stats("Naive", "AF_Age", 0.85),
                      calc_stats("Naive", "AF_Sex", 1.10),
                      calc_stats("Naive", "AF_READ", 0.90),
                      calc_stats("Naive", "AF_COADREAD", 0.95),
                      calc_stats("Naive", "Med_WT", True_Med_WT),
                      calc_stats("Naive", "S5_WT", True_S5_WT)),
      `Standard LT AFT` = c(calc_stats("LT", "AF_X", True_AF),
                            calc_stats("LT", "AF_Age", 0.85),
                            calc_stats("LT", "AF_Sex", 1.10),
                            calc_stats("LT", "AF_READ", 0.90),
                            calc_stats("LT", "AF_COADREAD", 0.95),
                            calc_stats("LT", "Med_WT", True_Med_WT),
                            calc_stats("LT", "S5_WT", True_S5_WT)),
      `Proposed (G-comp)` = c(paste0("<b>", calc_stats("Prop", "AF_X", True_AF), "</b>"),
                              paste0("<b>", calc_stats("Prop", "AF_Age", 0.85), "</b>"),
                              paste0("<b>", calc_stats("Prop", "AF_Sex", 1.10), "</b>"),
                              paste0("<b>", calc_stats("Prop", "AF_READ", 0.90), "</b>"),
                              paste0("<b>", calc_stats("Prop", "AF_COADREAD", 0.95), "</b>"),
                              paste0("<b>", calc_stats("Prop", "Med_WT", True_Med_WT), "</b>"),
                              paste0("<b>", calc_stats("Prop", "S5_WT", True_S5_WT), "</b>"))
    )
  }, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
})
