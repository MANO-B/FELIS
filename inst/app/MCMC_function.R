# メイン関数：指定した2セットのIDで生存曲線を描画（reactiveValues対応）
run_custom_survival_analysis <- function(Data_survival, id_set_1, id_set_2, label_1 = "Group 1", label_2 = "Group 2") {
  tryCatch({
    # データフレームに変換
    Data_survival <- as.data.frame(Data_survival)
    # データ作成
    Data_survival_1 <- Data_survival %>%
      dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% id_set_1) %>%
      dplyr::mutate(group = label_1)
    Data_survival_2 <- Data_survival %>%
      dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% id_set_2) %>%
      dplyr::mutate(group = label_2)
    Data_survival <- rbind(Data_survival_1, Data_survival_2)
    # 従来の生存解析
    traditional_fit <- survfit(Surv(event = censor, time = time_palliative_final) ~ group, data = Data_survival)
    # 2群存在チェック
    if (length(unique(Data_survival$group)) != 2) {
      warning("2群が存在しません")
      return(NULL)
    }
    # 群名の設定
    groups <- sort(unique(Data_survival$group))
    name_1 <- groups[1]
    name_2 <- groups[2]
    name_1_short <- name_1
    name_2_short <- name_2
    # group_binary列の追加（ベイズ解析用）
    Data_survival$group <- Data_survival$group == groups[2]
    Data_survival_curve <- Data_survival
    # 時間窓処理とデータ準備（元のコードの複雑な処理を関数化）
    processed_data <- process_time_windows_original(Data_survival_curve)
    if (is.null(processed_data)) {
      warning("データ処理に失敗しました")
      return(NULL)
    }
    Data_drug_survival_1 <- processed_data$survival_data
    Data_drug_survival_1_2 <- processed_data$delay_data
    Median_entry_pos <- processed_data$median_pos
    Median_entry_neg <- processed_data$median_neg
    
    # delay_1カテゴリの追加
    Data_drug_survival_1 <- Data_drug_survival_1 %>%
      dplyr::mutate(delay_1 = case_when(
        time_palliative_enroll <= Median_entry_pos & group == TRUE ~ "SPT to entry early",
        time_palliative_enroll > Median_entry_pos & group == TRUE ~ "SPT to entry later",
        time_palliative_enroll <= Median_entry_neg & group == FALSE ~ "SPT to entry early",
        time_palliative_enroll > Median_entry_neg & group == FALSE ~ "SPT to entry later"
      ))
    
    # ベイズ解析実行条件チェック
    bayesian_conditions <- check_bayesian_conditions_original(Data_drug_survival_1)
    
    if (bayesian_conditions) {
      # Stanデータの準備
      stan_data <- prepare_stan_data_original(Data_drug_survival_1, Data_drug_survival_1_2, Median_entry_pos, Median_entry_neg)
      
      # Stanモデルの実行
      stan_fit <- run_stan_model_original(stan_data)
      
      # 結果の処理（元のコードと同じ処理）
      results <- process_stan_results_original(stan_fit, Data_survival_curve, traditional_fit, name_1, name_2, name_1_short, name_2_short)
      OUTPUT_DATA$figure_surv_Bayes_compare <- results$plot
    } else {
      warning("ベイズ解析の実行条件が満たされていません")
      return(NULL)
    }
    
  }, error = function(e) {
    warning(paste("遺伝子変異生存解析でエラーが発生しました:", e$message))
    return(NULL)
  })
}

# 時間窓処理とデータ準備（元のコードの複雑な処理を再現）
process_time_windows_original <- function(Data_survival_curve) {
  
  tryCatch({
    # 変異ありグループの処理
    Data_survival_tmp_pos <- Data_survival_curve %>% dplyr::filter(group == TRUE)
    processed_pos <- process_single_time_window(Data_survival_tmp_pos, "time_enroll_final", "censor")
    
    Data_survival_tmp_pos_2 <- Data_survival_curve %>% dplyr::filter(group == TRUE)
    processed_pos_2 <- process_single_time_window(Data_survival_tmp_pos_2, "time_palliative_enroll", "enroll_censor")
    
    # 変異なしグループの処理
    Data_survival_tmp_neg <- Data_survival_curve %>% dplyr::filter(group == FALSE)
    processed_neg <- process_single_time_window(Data_survival_tmp_neg, "time_enroll_final", "censor")
    
    Data_survival_tmp_neg_2 <- Data_survival_curve %>% dplyr::filter(group == FALSE)
    processed_neg_2 <- process_single_time_window(Data_survival_tmp_neg_2, "time_palliative_enroll", "enroll_censor")
    
    # データの結合
    Data_drug_survival_1 <- rbind(processed_pos$adjusted_data, processed_neg$adjusted_data)
    Data_drug_survival_1_2 <- rbind(processed_pos_2$adjusted_data, processed_neg_2$adjusted_data)
    
    # 中央値の計算
    Median_entry_pos <- median(processed_pos_2$adjusted_data$time_palliative_enroll[processed_pos_2$adjusted_data$group == TRUE], na.rm = TRUE)
    Median_entry_neg <- median(processed_neg_2$adjusted_data$time_palliative_enroll[processed_neg_2$adjusted_data$group == FALSE], na.rm = TRUE)
    
    return(list(
      survival_data = Data_drug_survival_1,
      delay_data = Data_drug_survival_1_2,
      median_pos = Median_entry_pos,
      median_neg = Median_entry_neg
    ))
    
  }, error = function(e) {
    warning(paste("時間窓処理でエラー:", e$message))
    return(NULL)
  })
}

# 単一時間窓の処理（元のコードの処理を関数化）
process_single_time_window <- function(data, time_col, event_col) {
  
  if (nrow(data) == 0) {
    return(list(adjusted_data = data.frame()))
  }
  
  if (sum(data[[event_col]], na.rm = TRUE) > 0) {
    fit <- survfit(as.formula(paste("Surv(", time_col, ",", event_col, ") ~ 1")), data = data)
    time_10 <- max(90, quantile(fit, 0.10)$quantile[1], na.rm = TRUE)
    survival_rate_10 <- summary(fit, times = time_10)[[6]]
    survival_rate_90 <- summary(fit, times = max(data[[time_col]], na.rm = TRUE) * 0.9)[[6]]
  } else {
    time_10 <- 90
    survival_rate_10 <- 1
    survival_rate_90 <- 0
  }
  
  time_90 <- max(time_10 + 1, max(data[[time_col]], na.rm = TRUE) * 0.9)
  
  # 中間期間のデータ
  data_middle <- data %>% dplyr::filter(
    .data[[time_col]] >= time_10 & .data[[time_col]] <= time_90
  )
  data_middle[[time_col]] <- data_middle[[time_col]] - (time_10 - 1)
  
  # 後期データの処理
  data_late <- data %>% dplyr::filter(.data[[time_col]] > time_90) %>% dplyr::arrange(.data[[time_col]])
  data_late[[time_col]] <- time_90 - (time_10 - 1)
  
  if (nrow(data_late) > 0) {
    idx <- 1:nrow(data_late)
    sample_size <- floor(max(idx) * (1 - survival_rate_10))
    if (sample_size > 0 && sample_size < nrow(data_late)) {
      data_late <- data_late[idx[idx < sample_size], ]
    }
  }
  
  adjusted_data <- rbind(data_middle, data_late)
  
  return(list(adjusted_data = adjusted_data))
}

# ベイズ解析実行条件チェック（元のコードと同じ条件）
check_bayesian_conditions_original <- function(Data_drug_survival_1) {
  
  if (!"delay_1" %in% colnames(Data_drug_survival_1)) {
    return(FALSE)
  }
  
  conditions <- c(
    sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE) >= 1,
    sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE) >= 1,
    sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE) >= 1,
    sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE) >= 1
  )
  
  return(all(conditions))
}


# Stanデータの準備（元のコードと同じ構造）
prepare_stan_data_original <- function(Data_drug_survival_1, Data_drug_survival_1_2, Median_entry_pos, Median_entry_neg) {
  
  stan_data <- list(
    Median_pos = as.integer(Median_entry_pos),
    Median_neg = as.integer(Median_entry_neg),
    Nobs_early_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE),
    Nobs_late_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE),
    Ncen_early_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE),
    Ncen_late_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group != FALSE, na.rm = TRUE),
    Nobs_early_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE),
    Nobs_late_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE),
    Ncen_early_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE),
    Ncen_late_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group == FALSE, na.rm = TRUE),
    Nexp_pos = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$group != FALSE]),
    Nexp_neg = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$group == FALSE]),
    ybef_exp_pos = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$group != FALSE],
    ybef_exp_neg = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$group == FALSE],
    yobs_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group != FALSE],
    yobs_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group != FALSE],
    ycen_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group != FALSE],
    ycen_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group != FALSE],
    yobs_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group == FALSE],
    yobs_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group == FALSE],
    ycen_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$group == FALSE],
    ycen_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$group == FALSE]
  )
  
  return(stan_data)
}

# Stanモデルの実行（元のコードと同じ）
run_stan_model_original <- function(stan_data) {
  if (DOCKER) {
    stan_model_compiled <- readRDS(stan_model_factor_save_path_cmdstan)
    stan_fit <- stan_model_compiled$sample(
      data = stan_data,
      seed = ANALYSIS_CONFIG$seed,
      chains = ANALYSIS_CONFIG$chains,
      parallel_chains = ANALYSIS_CONFIG$parallel_chains,
      iter_sampling = ANALYSIS_CONFIG$iter_sampling,
      iter_warmup = ANALYSIS_CONFIG$iter_warmup,
      adapt_delta = ANALYSIS_CONFIG$adapt_delta,
      max_treedepth = ANALYSIS_CONFIG$max_treedepth,
      thin = ANALYSIS_CONFIG$thin,
      refresh = ANALYSIS_CONFIG$refresh
    )
  } else {
    mod <- readRDS(stan_model_factor_save_path_rstan)
    stan_fit <- rstan::sampling(object = mod,
      data = stan_data,
      control = list(adapt_delta = ANALYSIS_CONFIG$adapt_delta, max_treedepth = ANALYSIS_CONFIG$max_treedepth),
      seed = ANALYSIS_CONFIG$seed,
      chains = ANALYSIS_CONFIG$chains,
      iter = as.integer(ANALYSIS_CONFIG$iter_sampling * 3 / 2),
      warmup = ANALYSIS_CONFIG$iter_warmup,
      thin = ANALYSIS_CONFIG$thin,
      refresh = ANALYSIS_CONFIG$refresh,
      verbose = FALSE
    )
  }
  
  return(stan_fit)
}

# Stan結果の処理（元のコードと同じ処理）
process_stan_results_original <- function(stan_fit, Data_survival_curve, traditional_fit, name_1, name_2, name_1_short, name_2_short) {
  
  # MCMCサンプルの抽出（元のコードと同じ）
  stan_draws <- tidybayes::tidy_draws(stan_fit)
  
  stan_draws_total <-
    stan_draws %>%
    dplyr::select(.chain, .iteration, .draw, starts_with("yhat_total")) %>%
    tidyr::gather(key = key, value = yhat_total, starts_with("yhat_total")) %>%
    dplyr::select(-key) 
  stan_draws_total_exp <-
    stan_draws %>%
    dplyr::select(.chain, .iteration, .draw, starts_with("yhat_factor_total")) %>%
    tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_factor_total")) %>%
    dplyr::select(-key) 

  # 中央値の計算（元のコードと同じ）
  median_os_list_weibull_total_exp <- matrix(stan_draws_total_exp$yhat_total_exp, nrow = 4 * ITER)
  median_os_list_weibull_total <- matrix(stan_draws_total$yhat_total, nrow = 4 * ITER)
  median_os_list_weibull_total_exp <- rowMedians(median_os_list_weibull_total_exp, na.rm = TRUE)
  median_os_list_weibull_total <- rowMedians(median_os_list_weibull_total, na.rm = TRUE)
  median_os_list_weibull_bias <- median_os_list_weibull_total - median_os_list_weibull_total_exp
  
  # 信頼区間の計算
  median_os_summary_weibull_total <- quantile(median_os_list_weibull_total, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  median_os_summary_weibull_total_exp <- quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  median_os_summary_weibull_bias <- quantile(median_os_list_weibull_bias, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  # ソート済み推定値
  yhat_weibull_sort_total <- sort(stan_draws_total$yhat_total)
  yhat_weibull_sort_total_exp <- sort(stan_draws_total_exp$yhat_total_exp)
  yhat_weibull_sort_average_total <- yhat_weibull_sort_total[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total)/length(Data_survival_curve$censor))]
  yhat_weibull_sort_average_total_exp <- yhat_weibull_sort_total_exp[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total_exp)/length(Data_survival_curve$censor))]
  
  # 調整済み時間の設定
  Data_survival_curve$factor_neg <- yhat_weibull_sort_average_total
  Data_survival_curve$factor_pos <- yhat_weibull_sort_average_total_exp
  Data_survival_curve$factor_neg[Data_survival_curve$factor_neg > 10000] <- 10000
  Data_survival_curve$factor_pos[Data_survival_curve$factor_pos > 10000] <- 10000
  
  # 生存曲線のフィッティング
  traditional_fit_gene <- survfit(Surv(event = censor, time = time_palliative_final) ~ group, data = Data_survival_curve)
  simulation_fit_neg <- survfit(Surv(event = rep(1, length(Data_survival_curve$factor_neg)), time = factor_neg) ~ 1, data = Data_survival_curve)
  simulation_fit_pos <- survfit(Surv(event = rep(1, length(Data_survival_curve$factor_pos)), time = factor_pos) ~ 1, data = Data_survival_curve)
  
  # プロットの作成（元のコードと同じ）
  plot_obj <- ggsurvplot(
    fit = list(traditional_fit = traditional_fit_gene,
               simulation_fit_neg = simulation_fit_neg,
               simulation_fit_pos = simulation_fit_pos),
    combine = TRUE,
    data = Data_survival_curve,
    xlab = "Time from CTx induction to final observation (months)",
    ylab = "Survival Probability",
    censor = TRUE,
    font.title = 8,
    font.subtitle = 8,
    font.main = 8,
    font.submain = 8,
    font.caption = 8,
    font.legend = 8,
    surv.median.line = "hv",
    conf.int = FALSE,
    pval = FALSE,
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    tables.theme = clean_theme(),
    legend = c(0.8, 0.8),
    xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
    cumevents = FALSE,
    cumcensor = FALSE,
    break.x.by = 365.25 * 2.5,
    xscale = "d_m",
    legend.labs = c(paste0(name_1, ", Unadjusted"),
                    paste0(name_2, ", Unadjusted"),
                    paste0(name_1, ", Adjusted for Delayed Entry"),
                    paste0(name_2, ", Adjusted for Delayed Entry"))
  ) + 
    labs(title = paste(sum(traditional_fit_gene$n), " patients ",
                       "Median OS ",
                       format_p(summary(traditional_fit_gene)$table[[13]] / 365.25 * 12, digits = 1), " (",
                       format_p(summary(traditional_fit_gene)$table[[15]] / 365.25 * 12, digits = 1), "-",
                       format_p(summary(traditional_fit_gene)$table[[17]] / 365.25 * 12, digits = 1), ") (",
                       name_1_short, ", unadj.)/",
                       format_p(summary(traditional_fit_gene)$table[[14]] / 365.25 * 12, digits = 1), " (",
                       format_p(summary(traditional_fit_gene)$table[[16]] / 365.25 * 12, digits = 1), "-",
                       format_p(summary(traditional_fit_gene)$table[[18]] / 365.25 * 12, digits = 1), ") (",
                       name_2_short, ", unadj.)/",
                       format_p(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1), 
                       " (", format_p(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1), "-",
                       format_p(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1), ") (",
                       name_1_short, ", adj.)/",
                       format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                       format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                       format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                       ") (", name_2_short, ", adj.) months",
                       sep=""),
         subtitle = paste("Survival difference, ", name_1_short, "-", name_2_short, ": ",
                          format_p(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 1), " (",
                          format_p(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 1), "-", 
                          format_p(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 1), 
                          ") months, median (95% CI)",
                          sep=""))
  
  plot_obj$table <- plot_obj$table + theme(plot.title = element_blank(), plot.subtitle = element_blank())
  
  return(list(
    plot = plot_obj,
    median_os_summary_weibull_total = median_os_summary_weibull_total,
    median_os_summary_weibull_total_exp = median_os_summary_weibull_total_exp,
    median_os_summary_weibull_bias = median_os_summary_weibull_bias
  ))
}
