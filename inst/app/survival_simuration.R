observeEvent(input$survival_simuration, {
  analysis_env <- new.env()
  clear_reactive_data(OUTPUT_DATA)
  gc()
  req(input$shape_llogis, input$scale_llogis, input$shape_CGP, input$scale_CGP)
  with(analysis_env, {
    withProgress(message = sample(nietzsche)[1], {
      shape_llogis = input$shape_llogis
      median_os = 365.25*input$scale_llogis
      shape_CGP = input$shape_CGP
      scale_CGP = 365.25*input$scale_CGP
      shape_followup = input$shape_weibull
      scale_followup = 365.25*input$scale_weibull
      patient_no = input$test_patients
      censor_rate = input$censor_rate
      immortal_length = 30
      
      Data_survival = as.data.frame(paste0("Patient_", 1:patient_no))
      colnames(Data_survival) = "ID"
      
      Data_survival$true_OS = qllogis(((patient_no-1):0)/patient_no, shape = shape_llogis, scale = median_os, log = FALSE)
      
      Data_survival$true_OS = ifelse(Data_survival$true_OS>=365.25*10, 365.25*10, Data_survival$true_OS)
      Data_survival$true_censor = 1
      
      Data_survival$followup_date = rweibull(patient_no, shape = shape_followup, scale = scale_followup/(log(2)^(1/shape_followup)))
      Data_survival$time_palliative_enroll = Data_survival$true_OS - rllogis(patient_no, shape = shape_CGP, scale = scale_CGP)
      Data_survival$time_enroll_final = Data_survival$true_OS - Data_survival$time_palliative_enroll
      
      Data_survival$censor_date = runif(patient_no) * Data_survival$time_enroll_final
      Data_survival$censor = as.integer( runif(patient_no, min = 0, max = 1) + censor_rate)
      Data_survival$censor = 1 - ifelse(Data_survival$censor > 1, 1, Data_survival$censor)
      Data_survival$time_enroll_final = ifelse(Data_survival$censor == 0, Data_survival$censor_date, Data_survival$time_enroll_final)
      
      Data_survival$CGP_undergo = ifelse(Data_survival$time_enroll_final > immortal_length &
                                           Data_survival$time_enroll_final < Data_survival$followup_date &
                                           Data_survival$time_palliative_enroll > 0 &
                                           Data_survival$time_palliative_enroll < Data_survival$true_OS, 1, 0)
      
      Data_survival$observed_OS = Data_survival$time_palliative_enroll + Data_survival$time_enroll_final
      
      
      Data_survival_curve = Data_survival %>% dplyr::filter(CGP_undergo == 1 & observed_OS > immortal_length)
      if(sum(Data_survival_curve$censor) > 0 && FALSE){
        fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_curve)
        time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
        survival_rate_10 = summary(fit,times = time_10)[[6]]
        survival_rate_90 = summary(fit,times = max(Data_survival_curve$time_enroll_final)*0.9)[[6]]
      } else {
        time_10 = 90
        survival_rate_10 = 1
        survival_rate_90 = 0
      }
      time_90 = max(time_10 + 1, max(Data_survival_curve$time_enroll_final)*0.9)
      Data_survival_tmp = Data_survival_curve %>% dplyr::filter(
        time_enroll_final >= time_10 &
          time_enroll_final <= time_90)
      Data_survival_tmp$time_enroll_final =
        Data_survival_tmp$time_enroll_final - (time_10 - 1)
      Data_survival_tmp2 = Data_survival_curve %>% dplyr::filter(
        time_enroll_final > time_90) %>%
        dplyr::arrange(time_enroll_final)
      Data_survival_tmp2$time_enroll_final = time_90 - (time_10 - 1)
      idx <- 1:nrow(Data_survival_tmp2)
      Data_drug_survival_1 = rbind(Data_survival_tmp, Data_survival_tmp2[idx[idx < max(idx)*(1-survival_rate_10)],])
      
      if(sum(Data_survival_curve$censor) > 0 && FALSE){
        fit = survfit(Surv(time_palliative_enroll, true_censor) ~ 1, data = Data_survival_curve)
        time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
        survival_rate_10 = summary(fit,times = time_10)[[6]]
        survival_rate_90 = summary(fit,times = max(Data_survival_curve$time_palliative_enroll)*0.9)[[6]]
      } else {
        time_10 = 90
        survival_rate_10 = 1
        survival_rate_90 = 0
      }
      time_90 = max(time_10 + 1, max(Data_survival_curve$time_palliative_enroll)*0.9)
      Data_survival_tmp5 = Data_survival_curve %>% dplyr::filter(
        time_palliative_enroll >= time_10 &
          time_palliative_enroll <= time_90)
      Data_survival_tmp5$time_palliative_enroll =
        Data_survival_tmp5$time_palliative_enroll - (time_10 - 1)
      Data_survival_tmp6 = Data_survival_curve %>% dplyr::filter(
        time_palliative_enroll > time_90) %>%
        dplyr::arrange(time_palliative_enroll)
      Data_survival_tmp6$time_palliative_enroll = time_90 - (time_10 - 1)
      idx <- 1:nrow(Data_survival_tmp6)
      Data_drug_survival_2 = rbind(Data_survival_tmp5, Data_survival_tmp6[idx[idx < max(idx)*(1-survival_rate_10)],])
      
      Median_entry = median(Data_survival_tmp5$time_palliative_enroll)
      
      Data_drug_survival_1 = Data_drug_survival_1 %>%
        dplyr::mutate(delay_1 = case_when(
          Data_drug_survival_1$time_palliative_enroll <= Median_entry ~ "SPT to entry early",
          TRUE ~ "SPT to entry later"))
      
      stan_weibull_survival_model_data <-
        list(
          Median = as.integer(Median_entry),
          Nobs_early = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later"),
          Nobs_late = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later"),
          Ncen_early = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later"),
          Ncen_late = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later"),
          M_bg = 1,
          Nexp = nrow(Data_drug_survival_2),
          ybef_exp = Data_drug_survival_2$time_palliative_enroll,
          yobs_early = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later"],
          yobs_late = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later"],
          ycen_early = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later"],
          ycen_late = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later"]
        )
      
      if(DOCKER){
        stan_model_compiled <- readRDS(stan_model_simple_save_path_cmdstan)
        stan_weibull_survival_model_fit <- stan_model_compiled$sample(
          data = stan_weibull_survival_model_data,
          seed = 1234,
          chains = 4,
          parallel_chains = PARALLEL,
          iter_sampling = 1000,  # (iter-warmup)
          iter_warmup = 500,adapt_delta=0.99, max_treedepth = 10,
          thin = 1,
          refresh = 500
        )
      } else {
        mod <- readRDS(stan_model_simple_save_path_rstan)
        stan_weibull_survival_model_fit <-
          rstan::sampling(object = mod,
                      data = stan_weibull_survival_model_data,
                      control=list(adapt_delta=0.99, max_treedepth = 10),
                      seed=1234, chains=4, iter=1500, warmup=500, thin=1, refresh=500, verbose = FALSE)
      }
      incProgress(1 / 3)
      
      stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_weibull_survival_model_fit)
      stan_weibull_survival_model_draws_total_exp <-
        stan_weibull_survival_model_draws %>%
        dplyr::select(.chain, .iteration, .draw, starts_with("yhat_exp_total")) %>%
        tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_exp_total")) %>%
        dplyr::select(-key) 
      stan_weibull_survival_model_draws_bef_exp <-
        stan_weibull_survival_model_draws %>%
        dplyr::select(.chain, .iteration, .draw, starts_with("y_exptmp")) %>%
        tidyr::gather(key = key, value = yhat_bef_exp, starts_with("y_exptmp")) %>%
        dplyr::select(-key) 
      stan_weibull_survival_model_draws_aft_exp <-
        stan_weibull_survival_model_draws %>%
        dplyr::select(.chain, .iteration, .draw, starts_with("yhat_exp_uncens")) %>%
        tidyr::gather(key = key, value = yhat_aft_exp, starts_with("yhat_exp_uncens")) %>%
        dplyr::select(-key) 
      stan_weibull_survival_model_draws_ear <-
        stan_weibull_survival_model_draws %>%
        dplyr::select(.chain, .iteration, .draw, starts_with("yhat_early")) %>%
        tidyr::gather(key = key, value = yhat_ear, starts_with("yhat_early")) %>%
        dplyr::select(-key) 
      stan_weibull_survival_model_draws_lat <-
        stan_weibull_survival_model_draws %>%
        dplyr::select(.chain, .iteration, .draw, starts_with("yhat_late")) %>%
        tidyr::gather(key = key, value = yhat_lat, starts_with("yhat_late")) %>%
        dplyr::select(-key) 
      
      median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=8000)
      median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)
      
      median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
      yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
      yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:nrow(Data_survival_curve) * as.integer(length(yhat_weibull_sort_total_exp)/nrow(Data_survival_curve))]
      
      Data_survival_curve$left_adj = yhat_weibull_sort_average_total_exp
      Data_survival_curve$left_adj[Data_survival_curve$left_adj > 10000] = 10000
      
      independence_check = with(Data_survival_curve, cKendall(
        trun = time_palliative_enroll,
        obs = observed_OS,
        delta = censor,
        trans = "linear"))
      
      simulation_fit = survfit(Surv(event = rep(1,nrow(Data_survival_curve)),
                                    time = left_adj) ~ 1, 
                               data = Data_survival_curve)
      survival_fit_true = survfit(Surv(true_OS, true_censor) ~ 1, data = Data_survival %>% dplyr::filter(true_OS > immortal_length))
      #survival_fit_observed = survfit(Surv(observed_OS, censor) ~ 1, data = Data_survival %>% dplyr::filter(true_OS > immortal_length & observed_OS > 0))
      survival_fit_CGP = survfit(Surv(observed_OS, censor) ~ 1, data = Data_survival %>% dplyr::filter(CGP_undergo == 1 & observed_OS > immortal_length))
      survival_fit_pre_CGP = survfit(Surv(time_palliative_enroll, true_censor) ~ 1, data = Data_survival %>% dplyr::filter(CGP_undergo == 1 & observed_OS > immortal_length))
      survival_fit_post_CGP = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival %>% dplyr::filter(CGP_undergo == 1 & observed_OS > immortal_length))
      survival_fit_CGP_risk_set = survfit(Surv(time = time_palliative_enroll,time2 = observed_OS, censor) ~ 1, data = Data_survival %>% dplyr::filter(CGP_undergo == 1 & observed_OS > immortal_length))
      g = ggsurvplot(
        fit = list(survival_fit_true = survival_fit_true,
                   simulation_fit = simulation_fit,
                   survival_fit_CGP = survival_fit_CGP,
                   survival_fit_pre_CGP = survival_fit_pre_CGP,
                   survival_fit_post_CGP = survival_fit_post_CGP,
                   survival_fit_CGP_risk_set = survival_fit_CGP_risk_set),
        combine = TRUE,
        # data = Data_survival,
        xlab = "simulation data, time from treatment start date (months)",
        ylab = "Survival Probability",
        censor = TRUE,
        conf.int = FALSE,
        surv.scale = "percent",
        font.title = 8,
        font.subtitle = 8,
        font.main = 8,
        font.submain = 8,
        font.caption = 8,
        font.legend = 8,
        cumevents = FALSE,
        cumcensor = FALSE,
        xscale = "d_m",
        pval = FALSE,
        surv.median.line = "hv",
        palette = "Dark2",
        risk.table = TRUE,
        risk.table.y.text = FALSE,
        tables.theme = clean_theme(),
        legend = c(0.8,0.90),
        xlim = c(0,2000),
        break.x.by = 365.25 * 2.5,
        legend.labs = c("Ideal data: all patients without censoring",
                        "Estimated survival with bias adjustment method",
                        "Conventional Kaplan-Meier estimater",
                        "Suvival before CGP",
                        "Suvival after CGP",
                        "Risk-set adjusted survival")
      ) +   
        labs(
          title = paste0(
            "Toy survival data with left-truncation bias, cKendall tau=",
            format_p(independence_check$PE, digits = 3),
            ", Median OS: Ideal data, ",
            format_p(digits = 1, summary(survival_fit_true)$table[[7]] / 365.25 * 12), " (",
            format_p(digits = 1, summary(survival_fit_true)$table[[8]] / 365.25 * 12), "-",
            format_p(digits = 1, summary(survival_fit_true)$table[[9]] / 365.25 * 12), "), months (95%CI)"
          ),
          subtitle = paste0(
            "bias-corrected, ",
            format_p(digits = 1, summary(simulation_fit)$table[[7]] / 365.25 * 12), " (",
            format_p(digits = 1, summary(simulation_fit)$table[[8]] / 365.25 * 12), "-",
            format_p(digits = 1, summary(simulation_fit)$table[[9]] / 365.25 * 12), "), ",
            "risk-set adjusted, ",
            format_p(digits = 1, summary(survival_fit_CGP_risk_set)$table[[7]] / 365.25 * 12), " (",
            format_p(digits = 1, summary(survival_fit_CGP_risk_set)$table[[8]] / 365.25 * 12), "-",
            format_p(digits = 1, summary(survival_fit_CGP_risk_set)$table[[9]] / 365.25 * 12), "), ",
            "conventional KM estimator, ",
            format_p(digits = 1, summary(survival_fit_CGP)$table[[7]] / 365.25 * 12), " (",
            format_p(digits = 1, summary(survival_fit_CGP)$table[[8]] / 365.25 * 12), "-",
            format_p(digits = 1, summary(survival_fit_CGP)$table[[9]] / 365.25 * 12), ")"
          )
        )
      g$table <- g$table + theme(plot.title = element_blank(),
                                 plot.subtitle = element_blank())
      incProgress(1 / 3)
      OUTPUT_DATA$figure_bias_correction_simulation_g = g
    })
  })
  rm(analysis_env)
  gc()
})

output$figure_bias_correction_simulation = renderPlot({
  req(OUTPUT_DATA$figure_bias_correction_simulation_g)
  OUTPUT_DATA$figure_bias_correction_simulation_g
})

output$text_bias_correction_simulation = renderText({
  req(input$scale_llogis)
  paste0("Median survival after chemotherapy induction is ", input$scale_llogis, " years",
         "\nSurvival after induction of chemotherapy follows a log-logistic distribution with shape=", input$shape_llogis,", scale=", input$scale_llogis, " years",
         "\nSurvival after CGP follows a log-logistic distribution with shape=", input$shape_CGP, ", scale=", input$scale_CGP, " years",
         "\nMaximum follow-up period after CGP testing follows a Weibull distribution with shape=", input$shape_weibull, ", scale=", input$scale_weibull, " years for remaining survival",
         "\nPatients surviving less than 1 months are not offered chemotherapy",
         "\nPatients with survival less than 1 months do not undergo CGP testing\n",
         input$censor_rate*100, "% of cases are right-censored at random timing after CGP testing"
  )
})
