survival_CTx_analysis_logic <- function() {
  analysis_env <- new.env()
  clear_reactive_data(OUTPUT_DATA)
  gc()
  with(analysis_env, {
    withProgress(message = sample(nietzsche)[1], {
      Data_case_target = Data_case()
      if(!is.null(input$gene_group_analysis) && input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed", "Only cases without mutations in the gene set are analyzed")){
        Data_case_target = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                          Data_report()$Tumor_Sample_Barcode)
      }
      Data_case_target = Data_case_target %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                      症例.基本情報.年齢,
                      Lymph_met,
                      Brain_met,
                      Lung_met,
                      Bone_met,
                      Liver_met,
                      Other_met,
                      EP_option,
                      EP_treat,
                      YoungOld,
                      症例.基本情報.性別.名称.,
                      症例.基本情報.がん種.OncoTree.,
                      症例.基本情報.がん種.OncoTree..名称.,
                      症例.基本情報.がん種.OncoTree.LEVEL1.,
                      症例.検体情報.パネル.名称.,
                      症例.背景情報.ECOG.PS.名称.,
                      症例.背景情報.喫煙歴有無.名称.,
                      症例.背景情報.アルコール多飲有無.名称.,
                      症例.背景情報.重複がん有無.異なる臓器..名称.,
                      症例.背景情報.多発がん有無.同一臓器..名称.,
                      final_observe,
                      censor,
                      CTx_lines_before_CGP,
                      pre_CGP_best_RECIST,
                      treat_group,
                      treat_group_2,
                      treat_group_3,
                      time_enroll_final,
                      time_palliative_final,
                      time_palliative_enroll,
                      time_2L_final,
                      time_diagnosis_enroll,
                      time_diagnosis_final,
                      time_2L_enroll
        )
      incProgress(1 / 13)

      Data_case_target$Cancers = Data_case_target$症例.基本情報.がん種.OncoTree.
      Data_case_target$enroll_censor = 1
      Data_case_target = Data_case_target %>%
        dplyr::distinct(.keep_all = TRUE, C.CAT調査結果.基本項目.ハッシュID)
      Data_drug = Data_drug_raw()
      OUTPUT_DATA$figure_surv_Bayes_Data_drug = Data_drug

      Data_MAF = Data_report() %>%
        dplyr::filter(
          !str_detect(Hugo_Symbol, ",") &
            Hugo_Symbol != "" &
            Evidence_level %in% c("","A","B","C","D","E","F") &
            Variant_Classification != "expression"
        ) %>%
        dplyr::arrange(desc(Evidence_level)) %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        Hugo_Symbol,
                        Start_Position,
                        .keep_all = TRUE)
      if(length(Data_MAF[Data_MAF$TMB > 30,]$TMB) > 0){
        Data_MAF[Data_MAF$TMB > 30,]$TMB = 30
      }

      Data_report_TMB = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
        dplyr::select(Tumor_Sample_Barcode, TMB) %>%
        dplyr::distinct(Tumor_Sample_Barcode, .keep_all=T)
      colnames(Data_report_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB")
      Data_case_target = left_join(Data_case_target, Data_report_TMB,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_cluster_ID_list = Data_cluster_ID() %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster)
      Data_case_target = left_join(Data_case_target,
                                   Data_cluster_ID_list,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_case_target$cluster[is.na(Data_case_target$cluster)] = max(Data_case_target$cluster, na.rm = T) + 1
      Data_MAF_target = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
      if(input$patho == "Only pathogenic muts"){
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Evidence_level == "F")
      }
      Data_MAF_target = Data_MAF_target %>%
        dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, .keep_all = T)
      incProgress(1 / 13)
      # Gene_list = unique(names(sort(table(Data_MAF_target$Hugo_Symbol),
      #                               decreasing = T)))
      Data_case_target$diagnosis = Data_case_target$Cancers
      if(nrow(Data_case_target)>0){
        if(length(Data_case_target$censor[is.na(Data_case_target$censor)])> 0){
          Data_case_target$censor[is.na(Data_case_target$censor)] = 0
        }
        Data_survival_interactive = Data_case_target
        candidate_genes = sort(unique(c(Data_MAF_target$Hugo_Symbol,
                                        paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                                        paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                                        paste0(input$special_gene, "_NOS"))))
        candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
        candidate_genes = candidate_genes[!is.na(candidate_genes)]
        Top_gene = unique(names(sort(table(Data_MAF_target$Hugo_Symbol), decreasing = T)))
        Top_gene = Top_gene[Top_gene %in% candidate_genes]
        candidate_drugs = sort(unique(c(Data_drug$Drug)))
        OUTPUT_DATA$figure_surv_interactive_Top_gene_Bayes = Top_gene
        OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes = candidate_genes
        OUTPUT_DATA$figure_surv_interactive_candidate_drugs_Bayes = candidate_drugs[!is.na(candidate_drugs)]
        OUTPUT_DATA$figure_surv_interactive_candidate_lines_Bayes = sort(unique(Data_survival_interactive$CTx_lines_before_CGP))
        OUTPUT_DATA$figure_surv_interactive_candidate_RECIST_Bayes = sort(unique(Data_survival_interactive$pre_CGP_best_RECIST))
        OUTPUT_DATA$figure_surv_interactive_candidate_Age_Bayes = sort(unique(Data_survival_interactive$YoungOld))
        OUTPUT_DATA$figure_surv_interactive_candidate_Sex_Bayes = sort(unique(Data_survival_interactive$症例.基本情報.性別.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_Histology_Bayes = sort(unique(Data_survival_interactive$Cancers))
        OUTPUT_DATA$figure_surv_interactive_candidate_cluster_Bayes = sort(unique(Data_survival_interactive$cluster))
        OUTPUT_DATA$figure_surv_interactive_candidate_PS_Bayes = sort(unique(Data_survival_interactive$症例.背景情報.ECOG.PS.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_meta_Bayes = c('Lymph_met','Brain_met','Lung_met','Bone_met','Liver_met')
        OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target = Data_MAF_target
        OUTPUT_DATA$figure_surv_Bayes_Data_case_target = Data_case_target
      }
    })
  })
  rm(analysis_env)
  gc()
}

observeEvent(input$start_Bayes_compare_prediction, {
  # analysis for common oncogenic mutations
  req(OUTPUT_DATA$figure_surv_Bayes_Data_case_target,
      OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target,
      OUTPUT_DATA$figure_surv_Bayes_Data_drug)
  Data_case_target = OUTPUT_DATA$figure_surv_Bayes_Data_case_target
  Data_MAF_target = OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target
  Data_drug = OUTPUT_DATA$figure_surv_Bayes_Data_drug
  Bayes_env <- new.env()
  gc()
  with(Bayes_env, {
    withProgress(message = sample(nietzsche)[1], {
      if(input$Bayes_prediction_var == "diagnosis"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_diagnosis_enroll
      } else if(input$Bayes_prediction_var == "1L"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_palliative_enroll
      } else if(input$Bayes_prediction_var == "2L"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_2L_enroll
      }
      Data_survival_interactive = Data_case_target %>% dplyr::filter(
        !is.na(time_palliative_enroll) &
          !is.na(time_enroll_final) &
          is.finite(time_palliative_enroll) &
          is.finite(time_enroll_final) &
          time_palliative_enroll > 0 &
          time_enroll_final > 0 &
          !is.na(censor) &
          is.finite(censor)
      )
      max_samples = isolate(input$figure_Bayes_max_samples)
      max_samples = ifelse(is.null(max_samples), 2000, max_samples)
      if (nrow(Data_survival_interactive) > max_samples){
        sample_indices = sample(nrow(Data_survival_interactive),max_samples)
        Data_survival_interactive = Data_survival_interactive[sample_indices, ]
      } else {
        Data_survival_interactive = Data_survival_interactive
      }
      # 転移部位のマッピングを定義
      metastasis_mapping <- c(
        "Lymph_met" = "Lymph_met",
        "Brain_met" = "Brain_met",
        "Lung_met" = "Lung_met",
        "Bone_met" = "Bone_met",
        "Liver_met" = "Liver_met"
      )

      # ID抽出関数を定義
      extract_group_ids <- function(group_num) {
        # 初期IDセット
        IDs <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)

        # 動的に入力名を構築
        input_prefix <- paste0("gene_survival_interactive_", group_num, "_")
        input_postfix <- "_Bayes"
        # 遺伝子フィルタ（P_1: 必須変異1）
        p1_input <- input[[paste0(input_prefix, "P_1", input_postfix)]]
        if(!all(is.null(p1_input))) {
          IDs <- intersect(IDs, (Data_MAF_target %>%
                                   dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
        }

        # 遺伝子フィルタ（P_2: 必須変異2）
        p2_input <- input[[paste0(input_prefix, "P_2", input_postfix)]]
        if(!all(is.null(p2_input))) {
          IDs <- intersect(IDs, (Data_MAF_target %>%
                                   dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
        }

        # 遺伝子除外（W: 除外変異）
        w_input <- input[[paste0(input_prefix, "W", input_postfix)]]
        if(!all(is.null(w_input))) {
          IDs <- setdiff(IDs, (Data_MAF_target %>%
                                 dplyr::filter(Hugo_Symbol %in% w_input))$Tumor_Sample_Barcode)
        }

        # 臨床データフィルタ
        clinical_filters <- list(
          L = "CTx_lines_before_CGP",
          R = "pre_CGP_best_RECIST",
          A = "YoungOld",
          S = "症例.基本情報.性別.名称.",
          H = "Cancers",
          C = "cluster"
        )

        for(filter_key in names(clinical_filters)) {
          filter_input <- input[[paste0(input_prefix, filter_key, input_postfix)]]
          if(!all(is.null(filter_input))) {
            column_name <- clinical_filters[[filter_key]]
            filter_expr <- paste0(column_name, " %in% filter_input")
            IDs <- intersect(IDs, (Data_survival_interactive %>%
                                     dplyr::filter(!!parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
          }
        }

        # 転移部位フィルタ（M: 複数選択可能、OR条件）
        m_input <- input[[paste0(input_prefix, "M", input_postfix)]]
        if(!all(is.null(m_input))) {
          met_ids <- NULL
          for(met_type in m_input) {
            if(met_type %in% names(metastasis_mapping)) {
              column_name <- metastasis_mapping[[met_type]]
              filter_expr <- paste0(column_name, " == 'Yes'")
              met_ids <- union(met_ids, (Data_survival_interactive %>%
                                           dplyr::filter(!!parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
            }
          }
          IDs <- intersect(IDs, met_ids)
        }

        # 薬剤フィルタ（D）
        d_input <- input[[paste0(input_prefix, "D", input_postfix)]]
        if(!all(is.null(d_input))) {
          IDs <- intersect(IDs, (Data_drug %>%
                                   dplyr::filter(Drug %in% d_input))$ID)
        }
        nd_input <- input[[paste0(input_prefix, "ND", input_postfix)]]
        if(!all(is.null(nd_input))) {
          IDs <- setdiff(IDs, (Data_drug %>%
                                   dplyr::filter(Drug %in% nd_input))$ID)
        }

        return(IDs)
      }

      # Group1とGroup2のIDを取得
      ID_1 <- extract_group_ids(1)
      ID_2 <- extract_group_ids(2)

      run_custom_survival_analysis(
        Data_survival = Data_survival_interactive,
        id_set_1 = ID_1,
        id_set_2 = ID_2,
        label_1 = "Group 1",
        label_2 = "Group 2"
      )
    })
  })
  rm(Bayes_env)
  gc()
})



observeEvent(input$start_Bayes_basic_prediction, {
  # analysis for common oncogenic mutations
  req(OUTPUT_DATA$figure_surv_Bayes_Data_case_target,
      OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target)
  Data_case_target = OUTPUT_DATA$figure_surv_Bayes_Data_case_target
  Data_MAF_target = OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target
  Bayes_env <- new.env()
  gc()
  with(Bayes_env, {
    withProgress(message = sample(nietzsche)[1], {
      if(input$Bayes_basic_var == "diagnosis"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_diagnosis_enroll
      } else if(input$Bayes_basic_var == "1L"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_palliative_enroll
      } else if(input$Bayes_basic_var == "2L"){
        Data_case_target$time_palliative_enroll = Data_case_target$time_2L_enroll
      }
      Data_survival = Data_case_target %>% dplyr::filter(
        !is.na(time_palliative_enroll) &
          !is.na(time_enroll_final) &
          is.finite(time_palliative_enroll) &
          is.finite(time_enroll_final) &
          time_palliative_enroll > 0 &
          time_enroll_final > 0 &
          !is.na(censor) &
          is.finite(censor)
      )
      max_samples = isolate(input$figure_Bayes_max_samples)
      max_samples = ifelse(is.null(max_samples), 2000, max_samples)
      if (nrow(Data_survival) > max_samples){
        sample_indices = sample(nrow(Data_survival),max_samples)
        Data_survival = Data_survival[sample_indices, ]
      } else {
        Data_survival = Data_survival
      }
      Data_survival_curve = Data_survival
      if(sum(Data_survival_curve$censor) > 0){
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

      if(sum(Data_survival_curve$censor) > 0 ){
        fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_curve)
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
      Data_MAF_target = Data_MAF_target %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_survival_curve$C.CAT調査結果.基本項目.ハッシュID)

      Data_drug_survival_1 = Data_drug_survival_1 %>%
        dplyr::mutate(delay_1 = case_when(
          Data_drug_survival_1$time_palliative_enroll <= Median_entry ~ "SPT to entry early",
          TRUE ~ "SPT to entry later"))
      Cancername = sort(unique(Data_survival_curve$Cancers))
      Total_pts = as.vector(table((Data_survival_curve %>%
                                     dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = T))$Cancers))

      incProgress(1 / 13)
      stan_weibull_survival_model_data <-
        list(
          Median = as.integer(Median_entry),
          Nobs_early = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later"),
          Nobs_late = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later"),
          Ncen_early = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later"),
          Ncen_late = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later"),
          Nexp = length(Data_drug_survival_2$time_palliative_enroll),
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
          iter_sampling = ITER,  # (iter-warmup)
          iter_warmup = as.integer(ITER / 2),adapt_delta=0.99, max_treedepth = 10,
          thin = 1,
          refresh = 500
        )
      } else {
        mod <- readRDS(stan_model_simple_save_path_rstan)
        stan_weibull_survival_model_fit <-
          rstan::sampling(object = mod,
                      data = stan_weibull_survival_model_data,
                      control=list(adapt_delta=0.99, max_treedepth = 10),
                      seed=1234, chains=4, iter=as.integer(ITER * 3 / 2),warmup=as.integer(ITER / 2), thin=1, refresh=500, verbose = FALSE)
      }
      incProgress(1 / 13)

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

      median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=4 * ITER)
      median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)

      median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
      yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
      yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total_exp)/length(Data_survival_curve$censor))]

      Data_survival_curve$left_adj = yhat_weibull_sort_average_total_exp
      Data_survival_curve$left_adj[Data_survival_curve$left_adj > 10000] = 10000


      independence_check = with(Data_survival_curve, cKendall(
        trun = time_palliative_enroll,
        obs = time_palliative_final,
        delta = censor,
        trans = "linear"))


      traditional_fit = survfit(Surv(event = censor,
                                     time = time_palliative_final) ~ 1,
                                data = Data_survival_curve)
      delayed_entry_fit = survfit(Surv(event = censor,
                                       time = time_palliative_enroll,
                                       time2 = time_palliative_final) ~ 1,
                                  data = Data_survival_curve)
      simulation_fit_2 = survfit(Surv(event = rep(1,length(Data_survival_curve$left_adj)),
                                      time = left_adj) ~ 1,
                                 data = Data_survival_curve)
      hsb_base = Data_survival_curve
      hsb_base$diagnosis = " ALL"
      OUTPUT_DATA$figure_surv_Bayes_hsb_base = hsb_base
      OUTPUT_DATA$figure_surv_Bayes_g = list()
      OUTPUT_DATA$figure_surv_Bayes_g[[1]] = ggsurvplot(
        fit = list(traditional_fit = traditional_fit,
                   delayed_entry_fit = delayed_entry_fit,
                   simulation_fit_2 = simulation_fit_2
        ),
        combine = TRUE,
        data = Data_survival_curve,
        xlab = "Time from CTx induction to final observation (months)",
        ylab = "Survival Probability",
        censor = TRUE,
        conf.int = FALSE,
        font.title = 8,
        font.subtitle = 8,
        font.main = 8,
        font.submain = 8,
        font.caption = 8,
        font.legend = 8,
        pval = FALSE,
        surv.median.line = "hv",
        risk.table = TRUE,
        risk.table.y.text = FALSE,
        tables.theme = clean_theme(),
        legend = c(0.8,0.90),
        cumevents = FALSE,
        cumcensor = FALSE,
        break.x.by = 365.25 * 2.5,
        xscale = "d_m",
        xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
        #break.x.by = ceiling(max(Data_survival_curve$time_palliative_final) / 500) * 100,
        legend.labs = c("All patients, Unadjusted",
                        "All patients, Number at risk adjusted",
                        "All patients, Adj. right cens. & left trunc.")
      ) +
        labs(title = paste(traditional_fit[[1]], " patients; Median OS ",
                           format_p(summary(traditional_fit)$table[[7]] / 365.25 * 12, digits = 1)," (",
                           format_p(summary(traditional_fit)$table[[8]] / 365.25 * 12, digits = 1),"-",
                           format_p(summary(traditional_fit)$table[[9]] / 365.25 * 12, digits = 1),") (unadj.)/",
                           format_p(summary(delayed_entry_fit)$table[[7]] / 365.25 * 12, digits = 1), " (",
                           format_p(summary(delayed_entry_fit)$table[[8]] / 365.25 * 12, digits = 1),"-",
                           format_p(summary(delayed_entry_fit)$table[[9]] / 365.25 * 12, digits = 1),") (risk adj.)/",
                           "\n",
                           format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                           format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                           format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                           ") (adj. right & left) months",
                           sep=""),
             subtitle = paste(2, "-year rate ",
                              format_p(summary(traditional_fit, time = 2 * 365.25)$surv * 100,digits=1),
                              "% (unadj.)/",
                              format_p(summary(delayed_entry_fit, time = 2 * 365.25)$surv * 100,digits=1),
                              "% (risk adj.)/",
                              format_p(summary(simulation_fit_2, time = 2 * 365.25)$surv * 100,digits=1),
                              "% (adj. right & left), ",
                              "cKendall tau= ", format_p(independence_check$PE,digits=2),
                              ", SE= ", format_p(independence_check$SE,digits=2),
                              ", p= ", format_p(independence_check$p.value,digits=2), sep=""))
      OUTPUT_DATA$figure_surv_Bayes_g[[1]]$table <- OUTPUT_DATA$figure_surv_Bayes_g[[1]]$table + theme(plot.title = element_blank(),
                                           plot.subtitle = element_blank())
      OUTPUT_DATA$figure_surv_Bayes_g[[2]] = ggsurvplot(
        fit = list(traditional_fit = traditional_fit,
                   # delayed_entry_fit = delayed_entry_fit,
                   simulation_fit_2 = simulation_fit_2
        ),
        combine = TRUE,
        data = Data_survival_curve,
        xlab = "Time from CTx induction to final observation (months)",
        ylab = "Survival Probability",
        censor = TRUE,
        conf.int = FALSE,
        font.title = 8,
        font.subtitle = 8,
        font.main = 8,
        font.submain = 8,
        font.caption = 8,
        font.legend = 8,
        pval = FALSE,
        surv.median.line = "hv",
        risk.table = TRUE,
        risk.table.y.text = FALSE,
        tables.theme = clean_theme(),
        legend = c(0.8,0.90),
        cumevents = FALSE,
        cumcensor = FALSE,
        break.x.by = 365.25 * 2.5,
        xscale = "d_m",
        xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
        #break.x.by = ceiling(max(Data_survival_curve$time_palliative_final) / 500) * 100,
        legend.labs = c("All patients, Unadjusted",
                        # "All patients, Number at risk adjusted",
                        "All patients, Adj. right cens. & left trunc.")
      ) +
        labs(title = paste(traditional_fit[[1]], " patients; Median OS ",
                           format_p(summary(traditional_fit)$table[[7]] / 365.25 * 12, digits = 1)," (",
                           format_p(summary(traditional_fit)$table[[8]] / 365.25 * 12, digits = 1),"-",
                           format_p(summary(traditional_fit)$table[[9]] / 365.25 * 12, digits = 1),") (unadj.)/",
                           "\n",
                           format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                           format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                           format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                           ") (adj. right & left) months",
                           sep=""),
             subtitle = paste(2, "-year rate ",
                              format_p(summary(traditional_fit, time = 2 * 365.25)$surv * 100,digits=1),
                              "% (unadj.)/",
                              format_p(summary(simulation_fit_2, time = 2 * 365.25)$surv * 100,digits=1),
                              "% (adj. right & left), ",
                              "cKendall tau= ", format_p(independence_check$PE,digits=2),
                              ", SE= ", format_p(independence_check$SE,digits=2),
                              ", p= ", format_p(independence_check$p.value,digits=2), sep=""))
      OUTPUT_DATA$figure_surv_Bayes_g[[2]]$table <- OUTPUT_DATA$figure_surv_Bayes_g[[2]]$table + theme(plot.title = element_blank(),
                                           plot.subtitle = element_blank())
      OUTPUT_DATA$figure_surv_Bayes_text_median_base = paste0(format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                                format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                                format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                                ")")
      incProgress(1 / 13)
      # Survival analysis for mutation
      mut_gene = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
      if(length(mut_gene)==0){
        mut_gene = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
      }
      ID_mutation = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene))$Tumor_Sample_Barcode
      Data_survival = Data_survival %>% dplyr::mutate(
        mutation = case_when(
          C.CAT調査結果.基本項目.ハッシュID %in% ID_mutation ~ paste(mut_gene, collapse = ";"),
          TRUE ~ "no-mut"
        )
      )
      traditional_fit = survfit(Surv(event = censor,
                                     time = time_palliative_final) ~ mutation,
                                data = Data_survival)
      if(length(Data_survival[,1]) != traditional_fit[[1]][[1]]){
        name_1 = paste("mutation", sort(unique(Data_survival$mutation))[[1]])
        name_2 = paste("mutation", sort(unique(Data_survival$mutation))[[2]])
        name_1_short = sort(unique(Data_survival$mutation))[[1]]
        name_2_short = sort(unique(Data_survival$mutation))[[2]]

        Data_survival$gene_mut = Data_survival$mutation == sort(unique(Data_survival$mutation))[[2]]
        Data_survival_curve = Data_survival
        Data_survival_tmp_pos = Data_survival_curve %>%
          dplyr::filter(gene_mut == TRUE)
        if(sum(Data_survival_tmp_pos$censor) > 0 ){
          fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_pos)
          time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_enroll_final)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 1
          survival_rate_90 = 0
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_enroll_final)*0.9)
        Data_survival_tmp = Data_survival_tmp_pos %>% dplyr::filter(
          time_enroll_final >= time_10 &
            time_enroll_final <= time_90)
        Data_survival_tmp$time_enroll_final =
          Data_survival_tmp$time_enroll_final - (time_10 - 1)
        Data_survival_tmp2 = Data_survival_tmp_pos %>% dplyr::filter(
          time_enroll_final > time_90) %>%
          dplyr::arrange(time_enroll_final)

        Data_survival_tmp2$time_enroll_final = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp2)
        Data_survival_tmp = rbind(Data_survival_tmp, Data_survival_tmp2[idx[idx < max(idx)*(1-survival_rate_10)],])

        Data_survival_tmp_neg = Data_survival_curve %>%
          dplyr::filter(gene_mut == FALSE)
        if(sum(Data_survival_tmp_neg$censor) > 0 ){
          fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_neg)
          time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_enroll_final)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 1
          survival_rate_90 = 0
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_enroll_final)*0.9)
        Data_survival_tmp3 = Data_survival_tmp_neg %>% dplyr::filter(
          time_enroll_final >= time_10 &
            time_enroll_final <= time_90)
        Data_survival_tmp3$time_enroll_final =
          Data_survival_tmp3$time_enroll_final - (time_10 - 1)
        Data_survival_tmp4 = Data_survival_tmp_neg %>% dplyr::filter(
          time_enroll_final > time_90) %>%
          dplyr::arrange(time_enroll_final)

        Data_survival_tmp4$time_enroll_final = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp4)
        Data_survival_tmp3 = rbind(Data_survival_tmp3, Data_survival_tmp4[idx[idx < max(idx)*(1-survival_rate_10)],])
        Data_drug_survival_1 = rbind(Data_survival_tmp, Data_survival_tmp3)

        Data_survival_tmp_pos = Data_survival_curve %>%
          dplyr::filter(gene_mut == TRUE)
        if(sum(Data_survival_tmp_pos$censor) > 0 ){
          fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_pos)
          time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 1
          survival_rate_90 = 0
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)
        Data_survival_tmp5 = Data_survival_tmp_pos %>% dplyr::filter(
          time_palliative_enroll >= time_10 &
            time_palliative_enroll <= time_90)
        Data_survival_tmp5$time_palliative_enroll =
          Data_survival_tmp5$time_palliative_enroll - (time_10 - 1)
        Data_survival_tmp6 = Data_survival_tmp_pos %>% dplyr::filter(
          time_palliative_enroll > time_90) %>%
          dplyr::arrange(time_palliative_enroll)
        Data_survival_tmp6$time_palliative_enroll = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp6)
        Data_survival_tmp5 = rbind(Data_survival_tmp5, Data_survival_tmp6[idx[idx < max(idx)*(1-survival_rate_10)],])

        Data_survival_tmp_neg = Data_survival_curve %>%
          dplyr::filter(gene_mut == FALSE)

        if(sum(Data_survival_tmp_neg$censor) > 0 ){
          fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_neg)
          time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 1
          survival_rate_90 = 0
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)
        Data_survival_tmp7 = Data_survival_tmp_neg %>% dplyr::filter(
          time_palliative_enroll >= time_10 &
            time_palliative_enroll <= time_90)
        Data_survival_tmp7$time_palliative_enroll =
          Data_survival_tmp7$time_palliative_enroll - (time_10 - 1)
        Data_survival_tmp8 = Data_survival_tmp_neg %>% dplyr::filter(
          time_palliative_enroll > time_90) %>%
          dplyr::arrange(time_palliative_enroll)
        Data_survival_tmp8$time_palliative_enroll = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp8)
        Data_survival_tmp7 = rbind(Data_survival_tmp7, Data_survival_tmp8[idx[idx < max(idx)*(1-survival_rate_10)],])

        Median_entry_pos = median(Data_survival_tmp5$time_palliative_enroll)
        Median_entry_neg = median(Data_survival_tmp7$time_palliative_enroll)
        Data_drug_survival_1_2 = rbind(Data_survival_tmp5, Data_survival_tmp7)

        Data_drug_survival_1 = Data_drug_survival_1 %>%
          dplyr::mutate(delay_1 = case_when(
            time_palliative_enroll <= Median_entry_pos &
              gene_mut == TRUE ~ "SPT to entry early",
            time_palliative_enroll > Median_entry_pos &
              gene_mut == TRUE ~ "SPT to entry later",
            time_palliative_enroll <= Median_entry_neg &
              gene_mut == FALSE ~ "SPT to entry early",
            time_palliative_enroll > Median_entry_neg &
              gene_mut == FALSE ~ "SPT to entry later"))

        if(sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 == "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 == "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 != "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 != "SPT to entry later"))$censor) >= 1){

          stan_weibull_survival_model_data <-
            list(
              Median_pos = as.integer(Median_entry_pos),
              Median_neg = as.integer(Median_entry_neg),
              Nobs_early_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Nobs_late_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Ncen_early_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Ncen_late_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Nobs_early_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Nobs_late_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Ncen_early_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Ncen_late_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Nexp_pos = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE]),
              Nexp_neg = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE]),
              ybef_exp_pos = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE],
              ybef_exp_neg = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE],
              yobs_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              yobs_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              ycen_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              ycen_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              yobs_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              yobs_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              ycen_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              ycen_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE]
            )

          if(DOCKER){
            stan_model_compiled <- readRDS(stan_model_factor_save_path_cmdstan)
            stan_weibull_survival_model_fit <- stan_model_compiled$sample(
              data = stan_weibull_survival_model_data,
              seed = 1234,
              chains = 4,
              parallel_chains = PARALLEL,
              iter_sampling = ITER,  # (iter-warmup)
              iter_warmup = as.integer(ITER / 2),adapt_delta=0.99, max_treedepth = 10,
              thin = 1,
              refresh = 500
            )
          } else {
            mod <- readRDS(stan_model_factor_save_path_rstan)
            stan_weibull_survival_model_fit <-
              rstan::sampling(object = mod,
                          data = stan_weibull_survival_model_data,
                          control=list(adapt_delta=0.99, max_treedepth = 10),
                          seed=1234, chains=4, iter=as.integer(ITER * 3 / 2),warmup=as.integer(ITER / 2), thin=1, refresh=500, verbose = FALSE)
          }

          stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_weibull_survival_model_fit)
          stan_weibull_survival_model_draws_total <-
            stan_weibull_survival_model_draws %>%
            dplyr::select(.chain, .iteration, .draw, starts_with("yhat_total")) %>%
            tidyr::gather(key = key, value = yhat_total, starts_with("yhat_total")) %>%
            dplyr::select(-key)
          stan_weibull_survival_model_draws_total_exp <-
            stan_weibull_survival_model_draws %>%
            dplyr::select(.chain, .iteration, .draw, starts_with("yhat_factor_total")) %>%
            tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_factor_total")) %>%
            dplyr::select(-key)

          median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=4 * ITER)
          median_os_list_weibull_total = matrix(stan_weibull_survival_model_draws_total$yhat_total, nrow=4 * ITER)
          median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)
          median_os_list_weibull_total = rowMedians(median_os_list_weibull_total,na.rm = T)
          median_os_list_weibull_bias = median_os_list_weibull_total - median_os_list_weibull_total_exp

          median_os_summary_weibull_total = quantile(median_os_list_weibull_total, probs = c(0.025, 0.5, 0.975))
          median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
          median_os_summary_weibull_bias = quantile(median_os_list_weibull_bias, probs = c(0.025, 0.5, 0.975))
          yhat_weibull_sort_total = sort(stan_weibull_survival_model_draws_total$yhat_total)
          yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
          yhat_weibull_sort_average_total = yhat_weibull_sort_total[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total)/length(Data_survival_curve$censor))]
          yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total_exp)/length(Data_survival_curve$censor))]

          Data_survival_curve$factor_neg = yhat_weibull_sort_average_total
          Data_survival_curve$factor_pos = yhat_weibull_sort_average_total_exp
          Data_survival_curve$factor_neg[Data_survival_curve$factor_neg > 10000] = 10000
          Data_survival_curve$factor_pos[Data_survival_curve$factor_pos > 10000] = 10000

          traditional_fit = survfit(Surv(event = censor,
                                         time = time_palliative_final) ~ gene_mut,
                                    data = Data_survival_curve)
          simulation_fit_neg = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_neg)),
                                            time = factor_neg) ~ 1,
                                       data = Data_survival_curve)
          simulation_fit_pos = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_pos)),
                                            time = factor_pos) ~ 1,
                                       data = Data_survival_curve)

          OUTPUT_DATA$figure_surv_Bayes_g[[3]] = ggsurvplot(
            fit = list(traditional_fit = traditional_fit,
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
            legend = c(0.8,0.8),
            xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
            cumevents = FALSE,
            cumcensor = FALSE,
            break.x.by = 365.25 * 2.5,
            xscale = "d_m",
            #break.x.by = ceiling(max(Data_survival_curve$time_palliative_final) / 500) * 100,
            legend.labs = c(paste0(name_1, ", Unadjusted"),
                            paste0(name_2, ", Unadjusted"),
                            paste0(name_1, ", Adjusted for Delayed Entry"),
                            paste0(name_2, ", Adjusted for Delayed Entry"))
          ) +
            labs(title = paste(sum(traditional_fit[[1]]), " patients ",
                               "Median OS ",
                               format_p(summary(traditional_fit)$table[[13]] / 365.25 * 12, digits = 1)," (",
                               format_p(summary(traditional_fit)$table[[15]] / 365.25 * 12, digits = 1),"-",
                               format_p(summary(traditional_fit)$table[[17]] / 365.25 * 12, digits = 1),") (",
                               name_1_short,", unadj.)/",
                               format_p(summary(traditional_fit)$table[[14]] / 365.25 * 12, digits = 1), " (",
                               format_p(summary(traditional_fit)$table[[16]] / 365.25 * 12, digits = 1),"-",
                               format_p(summary(traditional_fit)$table[[18]] / 365.25 * 12, digits = 1),") (",
                               name_2_short,", unadj.)/",
                               format_p(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1),
                               " (", format_p(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1), "-",
                               format_p(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1), ") (",
                               name_1_short,", adj.)/",
                               format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                               format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                               format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                               ") (", name_2_short,", adj.) months",
                               sep=""),
                 subtitle = paste("Survival difference, ",name_1_short, "-", name_2_short, ": ",
                                  format_p(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 1), " (",
                                  format_p(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 1), "-",
                                  format_p(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 1),
                                  ") months, median (95% CI)",
                                  sep=""))
          OUTPUT_DATA$figure_surv_Bayes_g[[3]]$table <- OUTPUT_DATA$figure_surv_Bayes_g[[3]]$table + theme(plot.title = element_blank(),
                                               plot.subtitle = element_blank())
        }
      }
    })
  })
  rm(Bayes_env)
  gc()
})

observeEvent(input$start_Bayes_mutation_prediction, {
  # analysis for common oncogenic mutations
  req(input$Bayes_mutation_no,
      input$Bayes_mutation_var,
      input$gene,
      OUTPUT_DATA$figure_surv_Bayes_Data_case_target,
      OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target)
  Data_case_target = OUTPUT_DATA$figure_surv_Bayes_Data_case_target
  Data_MAF_target = OUTPUT_DATA$figure_surv_Bayes_Data_MAF_target
  Bayes_env <- new.env()
  gc()
  with(Bayes_env, {
    if(input$Bayes_mutation_var == "diagnosis"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_diagnosis_enroll
    } else if(input$Bayes_mutation_var == "1L"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_palliative_enroll
    } else if(input$Bayes_mutation_var == "2L"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_2L_enroll
    }
    Data_survival = Data_case_target %>% dplyr::filter(
      !is.na(time_palliative_enroll) &
        !is.na(time_enroll_final) &
        is.finite(time_palliative_enroll) &
        is.finite(time_enroll_final) &
        time_palliative_enroll > 0 &
        time_enroll_final > 0 &
        !is.na(censor) &
        is.finite(censor)
    )
    max_samples = isolate(input$figure_Bayes_max_samples)
    max_samples = ifelse(is.null(max_samples), 2000, max_samples)
    if (nrow(Data_survival) > max_samples){
      sample_indices = sample(nrow(Data_survival),max_samples)
      Data_survival = Data_survival[sample_indices, ]
    } else {
      Data_survival = Data_survival
    }
    h=list()
    hs=list()
    k = 1
    oncogenic_genes = Data_MAF_target %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
      dplyr::distinct() %>%
      dplyr::count(Hugo_Symbol) %>%
      dplyr::arrange(-n)
    colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
    oncogenic_genes = rbind(oncogenic_genes %>% dplyr::filter(gene_mutation %in% input$gene),
                            oncogenic_genes %>% dplyr::filter(!gene_mutation %in% input$gene))
    oncogenic_genes = oncogenic_genes[1:min(length(oncogenic_genes$all_patients), input$Bayes_mutation_no),]

    gene_table = data.frame(oncogenic_genes$gene_mutation)
    colnames(gene_table) = c("Gene")
    gene_table$positive_patients = oncogenic_genes$all_patients
    gene_table$positive_freq = format_p(100 * gene_table$positive_patients / length(Data_survival$C.CAT調査結果.基本項目.ハッシュID), digits = 1)
    gene_table$positive_median = 0
    gene_table$positive_LL = 0
    gene_table$positive_UL = 0
    gene_table$negative_median = 0
    gene_table$negative_LL = 0
    gene_table$negative_UL = 0
    gene_table$diff_median = 0
    gene_table$diff_LL = 0
    gene_table$diff_UL = 0
    gene_table$positive_mortal_event = 0
    gene_table$negative_mortal_event = 0
    withProgress(message = "Oncogenic variant and overall survival", {
      for(i in 1:nrow(oncogenic_genes)){
        ID_mutation = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% oncogenic_genes$gene_mutation[i]))$Tumor_Sample_Barcode
        Data_survival = Data_survival %>% dplyr::mutate(
          gene_mut = case_when(
            C.CAT調査結果.基本項目.ハッシュID %in% ID_mutation ~ TRUE,
            TRUE ~ FALSE
          )
        )

        mut_gene = paste0(oncogenic_genes$gene_mutation[i],
                          ", top ", i, " gene")
        traditional_fit = survfit(Surv(event = censor,
                                       time = time_palliative_final) ~ gene_mut,
                                  data = Data_survival)
        if(length(Data_survival[,1]) != traditional_fit[[1]][[1]]){
          gene_table$negative_mortal_event[i] = summary(traditional_fit)$table[7]
          gene_table$positive_mortal_event[i] = summary(traditional_fit)$table[8]

          Data_survival_curve = Data_survival
          Data_survival_tmp_pos = Data_survival_curve %>%
            dplyr::filter(gene_mut == TRUE)
          if(sum(Data_survival_tmp_pos$censor) > 0 ){
            fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_pos)
            time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
            survival_rate_10 = summary(fit,times = time_10)[[6]]
            survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_enroll_final)*0.9)[[6]]
          } else {
            time_10 = 90
            survival_rate_10 = 1
            survival_rate_90 = 0
          }
          time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_enroll_final)*0.9)
          Data_survival_tmp = Data_survival_tmp_pos %>% dplyr::filter(
            time_enroll_final >= time_10 &
              time_enroll_final <= time_90)
          Data_survival_tmp$time_enroll_final =
            Data_survival_tmp$time_enroll_final - (time_10 - 1)
          Data_survival_tmp2 = Data_survival_tmp_pos %>% dplyr::filter(
            time_enroll_final > time_90) %>%
            dplyr::arrange(time_enroll_final)

          Data_survival_tmp2$time_enroll_final = time_90 - (time_10 - 1)
          idx <- 1:nrow(Data_survival_tmp2)
          Data_survival_tmp = rbind(Data_survival_tmp, Data_survival_tmp2[idx[idx < max(idx)*(1-survival_rate_10)],])

          Data_survival_tmp_neg = Data_survival_curve %>%
            dplyr::filter(gene_mut == FALSE)
          if(sum(Data_survival_tmp_neg$censor) > 0){
            fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_neg)
            time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
            survival_rate_10 = summary(fit,times = time_10)[[6]]
            survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_enroll_final)*0.9)[[6]]
          } else {
            time_10 = 90
            survival_rate_10 = 0
            survival_rate_90 = 1
          }
          time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_enroll_final)*0.9)
          Data_survival_tmp3 = Data_survival_tmp_neg %>% dplyr::filter(
            time_enroll_final >= time_10 &
              time_enroll_final <= time_90)
          Data_survival_tmp3$time_enroll_final =
            Data_survival_tmp3$time_enroll_final - (time_10 - 1)
          Data_survival_tmp4 = Data_survival_tmp_neg %>% dplyr::filter(
            time_enroll_final > time_90) %>%
            dplyr::arrange(time_enroll_final)

          Data_survival_tmp4$time_enroll_final = time_90 - (time_10 - 1)
          idx <- 1:nrow(Data_survival_tmp4)
          Data_survival_tmp3 = rbind(Data_survival_tmp3, Data_survival_tmp4[idx[idx < max(idx)*(1-survival_rate_10)],])
          Data_drug_survival_1 = rbind(Data_survival_tmp, Data_survival_tmp3)

          Data_survival_tmp_pos = Data_survival_curve %>%
            dplyr::filter(gene_mut == TRUE)
          if(sum(Data_survival_tmp_pos$censor) > 0 ){
            fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_pos)
            time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
            survival_rate_10 = summary(fit,times = time_10)[[6]]
            survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)[[6]]
          } else {
            time_10 = 90
            survival_rate_10 = 1
            survival_rate_90 = 0
          }
          time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)
          Data_survival_tmp5 = Data_survival_tmp_pos %>% dplyr::filter(
            time_palliative_enroll >= time_10 &
              time_palliative_enroll <= time_90)
          Data_survival_tmp5$time_palliative_enroll =
            Data_survival_tmp5$time_palliative_enroll - (time_10 - 1)
          Data_survival_tmp6 = Data_survival_tmp_pos %>% dplyr::filter(
            time_palliative_enroll > time_90) %>%
            dplyr::arrange(time_palliative_enroll)
          Data_survival_tmp6$time_palliative_enroll = time_90 - (time_10 - 1)
          idx <- 1:nrow(Data_survival_tmp6)
          Data_survival_tmp5 = rbind(Data_survival_tmp5, Data_survival_tmp6[idx[idx < max(idx)*(1-survival_rate_10)],])

          Data_survival_tmp_neg = Data_survival_curve %>%
            dplyr::filter(gene_mut == FALSE)

          if(sum(Data_survival_tmp_neg$censor) > 0 ){
            fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_neg)
            time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
            survival_rate_10 = summary(fit,times = time_10)[[6]]
            survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)[[6]]
          } else {
            time_10 = 90
            survival_rate_10 = 1
            survival_rate_90 = 0
          }
          time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)
          Data_survival_tmp7 = Data_survival_tmp_neg %>% dplyr::filter(
            time_palliative_enroll >= time_10 &
              time_palliative_enroll <= time_90)
          Data_survival_tmp7$time_palliative_enroll =
            Data_survival_tmp7$time_palliative_enroll - (time_10 - 1)
          Data_survival_tmp8 = Data_survival_tmp_neg %>% dplyr::filter(
            time_palliative_enroll > time_90) %>%
            dplyr::arrange(time_palliative_enroll)
          Data_survival_tmp8$time_palliative_enroll = time_90 - (time_10 - 1)
          idx <- 1:nrow(Data_survival_tmp8)
          Data_survival_tmp7 = rbind(Data_survival_tmp7, Data_survival_tmp8[idx[idx < max(idx)*(1-survival_rate_10)],])

          Median_entry_pos = median(Data_survival_tmp5$time_palliative_enroll)
          Median_entry_neg = median(Data_survival_tmp7$time_palliative_enroll)
          Data_drug_survival_1_2 = rbind(Data_survival_tmp5, Data_survival_tmp7)

          Data_drug_survival_1 = Data_drug_survival_1 %>%
            dplyr::mutate(delay_1 = case_when(
              time_palliative_enroll <= Median_entry_pos &
                gene_mut == TRUE ~ "SPT to entry early",
              time_palliative_enroll > Median_entry_pos &
                gene_mut == TRUE ~ "SPT to entry later",
              time_palliative_enroll <= Median_entry_neg &
                gene_mut == FALSE ~ "SPT to entry early",
              time_palliative_enroll > Median_entry_neg &
                gene_mut == FALSE ~ "SPT to entry later"))

          if(sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 == "SPT to entry later"))$censor) >= 1 &
             sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 == "SPT to entry later"))$censor) >= 1 &
             sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 != "SPT to entry later"))$censor) >= 1 &
             sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 != "SPT to entry later"))$censor) >= 1){

            stan_weibull_survival_model_data <-
              list(
                Median_pos = as.integer(Median_entry_pos),
                Median_neg = as.integer(Median_entry_neg),
                Nobs_early_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                Nobs_late_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                Ncen_early_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                Ncen_late_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                Nobs_early_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                Nobs_late_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                Ncen_early_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                Ncen_late_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                Nexp_pos = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE]),
                Nexp_neg = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE]),
                ybef_exp_pos = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE],
                ybef_exp_neg = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE],
                yobs_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                yobs_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                ycen_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                ycen_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                yobs_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                yobs_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                ycen_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                ycen_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE]
              )

            if(DOCKER){
              stan_model_compiled <- readRDS(stan_model_factor_save_path_cmdstan)
              stan_weibull_survival_model_fit <- stan_model_compiled$sample(
                data = stan_weibull_survival_model_data,
                seed = 1234,
                chains = 4,
                parallel_chains = PARALLEL,
                iter_sampling = ITER,
                iter_warmup = as.integer(ITER / 2),adapt_delta=0.99, max_treedepth = 10,
                thin = 1,
                refresh = 500
              )
            } else {
              mod <- readRDS(stan_model_factor_save_path_rstan)
              stan_weibull_survival_model_fit <-
                rstan::sampling(object = mod,
                            data = stan_weibull_survival_model_data,
                            control=list(adapt_delta=0.99, max_treedepth = 10),
                            seed=1234, chains=4, iter=as.integer(ITER * 3 / 2),warmup=as.integer(ITER / 2), thin=1, refresh=500, verbose = FALSE)
            }
            stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_weibull_survival_model_fit)
            stan_weibull_survival_model_draws_total <-
              stan_weibull_survival_model_draws %>%
              dplyr::select(.chain, .iteration, .draw, starts_with("yhat_total")) %>%
              tidyr::gather(key = key, value = yhat_total, starts_with("yhat_total")) %>%
              dplyr::select(-key)
            stan_weibull_survival_model_draws_total_exp <-
              stan_weibull_survival_model_draws %>%
              dplyr::select(.chain, .iteration, .draw, starts_with("yhat_factor_total")) %>%
              tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_factor_total")) %>%
              dplyr::select(-key)

            median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=4 * ITER)
            median_os_list_weibull_total = matrix(stan_weibull_survival_model_draws_total$yhat_total, nrow=4 * ITER)
            median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)
            median_os_list_weibull_total = rowMedians(median_os_list_weibull_total,na.rm = T)
            median_os_list_weibull_bias = median_os_list_weibull_total - median_os_list_weibull_total_exp

            median_os_summary_weibull_total = quantile(median_os_list_weibull_total, probs = c(0.025, 0.5, 0.975))
            median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
            median_os_summary_weibull_bias = quantile(median_os_list_weibull_bias, probs = c(0.025, 0.5, 0.975))
            yhat_weibull_sort_total = sort(stan_weibull_survival_model_draws_total$yhat_total)
            yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
            yhat_weibull_sort_average_total = yhat_weibull_sort_total[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total)/length(Data_survival_curve$censor))]
            yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total_exp)/length(Data_survival_curve$censor))]

            Data_survival_curve$factor_neg = yhat_weibull_sort_average_total
            Data_survival_curve$factor_pos = yhat_weibull_sort_average_total_exp
            Data_survival_curve$factor_neg[Data_survival_curve$factor_neg > 10000] = 10000
            Data_survival_curve$factor_pos[Data_survival_curve$factor_pos > 10000] = 10000

            traditional_fit = survfit(Surv(event = censor,
                                           time = time_palliative_final) ~ gene_mut,
                                      data = Data_survival_curve)
            simulation_fit_neg = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_neg)),
                                              time = factor_neg) ~ 1,
                                         data = Data_survival_curve)
            simulation_fit_pos = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_pos)),
                                              time = factor_pos) ~ 1,
                                         data = Data_survival_curve)
            hs[[k]] = ggsurvplot(
              fit = list(traditional_fit = traditional_fit,
                         simulation_fit_neg = simulation_fit_neg,
                         simulation_fit_pos = simulation_fit_pos),
              combine = TRUE,
              data = Data_survival_curve,
              xlab = "Time from CTx induction to final observation (months)",
              ylab = "Survival Probability",
              censor = TRUE,
              surv.scale = "percent",
              font.title = 8,
              font.subtitle = 8,
              font.main = 8,
              font.submain = 8,
              font.caption = 8,
              font.legend = 8,
              surv.median.line = "v",
              palette = "Dark2",
              conf.int = FALSE,
              pval = FALSE,
              risk.table = TRUE,
              risk.table.y.text = FALSE,
              tables.theme = clean_theme(),
              legend = c(0.8,0.8),
              xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
              cumevents = FALSE,
              cumcensor = FALSE,
              break.x.by = 365.25 * 2.5,
              xscale = "d_m",
              #break.x.by = ceiling(max(Data_survival_curve$time_palliative_final) / 500) * 100,
              legend.labs = c(paste(mut_gene, "(-), Unadjusted"),
                              paste(mut_gene, "(+), Unadjusted"),
                              paste(mut_gene, "(-), Adjusted for Delayed Entry"),
                              paste(mut_gene, "(+), Adjusted for Delayed Entry"))
            ) +
              labs(title = paste(sum(traditional_fit[[1]]), " patients ",
                                 "Median OS ",
                                 format_p(summary(traditional_fit)$table[[13]] / 365.25 * 12, digits = 1)," (",
                                 format_p(summary(traditional_fit)$table[[15]] / 365.25 * 12, digits = 1),"-",
                                 format_p(summary(traditional_fit)$table[[17]] / 365.25 * 12, digits = 1),") (mut(-), unadj.)/",
                                 format_p(summary(traditional_fit)$table[[14]] / 365.25 * 12, digits = 1), " (",
                                 format_p(summary(traditional_fit)$table[[16]] / 365.25 * 12, digits = 1),"-",
                                 format_p(summary(traditional_fit)$table[[18]] / 365.25 * 12, digits = 1),") (mut(+), unadj.)/\n",
                                 format_p(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1),
                                 " (", format_p(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1), "-",
                                 format_p(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1), ") (mut(-), adj.)/",
                                 format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                                 format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                                 format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                                 ") (mut(+), adj.) months",
                                 sep=""),
                   subtitle = paste("Survival difference with gene mutation: ",
                                    format_p(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 1), " (",
                                    format_p(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 1), "-",
                                    format_p(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 1),
                                    ") months, median (95% CI)",
                                    sep=""))
            hs[[k]]$table <- hs[[k]]$table + theme(plot.title = element_blank(),
                                                   plot.subtitle = element_blank())
            k = k + 1
            gene_table$positive_median[i] = round(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1)
            gene_table$positive_LL[i] = round(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1)
            gene_table$positive_UL[i] = round(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1)
            gene_table$negative_median[i] = round(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1)
            gene_table$negative_LL[i] = round(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1)
            gene_table$negative_UL[i] = round(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1)
            gene_table$diff_median[i] = -round(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 2)
            gene_table$diff_LL[i] = -round(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 2)
            gene_table$diff_UL[i] = -round(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 2)
          }
        }
        incProgress(1 / length(oncogenic_genes$all_patients))
      }
      for(kk in k:32){
        hs[[kk]] = ggsurvplot_empty()
      }
    })
    incProgress(1 / 13)

    Gene_arrange = gene_table$Gene
    gene_table = gene_table %>% dplyr::mutate(
      Gene = factor(gene_table$Gene, levels = Gene_arrange))
    p_mid <-
      gene_table |>
      ggplot(aes(y = fct_rev(Gene))) +
      theme_classic() +
      geom_point(aes(x=diff_median), shape=15, size=3) +
      geom_linerange(aes(xmin=diff_LL, xmax=diff_UL))  +
      geom_vline(xintercept = 0, linetype="dashed") +
      labs(x="Survival difference by mutation (months)", y="Genes with oncogenic alteration") +
      coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                      xlim=c(-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5,
                             max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)) +
      annotate("text",
               x = (-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5)/2,
               y = length(Gene_arrange) + 1,
               label = "mut=short survival") +
      annotate("text",
               x = (max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)/2,
               y = length(Gene_arrange) + 1,
               label = "mut=long survival") +
      theme(axis.line.y = element_blank(),
            axis.ticks.y= element_blank(),
            axis.text.y= element_blank(),
            axis.title.y= element_blank())

    gene_table <- gene_table %>%
      dplyr::mutate(
        estimate_lab_1 = paste0(format_p(positive_median, digits = 1), " (", format_p(positive_LL, digits = 1), "-", format_p(positive_UL, digits = 1), ")")) %>%
      dplyr::mutate(
        estimate_lab_2 = paste0(format_p(negative_median, digits = 1), " (", format_p(negative_LL, digits = 1), "-", format_p(negative_UL, digits = 1), ")")) %>%
      dplyr::mutate(
        estimate_lab_3 = paste0(format_p(diff_median, digits = 1), " (", format_p(diff_LL, digits = 1), "-", format_p(diff_UL, digits = 1), ")")) %>%
      dplyr::mutate(
        patients = paste0(positive_patients, " (", positive_freq, "%)"))
    gene_table =
      dplyr::bind_rows(
        data.frame(
          Gene = "Gene",
          estimate_lab_1 = "Mut (+) (months)",
          estimate_lab_2 = "Mut (-)",
          estimate_lab_3 = "Survival difference",
          patients = "Patients"
        ), gene_table)
    Gene_arrange = gene_table$Gene
    gene_table = gene_table %>% dplyr::mutate(
      Gene = factor(gene_table$Gene, levels = Gene_arrange))

    p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
    p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                 fontface = ifelse(gene_table$Gene == "Gene", "bold", "bold.italic"))
    p_left <- p_left + geom_text(aes(x = 1.8, label = estimate_lab_1), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_1 == "Mut (+) (months)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 3.3, label = estimate_lab_2), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_2 == "Mut (-)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 4.8, label = estimate_lab_3), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_3 == "Survival difference", "bold", "plain"))
    p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 7.5))

    # right side of plot - pvalues
    p_right <- gene_table  |> ggplot() +
      geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
      theme_void()

    # final plot arrangement
    layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 7.5),
                patchwork::area(t = 1, l = 7.5, b = 30, r = 12.5),
                patchwork::area(t = 0, l = 13.5, b = 30, r = 14))

    h[[1]] = p_left + p_mid + p_right + plot_layout(design = layout)
    OUTPUT_DATA$figure_surv_Bayes_h = h
    OUTPUT_DATA$figure_surv_Bayes_hs = hs
    OUTPUT_DATA$figure_surv_Bayes_gene_table = gene_table
  })
  rm(Bayes_env)
  gc()
})


observeEvent(input$start_Bayes_diagnosis_prediction, {
  req(input$Bayes_diagnosis_no,
      input$Bayes_diagnosis_var,
      OUTPUT_DATA$figure_surv_Bayes_Data_case_target)
  Data_case_target = OUTPUT_DATA$figure_surv_Bayes_Data_case_target
  Bayes_env <- new.env()
  gc()
  with(Bayes_env, {
    if(input$Bayes_diagnosis_var == "diagnosis"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_diagnosis_enroll
    } else if(input$Bayes_diagnosis_var == "1L"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_palliative_enroll
    } else if(input$Bayes_diagnosis_var == "2L"){
      Data_case_target$time_palliative_enroll = Data_case_target$time_2L_enroll
    }
    Data_survival = Data_case_target %>% dplyr::filter(
      !is.na(time_palliative_enroll) &
        !is.na(time_enroll_final) &
        is.finite(time_palliative_enroll) &
        is.finite(time_enroll_final) &
        time_palliative_enroll > 0 &
        time_enroll_final > 0 &
        !is.na(censor) &
        is.finite(censor)
    )
    max_samples = isolate(input$figure_Bayes_max_samples)
    max_samples = ifelse(is.null(max_samples), 2000, max_samples)
    if (nrow(Data_survival) > max_samples){
      sample_indices = sample(nrow(Data_survival),max_samples)
      Data_survival = Data_survival[sample_indices, ]
    } else {
      Data_survival = Data_survival
    }
    if(length(unique(Data_survival$diagnosis)) > 1){
      hb = list()
      hsb = list()
      hsba = list()
      kb = 1
      text_median = NULL
      oncogenic_genes = data.frame(names(table(Data_survival$diagnosis)),
                                   as.integer(table(Data_survival$diagnosis)))
      colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
      oncogenic_genes = oncogenic_genes %>%
        dplyr::arrange(desc(all_patients))
      oncogenic_genes = oncogenic_genes[1:min(nrow(oncogenic_genes), input$Bayes_diagnosis_no),]

      gene_table = data.frame(oncogenic_genes$gene_mutation)
      colnames(gene_table) = c("Gene")
      gene_table$positive_patients = oncogenic_genes$all_patients
      gene_table$positive_freq = format_p(100 * gene_table$positive_patients / length(Data_survival$C.CAT調査結果.基本項目.ハッシュID), digits = 1)
      gene_table$positive_median = 0
      gene_table$positive_LL = 0
      gene_table$positive_UL = 0
      gene_table$negative_median = 0
      gene_table$negative_LL = 0
      gene_table$negative_UL = 0
      gene_table$diff_median = 0
      gene_table$diff_LL = 0
      gene_table$diff_UL = 0
      gene_table$positive_mortal_event = 0
      gene_table$negative_mortal_event = 0
      withProgress(message = "Diagnosis and overall survival", {
        for(i in 1:nrow(oncogenic_genes)){
          ID_mutation = (Data_survival %>% dplyr::filter(diagnosis %in% oncogenic_genes$gene_mutation[i]))$C.CAT調査結果.基本項目.ハッシュID
          Data_survival = Data_survival %>% dplyr::mutate(
            gene_mut = case_when(
              C.CAT調査結果.基本項目.ハッシュID %in% ID_mutation ~ TRUE,
              TRUE ~ FALSE
            )
          )

          mut_gene = oncogenic_genes$gene_mutation[i]
          traditional_fit = survfit(Surv(event = censor,
                                         time = time_palliative_final) ~ gene_mut,
                                    data = Data_survival)
          if(length(Data_survival[,1]) != traditional_fit[[1]][[1]]){
            gene_table$negative_mortal_event[i] = summary(traditional_fit)$table[7]
            gene_table$positive_mortal_event[i] = summary(traditional_fit)$table[8]

            Data_survival_curve = Data_survival
            Data_survival_tmp_pos = Data_survival_curve %>%
              dplyr::filter(gene_mut == TRUE)
            if(sum(Data_survival_tmp_pos$censor) > 0 ){
              fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_pos)
              time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
              survival_rate_10 = summary(fit,times = time_10)[[6]]
              survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_enroll_final)*0.9)[[6]]
            } else {
              time_10 = 90
              survival_rate_10 = 1
              survival_rate_90 = 0
            }
            time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_enroll_final)*0.9)
            Data_survival_tmp = Data_survival_tmp_pos %>% dplyr::filter(
              time_enroll_final >= time_10 &
                time_enroll_final <= time_90)
            Data_survival_tmp$time_enroll_final =
              Data_survival_tmp$time_enroll_final - (time_10 - 1)
            Data_survival_tmp2 = Data_survival_tmp_pos %>% dplyr::filter(
              time_enroll_final > time_90) %>%
              dplyr::arrange(time_enroll_final)

            Data_survival_tmp2$time_enroll_final = time_90 - (time_10 - 1)
            idx <- 1:nrow(Data_survival_tmp2)
            Data_survival_tmp = rbind(Data_survival_tmp, Data_survival_tmp2[idx[idx < max(idx)*(1-survival_rate_10)],])

            Data_survival_tmp_neg = Data_survival_curve %>%
              dplyr::filter(gene_mut == FALSE)
            if(sum(Data_survival_tmp_neg$censor) > 0 ){
              fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_neg)
              time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
              survival_rate_10 = summary(fit,times = time_10)[[6]]
              survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_enroll_final)*0.9)[[6]]
            } else {
              time_10 = 90
              survival_rate_10 = 1
              survival_rate_90 = 0
            }
            time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_enroll_final)*0.9)
            Data_survival_tmp3 = Data_survival_tmp_neg %>% dplyr::filter(
              time_enroll_final >= time_10 &
                time_enroll_final <= time_90)
            Data_survival_tmp3$time_enroll_final =
              Data_survival_tmp3$time_enroll_final - (time_10 - 1)
            Data_survival_tmp4 = Data_survival_tmp_neg %>% dplyr::filter(
              time_enroll_final > time_90) %>%
              dplyr::arrange(time_enroll_final)

            Data_survival_tmp4$time_enroll_final = time_90 - (time_10 - 1)
            idx <- 1:nrow(Data_survival_tmp4)
            Data_survival_tmp3 = rbind(Data_survival_tmp3, Data_survival_tmp4[idx[idx < max(idx)*(1-survival_rate_10)],])
            Data_drug_survival_1 = rbind(Data_survival_tmp, Data_survival_tmp3)

            Data_survival_tmp_pos = Data_survival_curve %>%
              dplyr::filter(gene_mut == TRUE)
            if(sum(Data_survival_tmp_pos$censor) > 0 ){
              fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_pos)
              time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
              survival_rate_10 = summary(fit,times = time_10)[[6]]
              survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)[[6]]
            } else {
              time_10 = 90
              survival_rate_10 = 1
              survival_rate_90 = 0
            }
            time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)
            Data_survival_tmp5 = Data_survival_tmp_pos %>% dplyr::filter(
              time_palliative_enroll >= time_10 &
                time_palliative_enroll <= time_90)
            Data_survival_tmp5$time_palliative_enroll =
              Data_survival_tmp5$time_palliative_enroll - (time_10 - 1)
            Data_survival_tmp6 = Data_survival_tmp_pos %>% dplyr::filter(
              time_palliative_enroll > time_90) %>%
              dplyr::arrange(time_palliative_enroll)
            Data_survival_tmp6$time_palliative_enroll = time_90 - (time_10 - 1)
            idx <- 1:nrow(Data_survival_tmp6)
            Data_survival_tmp5 = rbind(Data_survival_tmp5, Data_survival_tmp6[idx[idx < max(idx)*(1-survival_rate_10)],])

            Data_survival_tmp_neg = Data_survival_curve %>%
              dplyr::filter(gene_mut == FALSE)

            if(sum(Data_survival_tmp_neg$censor) > 0 ){
              fit = survfit(Surv(time_palliative_enroll, enroll_censor) ~ 1, data = Data_survival_tmp_neg)
              time_10 = max(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
              survival_rate_10 = summary(fit,times = time_10)[[6]]
              survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)[[6]]
            } else {
              time_10 = 90
              survival_rate_10 = 1
              survival_rate_90 = 0
            }
            time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)
            Data_survival_tmp7 = Data_survival_tmp_neg %>% dplyr::filter(
              time_palliative_enroll >= time_10 &
                time_palliative_enroll <= time_90)
            Data_survival_tmp7$time_palliative_enroll =
              Data_survival_tmp7$time_palliative_enroll - (time_10 - 1)
            Data_survival_tmp8 = Data_survival_tmp_neg %>% dplyr::filter(
              time_palliative_enroll > time_90) %>%
              dplyr::arrange(time_palliative_enroll)
            Data_survival_tmp8$time_palliative_enroll = time_90 - (time_10 - 1)
            idx <- 1:nrow(Data_survival_tmp8)
            Data_survival_tmp7 = rbind(Data_survival_tmp7, Data_survival_tmp8[idx[idx < max(idx)*(1-survival_rate_10)],])

            Median_entry_pos = median(Data_survival_tmp5$time_palliative_enroll)
            Median_entry_neg = median(Data_survival_tmp7$time_palliative_enroll)
            Data_drug_survival_1_2 = rbind(Data_survival_tmp5, Data_survival_tmp7)

            Data_drug_survival_1 = Data_drug_survival_1 %>%
              dplyr::mutate(delay_1 = case_when(
                time_palliative_enroll <= Median_entry_pos &
                  gene_mut == TRUE ~ "SPT to entry early",
                time_palliative_enroll > Median_entry_pos &
                  gene_mut == TRUE ~ "SPT to entry later",
                time_palliative_enroll <= Median_entry_neg &
                  gene_mut == FALSE ~ "SPT to entry early",
                time_palliative_enroll > Median_entry_neg &
                  gene_mut == FALSE ~ "SPT to entry later"))

            if(sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 == "SPT to entry later"))$censor) >= 1 &
               sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 == "SPT to entry later"))$censor) >= 1 &
               sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 != "SPT to entry later"))$censor) >= 1 &
               sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 != "SPT to entry later"))$censor) >= 1){

              stan_weibull_survival_model_data <-
                list(
                  Median_pos = as.integer(Median_entry_pos),
                  Median_neg = as.integer(Median_entry_neg),
                  Nobs_early_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                  Nobs_late_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                  Ncen_early_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                  Ncen_late_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
                  Nobs_early_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                  Nobs_late_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                  Ncen_early_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                  Ncen_late_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
                  Nexp_pos = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE]),
                  Nexp_neg = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE]),
                  ybef_exp_pos = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE],
                  ybef_exp_neg = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE],
                  yobs_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                  yobs_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                  ycen_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                  ycen_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
                  yobs_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                  yobs_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                  ycen_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
                  ycen_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE]
                )

              if(DOCKER){
                stan_model_compiled <- readRDS(stan_model_factor_save_path_cmdstan)
                stan_weibull_survival_model_fit <- stan_model_compiled$sample(
                  data = stan_weibull_survival_model_data,
                  seed = 1234,
                  chains = 4,
                  parallel_chains = PARALLEL,
                  iter_sampling = ITER,  # (iter-warmup)
                  iter_warmup = as.integer(ITER / 2),adapt_delta=0.99, max_treedepth = 10,
                  thin = 1,
                  refresh = 500
                )
              } else {
                mod <- readRDS(stan_model_factor_save_path_rstan)
                stan_weibull_survival_model_fit <-
                  rstan::sampling(object = mod,
                              data = stan_weibull_survival_model_data,
                              control=list(adapt_delta=0.99, max_treedepth = 10),
                              seed=1234, chains=4, iter=as.integer(ITER * 3 / 2),warmup=as.integer(ITER / 2), thin=1, refresh=500, verbose = FALSE)
              }
              stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_weibull_survival_model_fit)
              stan_weibull_survival_model_draws_total <-
                stan_weibull_survival_model_draws %>%
                dplyr::select(.chain, .iteration, .draw, starts_with("yhat_total")) %>%
                tidyr::gather(key = key, value = yhat_total, starts_with("yhat_total")) %>%
                dplyr::select(-key)
              stan_weibull_survival_model_draws_total_exp <-
                stan_weibull_survival_model_draws %>%
                dplyr::select(.chain, .iteration, .draw, starts_with("yhat_factor_total")) %>%
                tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_factor_total")) %>%
                dplyr::select(-key)

              median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=4 * ITER)
              median_os_list_weibull_total = matrix(stan_weibull_survival_model_draws_total$yhat_total, nrow=4 * ITER)
              median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)
              median_os_list_weibull_total = rowMedians(median_os_list_weibull_total,na.rm = T)
              median_os_list_weibull_bias = median_os_list_weibull_total - median_os_list_weibull_total_exp


              median_os_summary_weibull_total = quantile(median_os_list_weibull_total, probs = c(0.025, 0.5, 0.975))
              median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
              median_os_summary_weibull_bias = quantile(median_os_list_weibull_bias, probs = c(0.025, 0.5, 0.975))
              yhat_weibull_sort_total = sort(stan_weibull_survival_model_draws_total$yhat_total)
              yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
              yhat_weibull_sort_average_total = yhat_weibull_sort_total[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total)/length(Data_survival_curve$censor))]
              yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:length(Data_survival_curve$censor) * as.integer(length(yhat_weibull_sort_total_exp)/length(Data_survival_curve$censor))]

              Data_survival_curve$factor_neg = yhat_weibull_sort_average_total
              Data_survival_curve$factor_pos = yhat_weibull_sort_average_total_exp
              Data_survival_curve$factor_neg[Data_survival_curve$factor_neg > 10000] = 10000
              Data_survival_curve$factor_pos[Data_survival_curve$factor_pos > 10000] = 10000

              traditional_fit = survfit(Surv(event = censor,
                                             time = time_palliative_final) ~ gene_mut,
                                        data = Data_survival_curve)
              simulation_fit_neg = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_neg)),
                                                time = factor_neg) ~ 1,
                                           data = Data_survival_curve)
              simulation_fit_pos = survfit(Surv(event = rep(1,length(Data_survival_curve$factor_pos)),
                                                time = factor_pos) ~ 1,
                                           data = Data_survival_curve)
              hsba[[kb]] = Data_survival_curve
              hsba[[kb]]$diagnosis = mut_gene
              hsb[[kb]] = ggsurvplot(
                fit = list(traditional_fit = traditional_fit,
                           simulation_fit_neg = simulation_fit_neg,
                           simulation_fit_pos = simulation_fit_pos),
                combine = TRUE,
                data = Data_survival_curve,
                xlab = "Time from CTx induction to final observation (months)",
                ylab = "Survival Probability",
                censor = TRUE,
                surv.scale = "percent",
                font.title = 8,
                font.subtitle = 8,
                font.main = 8,
                font.submain = 8,
                font.caption = 8,
                font.legend = 8,
                surv.median.line = "v",
                palette = "Dark2",
                conf.int = FALSE,
                pval = FALSE,
                risk.table = TRUE,
                risk.table.y.text = FALSE,
                tables.theme = clean_theme(),
                legend = c(0.8,0.8),
                xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
                cumevents = FALSE,
                cumcensor = FALSE,
                break.x.by = 365.25 * 2.5,
                xscale = "d_m",
                #break.x.by = ceiling(max(Data_survival_curve$time_palliative_final) / 500) * 100,
                legend.labs = c(paste0("Others, Unadjusted"),
                                paste0(mut_gene, ", Unadjusted"),
                                paste0("Others, Adjusted for Delayed Entry"),
                                paste0(mut_gene, ", Adjusted for Delayed Entry"))
              ) +
                labs(title = paste(sum(traditional_fit[[1]]), " patients ",
                                   "Median OS ",
                                   format_p(summary(traditional_fit)$table[[13]] / 365.25 * 12, digits = 1)," (",
                                   format_p(summary(traditional_fit)$table[[15]] / 365.25 * 12, digits = 1),"-",
                                   format_p(summary(traditional_fit)$table[[17]] / 365.25 * 12, digits = 1),") (others, unadj.)/",
                                   format_p(summary(traditional_fit)$table[[14]] / 365.25 * 12, digits = 1), " (",
                                   format_p(summary(traditional_fit)$table[[16]] / 365.25 * 12, digits = 1),"-",
                                   format_p(summary(traditional_fit)$table[[18]] / 365.25 * 12, digits = 1),") (", mut_gene, ", unadj.)/\n",
                                   format_p(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1),
                                   " (", format_p(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1), "-",
                                   format_p(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1), ") (others, adj.)/",
                                   format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                                   format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                                   format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                                   ") (", mut_gene, ", adj.) months",
                                   sep=""),
                     subtitle = paste("Survival difference with diagnosis: ",
                                      format_p(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 1), " (",
                                      format_p(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 1), "-",
                                      format_p(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 1),
                                      ") months, median (95% CI)",
                                      sep=""))
              hsb[[kb]]$table <- hsb[[kb]]$table + theme(plot.title = element_blank(),
                                                         plot.subtitle = element_blank())
              kb = kb + 1
              text_median = c(text_median, paste0(format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                                                  format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                                                  format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1), ")"))
              gene_table$positive_median[i] = round(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1)
              gene_table$positive_LL[i] = round(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1)
              gene_table$positive_UL[i] = round(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1)
              gene_table$negative_median[i] = round(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1)
              gene_table$negative_LL[i] = round(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1)
              gene_table$negative_UL[i] = round(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1)
              gene_table$diff_median[i] = -round(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 2)
              gene_table$diff_LL[i] = -round(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 2)
              gene_table$diff_UL[i] = -round(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 2)
            }
          }
          incProgress(1 / length(oncogenic_genes$all_patients))
        }
        if(kb>2){
          surv_data = NULL
          for(kb_i in 1:(kb - 1)){
            surv_data = rbind(surv_data, hsba[[kb_i]])
          }
          surv_data = surv_data %>% dplyr::select(factor_pos, diagnosis)
          surv_data$censor = 1
          simulation_fit_diagnosis = survfit(Surv(event = censor,
                                                  time = factor_pos) ~ diagnosis,
                                             data = surv_data)
          survdiff_diagnosis = surv_median(simulation_fit_diagnosis)
          text_median = paste0(paste(text_median, collapse = "/"), " months")
          hsb[[kb]] = ggsurvplot(
            fit = simulation_fit_diagnosis,
            combine = TRUE,
            data = surv_data,
            xlab = "Time from CTx induction to final observation (months)",
            ylab = "Survival Probability",
            censor = TRUE,
            surv.scale = "percent",
            font.title = 8,
            font.subtitle = 8,
            font.main = 8,
            font.submain = 8,
            font.caption = 8,
            font.legend = 8,
            surv.median.line = "v",
            palette = "Paired",
            conf.int = FALSE,
            pval = FALSE,
            risk.table = TRUE,
            risk.table.y.text = FALSE,
            tables.theme = clean_theme(),
            legend = c(0.8,0.8),
            xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
            cumevents = FALSE,
            cumcensor = FALSE,
            break.x.by = 365.25 * 2.5,
            xscale = "d_m"
          ) +
            labs(title = "Median OS by diagnosis, left-truncation bias adjusted",
                 subtitle = text_median
            )
          hsb[[kb]]$table <- hsb[[kb]]$table + theme(plot.title = element_blank(),
                                                     plot.subtitle = element_blank())
          kb = kb + 1
        }
        if(kb>1){
          if(!is.null(OUTPUT_DATA$figure_surv_Bayes_hsb_base)){
            surv_data = NULL
            for(kb_i in 1:(kb - 2)){
              surv_data = rbind(surv_data, hsba[[kb_i]])
            }
            surv_data = surv_data %>% dplyr::select(factor_pos, diagnosis)
            surv_data_tmp = OUTPUT_DATA$figure_surv_Bayes_hsb_base %>% dplyr::select(left_adj, diagnosis)
            colnames(surv_data_tmp) = c("factor_pos", "diagnosis")
            surv_data = rbind(surv_data, surv_data_tmp)
            surv_data$censor = 1
            simulation_fit_diagnosis = survfit(Surv(event = censor,
                                                    time = factor_pos) ~ diagnosis,
                                               data = surv_data)
            survdiff_diagnosis = surv_median(simulation_fit_diagnosis)
            text_median = paste0(OUTPUT_DATA$figure_surv_Bayes_text_median_base, "/", paste(text_median, collapse = "/"), " months")
            hsb[[kb]] = ggsurvplot(
              fit = simulation_fit_diagnosis,
              combine = TRUE,
              data = surv_data,
              xlab = "Time from CTx induction to final observation (months)",
              ylab = "Survival Probability",
              censor = TRUE,
              surv.scale = "percent",
              font.title = 8,
              font.subtitle = 8,
              font.main = 8,
              font.submain = 8,
              font.caption = 8,
              font.legend = 8,
              surv.median.line = "v",
              palette = "Dark2",
              conf.int = FALSE,
              pval = FALSE,
              risk.table = TRUE,
              risk.table.y.text = FALSE,
              tables.theme = clean_theme(),
              legend = c(0.8,0.8),
              xlim = c(0, max(Data_survival_curve$time_palliative_final) * 1.05),
              cumevents = FALSE,
              cumcensor = FALSE,
              break.x.by = 365.25 * 2.5,
              xscale = "d_m"
            ) +
              labs(title = "Median OS by diagnosis, left-truncation bias adjusted",
                   subtitle = text_median
              )
            hsb[[kb]]$table <- hsb[[kb]]$table + theme(plot.title = element_blank(),
                                                       plot.subtitle = element_blank())
            kb = kb + 1
          }
          table_ToT <-
            data.frame(Cohort = as.character(),
                       Duration_months = as.character(),
                       No_total = as.numeric(),
                       No_at_risk = as.numeric(),
                       No_event = as.numeric(),
                       No_censor = as.numeric(),
                       Survival_rate = as.numeric(),
                       Survival_rate_95LL = as.numeric(),
                       Survival_rate_95UL = as.numeric())
          row_no = 1
          for(cohort in 1:length(summary(simulation_fit_diagnosis, times=365.25 / 12, extend=TRUE)$time)){
            for(mths in (0:10)*30){
              dataset = summary(simulation_fit_diagnosis, times=365.25 / 12* mths, extend=TRUE)
              if(length(unique(surv_data$diagnosis)) > 1){
                table_ToT[row_no,"Cohort"] = rownames(dataset$table)[cohort]
              } else{
                table_ToT[row_no,"Cohort"] = unique(surv_data$diagnosis)
              }
              table_ToT[row_no,"Duration_months"] = mths
              table_ToT[row_no,"No_total"] = dataset$n[cohort]
              table_ToT[row_no,"No_at_risk"] = dataset$n.risk[cohort]
              table_ToT[row_no,"No_event"] = dataset$n.event[cohort]
              table_ToT[row_no,"No_censor"] = dataset$n.censor[cohort]
              table_ToT[row_no,"Survival_rate"] = dataset$surv[cohort]
              table_ToT[row_no,"Survival_rate_95LL"] = dataset$lower[cohort]
              table_ToT[row_no,"Survival_rate_95UL"] = dataset$upper[cohort]
              row_no = row_no + 1
            }
          }
          surv_data = NULL
          for(kb_i in 1:(kb - 3)){
            surv_data = rbind(surv_data, hsba[[kb_i]])
          }
          surv_data = surv_data  %>% dplyr::filter(gene_mut == TRUE) %>%
            dplyr::select(time_palliative_final, diagnosis, censor)
          if(!is.null(OUTPUT_DATA$figure_surv_Bayes_hsb_base)){
            surv_data_tmp = OUTPUT_DATA$figure_surv_Bayes_hsb_base %>% dplyr::select(time_palliative_final, diagnosis, censor)
            surv_data = rbind(surv_data, surv_data_tmp)
          }
          simulation_fit_diagnosis = survfit(Surv(event = censor,
                                                  time = time_palliative_final) ~ diagnosis,
                                             data = surv_data)
          for(cohort in 1:length(summary(simulation_fit_diagnosis, times=365.25 / 12, extend=TRUE)$time)){
            for(mths in (0:10)*30){
              dataset = summary(simulation_fit_diagnosis, times=365.25 / 12* mths, extend=TRUE)
              if(length(unique(surv_data$diagnosis)) > 1){
                table_ToT[row_no,"Cohort"] = paste0(rownames(dataset$table)[cohort], "_raw_data")
              } else{
                table_ToT[row_no,"Cohort"] = paste0(unique(surv_data$diagnosis), "_raw_data")
              }
              table_ToT[row_no,"Duration_months"] = mths
              table_ToT[row_no,"No_total"] = dataset$n[cohort]
              table_ToT[row_no,"No_at_risk"] = dataset$n.risk[cohort]
              table_ToT[row_no,"No_event"] = dataset$n.event[cohort]
              table_ToT[row_no,"No_censor"] = dataset$n.censor[cohort]
              table_ToT[row_no,"Survival_rate"] = dataset$surv[cohort]
              table_ToT[row_no,"Survival_rate_95LL"] = dataset$lower[cohort]
              table_ToT[row_no,"Survival_rate_95UL"] = dataset$upper[cohort]
              row_no = row_no + 1
            }
          }
        }
        if(kb < 17){
          for(kkb in kb:16){
            hsb[[kkb]] = ggsurvplot_empty()
          }
        }
        OUTPUT_DATA$figure_surv_Bayes_table_ToT = table_ToT
      })
      Gene_arrange = gene_table$Gene
      gene_table = gene_table %>% dplyr::mutate(
        Gene = factor(gene_table$Gene, levels = Gene_arrange))
      p_mid <-
        gene_table |>
        ggplot(aes(y = fct_rev(Gene))) +
        theme_classic() +
        geom_point(aes(x=diff_median), shape=15, size=3) +
        geom_linerange(aes(xmin=diff_LL, xmax=diff_UL)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Survival difference by diagnosis (months)", y="Diagnosis") +
        coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                        xlim=c(-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5,
                               max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)) +
        annotate("text",
                 x = (-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5)/2,
                 y = length(Gene_arrange) + 1,
                 label = "diagnosis=short survival") +
        annotate("text",
                 x = (max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)/2,
                 y = length(Gene_arrange) + 1,
                 label = "diagnosis=long survival") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y= element_blank(),
              axis.text.y= element_blank(),
              axis.title.y= element_blank())

      gene_table <- gene_table %>%
        dplyr::mutate(
          estimate_lab_1 = paste0(format_p(positive_median, digits = 1), " (", format_p(positive_LL, digits = 1), "-", format_p(positive_UL, digits = 1), ")")) %>%
        dplyr::mutate(
          estimate_lab_2 = paste0(format_p(negative_median, digits = 1), " (", format_p(negative_LL, digits = 1), "-", format_p(negative_UL, digits = 1), ")")) %>%
        dplyr::mutate(
          estimate_lab_3 = paste0(format_p(diff_median, digits = 1), " (", format_p(diff_LL, digits = 1), "-", format_p(diff_UL, digits = 1), ")")) %>%
        dplyr::mutate(
          patients = paste0(positive_patients, " (", positive_freq, "%)"))
      gene_table =
        dplyr::bind_rows(
          data.frame(
            Gene = "Diagnosis",
            estimate_lab_1 = "Survival (months)",
            estimate_lab_2 = "Other diagnosis",
            estimate_lab_3 = "Survival difference",
            patients = "Patients"
          ), gene_table)
      Gene_arrange = gene_table$Gene
      gene_table = gene_table %>% dplyr::mutate(
        Gene = factor(gene_table$Gene, levels = Gene_arrange))

      p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
      p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                   fontface = ifelse(gene_table$Gene == "Diagnosis", "bold", "bold.italic"))
      p_left <- p_left + geom_text(aes(x = 2.4, label = estimate_lab_1), hjust = 0,
                                   fontface = ifelse(gene_table$estimate_lab_1 == "Survival (months)", "bold", "plain"))
      p_left <- p_left + geom_text(aes(x = 3.9, label = estimate_lab_2), hjust = 0,
                                   fontface = ifelse(gene_table$estimate_lab_2 == "Other diagnosis", "bold", "plain"))
      p_left <- p_left + geom_text(aes(x = 5.4, label = estimate_lab_3), hjust = 0,
                                   fontface = ifelse(gene_table$estimate_lab_3 == "Survival difference", "bold", "plain"))
      p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 9.0))

      # right side of plot - pvalues
      p_right <- gene_table  |> ggplot() +
        geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                  fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
        theme_void()

      # final plot arrangement
      layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 9.0),
                  patchwork::area(t = 1, l = 7.5, b = 30, r = 14.0),
                  patchwork::area(t = 0, l = 13.5, b = 30, r = 15.5))

      hb[[1]] = p_left + p_mid + p_right + plot_layout(design = layout)
      OUTPUT_DATA$figure_surv_Bayes_gene_table_hb = gene_table
      OUTPUT_DATA$figure_surv_Bayes_hb = hb
      OUTPUT_DATA$figure_surv_Bayes_hsb = hsb

    }
  })
  rm(Bayes_env)
  gc()
})


output$figure_survival_CTx_5_table = DT::renderDataTable(server = FALSE,{
  req(OUTPUT_DATA$figure_surv_Bayes_table_ToT)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_Bayes_table_ToT, "FELIS downloaded survival curve data in Diagnosis and survival tab")
})

output$figure_survival_CTx_4_data = DT::renderDataTable(server = FALSE,{
  req(OUTPUT_DATA$figure_surv_Bayes_gene_table_hb)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_Bayes_gene_table_hb, "FELIS downloaded forest plot data in Diagnosis and survival tab")
})

output$figure_survival_CTx_2_data = DT::renderDataTable(server = FALSE,{
  req(OUTPUT_DATA$figure_surv_Bayes_gene_table)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_Bayes_gene_table, "FELIS downloaded raw data in Genetic variants and survival tab")
})

output$figure_survival_CTx_1 = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_g)
  arrange_ggsurvplots(OUTPUT_DATA$figure_surv_Bayes_g,print=TRUE,ncol=1,nrow=3)
})

output$figure_survival_CTx_compare = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_compare)
  OUTPUT_DATA$figure_surv_Bayes_compare
})

output$figure_survival_CTx_4 = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_hb)
  OUTPUT_DATA$figure_surv_Bayes_hb[[1]]
})

output$figure_survival_CTx_3 = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_hs)
  arrange_ggsurvplots(OUTPUT_DATA$figure_surv_Bayes_hs,print=TRUE,ncol=4,nrow=8,surv.plot.height = 0.75,risk.table.height = 0.25)
})


output$figure_survival_CTx_5 = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_hsb)
  arrange_ggsurvplots(OUTPUT_DATA$figure_surv_Bayes_hsb,print=TRUE,ncol=4,nrow=4,surv.plot.height = 0.75,risk.table.height = 0.25)
})

output$figure_survival_CTx_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_Bayes_h)
  OUTPUT_DATA$figure_surv_Bayes_h[[1]]
})
