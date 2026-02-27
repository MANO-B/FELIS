survival_CTx_analysis2_logic_control <- function() {
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
                      HER2_IHC,
                      MSI_PCR,
                      MMR_IHC,
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
      if(input$HER2 == "No"){
        Data_case_target = Data_case_target %>%
          dplyr::select(-HER2_IHC)
      } else {
        Data_case_target = Data_case_target %>%
          dplyr::filter(
            HER2_IHC != "Unknown"
          )
      }
      if(input$MSI == "No"){
        Data_case_target = Data_case_target %>%
          dplyr::select(-MSI_PCR)
      } else {
        Data_case_target = Data_case_target %>%
          dplyr::filter(
            MSI_PCR != "Unknown"
          )
      }
      if(input$MMR == "No"){
        Data_case_target = Data_case_target %>%
          dplyr::select(-MMR_IHC)
      } else {
        Data_case_target = Data_case_target %>%
          dplyr::filter(
            MMR_IHC != "Unknown"
          )
      }
      Data_case_target$Cancers = Data_case_target$症例.基本情報.がん種.OncoTree.
      incProgress(1 / 13)
      Data_drug = Data_drug_raw()
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control = Data_drug

      Data_case_target = Data_case_target %>%
        dplyr::distinct(.keep_all = TRUE, C.CAT調査結果.基本項目.ハッシュID)
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
      if(nrow(Data_case_target)>0){
        Data_report_TMB = Data_MAF %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
          dplyr::select(Tumor_Sample_Barcode, TMB) %>%
          dplyr::distinct(Tumor_Sample_Barcode, .keep_all=T)
        colnames(Data_report_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB")
        Data_case_target = left_join(Data_case_target, Data_report_TMB,
                                     by = "C.CAT調査結果.基本項目.ハッシュID")

        Data_MAF_target = Data_MAF %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
        if(input$patho == "Only pathogenic muts"){
          Data_MAF_target = Data_MAF_target %>%
            dplyr::filter(Evidence_level == "F")
        }
        Data_MAF_target = Data_MAF_target %>%
          dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, .keep_all = T)
        Gene_list = unique(names(sort(table(Data_MAF_target$Hugo_Symbol),
                                      decreasing = T)))
        if(length(Data_case_target$censor[is.na(Data_case_target$censor)])> 0){
          Data_case_target$censor[is.na(Data_case_target$censor)] = 0
        }
        incProgress(1 / 13)

        Data_case_target$time_pre = Data_case_target$time_diagnosis_enroll
        Data_case_target$time_all = Data_case_target$time_diagnosis_final
        Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)
        adjustment = TRUE
        OUTPUT_DATA$figure_surv_CTx_adjustment_control = adjustment

        Data_case_target = Data_case_target %>% dplyr::filter(
          !is.na(time_enroll_final) &
            is.finite(time_enroll_final) &
            time_enroll_final > 0 &
            !is.na(time_pre) &
            is.finite(time_pre) &
            !is.na(censor) &
            is.finite(censor)
        )
        Data_survival = Data_case_target
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_survival$C.CAT調査結果.基本項目.ハッシュID)
        OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control = Data_MAF_target
        incProgress(1 / 13)

        Data_survival_interactive = Data_survival
        OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control = Data_survival_interactive
        candidate_genes = sort(unique(c(Data_MAF_target$Hugo_Symbol,
                                        paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                                        paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                                        paste0(input$special_gene, "_NOS"))))
        candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
        candidate_genes = candidate_genes[!is.na(candidate_genes)]
        Top_gene = unique(names(sort(table(Data_MAF_target$Hugo_Symbol), decreasing = T)))
        Top_gene = Top_gene[Top_gene %in% candidate_genes]
        candidate_drugs = sort(unique(c(Data_drug$Drug)))
        OUTPUT_DATA$figure_surv_interactive_Top_gene_control = Top_gene
        OUTPUT_DATA$figure_surv_interactive_candidate_genes_control = candidate_genes
        OUTPUT_DATA$figure_surv_interactive_candidate_drugs_control = candidate_drugs[!is.na(candidate_drugs)]
        OUTPUT_DATA$figure_surv_interactive_candidate_Age_control = sort(unique(Data_survival_interactive$YoungOld))
        OUTPUT_DATA$figure_surv_interactive_candidate_Sex_control = sort(unique(Data_survival_interactive$症例.基本情報.性別.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_Histology_control = sort(unique(Data_survival_interactive$Cancers))
      }
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}


# =========================================================================
# GLOBAL HELPERS (Tamura & Ikegami Model Ver 2.3.2 実データ用)
# エラー完全回避版：モデル行列を安全に生成するサンプリング関数
# =========================================================================
gen_sim_times_robust <- function(fit, newdata, dist_type = c("weibull", "llogis")) {
  dist_type <- match.arg(dist_type)
  if (is.null(fit) || nrow(newdata) == 0) return(rep(NA_real_, nrow(newdata)))

  shape <- fit$res["shape", "est"]
  base_scale <- fit$res["scale", "est"]

  # 【修正】Symbolエラーを防ぐため、埋め込まれたフォーミュラを直接使用
  form <- fit$my_formula
  if (is.null(form)) {
    form <- tryCatch(eval(fit$call$formula, envir = parent.frame()), error = function(e) fit$call$formula)
  }

  # 共変量がない場合
  if (length(fit$coefficients) == 2) {
    scale_vec <- rep(base_scale, nrow(newdata))
  } else {
    # 応答変数を削除した右辺だけのフォーミュラ (~ X1 + X2) を安全に生成
    form_str <- as.character(form)
    rhs_str <- if (length(form_str) == 3) form_str[3] else form_str[2]
    form_right <- as.formula(paste("~", rhs_str))

    # 欠損値でエラーが出ないよう安全にモデル行列を作成
    mf <- model.frame(form_right, data = newdata, xlev = fit$xlevels, na.action = na.pass)
    mm <- model.matrix(form_right, mf, contrasts.arg = fit$contrasts)

    cov_names <- colnames(mm)[-1] # Interceptを除外
    cov_coefs <- fit$coefficients[cov_names]
    cov_coefs[is.na(cov_coefs)] <- 0 # 安全装置

    # 線形予測子（eta）から個別のscaleを計算
    eta <- as.numeric(mm[, -1, drop = FALSE] %*% cov_coefs)
    scale_vec <- base_scale * exp(eta)
  }

  u <- pmax(pmin(runif(nrow(newdata)), 0.9999), 0.0001)

  if (dist_type == "weibull") {
    return(scale_vec * (-log(u))^(1 / shape))
  } else {
    return(scale_vec * ((1 - u) / u)^(1 / shape))
  }
}

# =========================================================================
# STEP 1: Calibration Weighting
# =========================================================================
calculate_calibrated_iptw <- function(data, ref_surv_list) {
  data$iptw <- 1.0
  t_points_days <- (1:5) * 365.25

  for (ag in unique(data$age_class)) {
    if (!(ag %in% names(ref_surv_list))) next

    S_macro_target <- pmax(pmin(ref_surv_list[[ag]][1:5] / 100, 0.999), 0.001)
    ag_data_idx <- which(data$age_class == ag)
    ag_data <- data[ag_data_idx, , drop = FALSE]

    if (nrow(ag_data) < 5) next

    log_t1 <- log(pmax(ag_data$time_pre / 365.25, 1e-6))
    log_t1_sq <- log_t1^2

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t1 + theta[2] * log_t1_sq)
      w <- pmax(w, 1e-4)
      w <- w / mean(w)

      km <- tryCatch(survfit(Surv(time_all, censor) ~ 1, data = ag_data, weights = w), error = function(e) NULL)
      if (is.null(km)) return(1e6)

      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      S_est[is.na(S_est)] <- min(S_est, na.rm = TRUE)

      sum((S_est - S_macro_target)^2) + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(optim(c(0, 0), obj_func, method = "BFGS"), error = function(e) list(par = c(0, 0)))
    w_opt <- pmax(exp(opt$par[1] * log_t1 + opt$par[2] * log_t1_sq), 1e-4)

    lower_bound <- quantile(w_opt, 0.01, na.rm = TRUE)
    upper_bound <- quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lower_bound, pmin(w_opt, upper_bound))

    data$iptw[ag_data_idx] <- w_opt / mean(w_opt)
  }
  return(data)
}

# =========================================================================
# 第1部：生存曲線描画ロジック (figure_survival_CTx_interactive_1_control)
# =========================================================================
output$figure_survival_CTx_interactive_1_control = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control,
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control)

  lapply(1:2, function(i) {
    prefix <- paste0("gene_survival_interactive_", i, "_")
    list(
      input[[paste0(prefix, "P_1_control")]], input[[paste0(prefix, "P_2_control")]],
      input[[paste0(prefix, "W_control")]], input[[paste0(prefix, "A_control")]],
      input[[paste0(prefix, "S_control")]], input[[paste0(prefix, "H_control")]],
      input[[paste0(prefix, "D_control")]]
    )
  })

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  Data_drug = OUTPUT_DATA$figure_surv_CTx_Data_drug_control

  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = as.numeric(症例.基本情報.年齢),
      age_class = dplyr::case_when(
        age_num < 40 ~ "40未満", age_num < 50 ~ "40代", age_num < 60 ~ "50代",
        age_num < 70 ~ "60代", age_num < 80 ~ "70代", age_num >= 80 ~ "80以上",
        TRUE ~ "全年齢"
      )
    )

  ref_surv_list <- list()
  age_groups <- c("40未満", "40代", "50代", "60代", "70代", "80以上")

  if (!is.null(input$survival_data_source)) {
    if (input$survival_data_source == "registry") {
      if(exists("Data_age_survival_5_year") && input$registry_cancer_type %in% names(Data_age_survival_5_year)) {
        cancer_data <- Data_age_survival_5_year[[input$registry_cancer_type]]
        fallback_surv <- cancer_data[["全年齢"]]
        for (ag in age_groups) {
          if (!is.null(cancer_data[[ag]]) && length(cancer_data[[ag]]) == 5) {
            ref_surv_list[[ag]] <- as.numeric(cancer_data[[ag]])
          } else if (!is.null(fallback_surv) && length(fallback_surv) == 5) {
            ref_surv_list[[ag]] <- as.numeric(fallback_surv)
          }
        }
      }
    } else if (input$survival_data_source == "manual_all") {
      vals <- c(input$surv_all_1y, input$surv_all_2y, input$surv_all_3y, input$surv_all_4y, input$surv_all_5y)
      for (ag in age_groups) ref_surv_list[[ag]] <- vals
    } else if (input$survival_data_source == "manual_age") {
      ref_surv_list[["40未満"]] <- c(input$surv_u40_1y, input$surv_u40_2y, input$surv_u40_3y, input$surv_u40_4y, input$surv_u40_5y)
      ref_surv_list[["40代"]] <- c(input$surv_40s_1y, input$surv_40s_2y, input$surv_40s_3y, input$surv_40s_4y, input$surv_40s_5y)
      ref_surv_list[["50代"]] <- c(input$surv_50s_1y, input$surv_50s_2y, input$surv_50s_3y, input$surv_50s_4y, input$surv_50s_5y)
      ref_surv_list[["60代"]] <- c(input$surv_60s_1y, input$surv_60s_2y, input$surv_60s_3y, input$surv_60s_4y, input$surv_60s_5y)
      ref_surv_list[["70代"]] <- c(input$surv_70s_1y, input$surv_70s_2y, input$surv_70s_3y, input$surv_70s_4y, input$surv_70s_5y)
      ref_surv_list[["80以上"]] <- c(input$surv_80s_1y, input$surv_80s_2y, input$surv_80s_3y, input$surv_80s_4y, input$surv_80s_5y)
    }
  }

  if (length(ref_surv_list) == 0) {
    for (ag in age_groups) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
  } else {
    for (ag in age_groups) { if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- ref_surv_list[[1]] }
  }

  if (length(ref_surv_list) > 0) {
    Data_survival_interactive <- calculate_calibrated_iptw(Data_survival_interactive, ref_surv_list)
  } else {
    Data_survival_interactive$iptw <- 1.0
  }
  Data_survival_interactive$iptw <- ifelse(is.na(Data_survival_interactive$iptw) | Data_survival_interactive$iptw <= 0, 1.0, Data_survival_interactive$iptw)

  extract_group_ids <- function(group_num) {
    IDs <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)
    input_prefix <- paste0("gene_survival_interactive_", group_num, "_")
    p1_input <- input[[paste0(input_prefix, "P_1_control")]]
    if(!all(is.null(p1_input))) IDs <- intersect(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
    p2_input <- input[[paste0(input_prefix, "P_2_control")]]
    if(!all(is.null(p2_input))) IDs <- intersect(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
    w_input <- input[[paste0(input_prefix, "W_control")]]
    if(!all(is.null(w_input))) IDs <- setdiff(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% w_input))$Tumor_Sample_Barcode)
    clinical_filters <- list(A_control = "YoungOld", S_control = "症例.基本情報.性別.名称.", H_control = "Cancers")
    for(filter_key in names(clinical_filters)) {
      filter_input <- input[[paste0(input_prefix, filter_key)]]
      if(!all(is.null(filter_input))) {
        column_name <- clinical_filters[[filter_key]]
        filter_expr <- paste0(column_name, " %in% filter_input")
        IDs <- intersect(IDs, (Data_survival_interactive %>% dplyr::filter(!!rlang::parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
      }
    }
    d_input <- input[[paste0(input_prefix, "D_control")]]
    if(!all(is.null(d_input))) IDs <- intersect(IDs, (Data_drug %>% dplyr::filter(Drug %in% d_input))$ID)
    nd_input <- input[[paste0(input_prefix, "ND_control")]]
    if(!all(is.null(nd_input))) IDs <- setdiff(IDs, (Data_drug %>% dplyr::filter(Drug %in% nd_input))$ID)
    return(IDs)
  }

  ID_1 <- extract_group_ids(1)
  ID_2 <- extract_group_ids(2)

  Data_survival_1 <- Data_survival_interactive %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>% dplyr::mutate(Group = "1")
  Data_survival_2 <- Data_survival_interactive %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>% dplyr::mutate(Group = "2")
  Data_EP0 <- Data_survival_interactive %>% dplyr::filter(EP_treat == 0) %>% dplyr::mutate(Group = "EP0")

  Data_model <- rbind(Data_EP0, Data_survival_1, Data_survival_2) %>%
    dplyr::mutate(
      time_t2 = pmax(time_all - time_pre, 0.1),
      time_pre_event = 1,
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.factor(Cancers)
    ) %>%
    dplyr::filter(time_pre > 0, time_all > time_pre)

  shiny::validate(shiny::need(nrow(Data_model) > 0, "No valid data to fit the model."))
  Data_model$Group <- factor(Data_model$Group, levels = c("EP0", "1", "2"))

  library(splines)
  Data_model$logT1_scale <- log(pmax(Data_model$time_pre / 365.25, 1e-6))
  ns_obj <- ns(Data_model$logT1_scale, df = 3)
  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for(j in seq_len(ncol(ns_mat))) Data_model[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  potential_covariates <- c("age_num", "Sex", "Cancers")
  valid_covariates <- c("Group")
  for (cov in potential_covariates) {
    if (cov %in% colnames(Data_model) && length(unique(na.omit(Data_model[[cov]]))) > 1) {
      valid_covariates <- c(valid_covariates, cov)
    }
  }

  req(requireNamespace("flexsurv", quietly = TRUE))

  # 【修正】フォーミュラをモデル内に強制保存 (my_formula)
  form_t1 <- as.formula(paste("Surv(time_pre, time_pre_event) ~", paste(valid_covariates, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurv::flexsurvreg(form_t1, data = Data_model, weights = iptw, dist = "weibull"), error = function(e) NULL)
  if (!is.null(fit_t1)) fit_t1$my_formula <- form_t1

  form_t2 <- as.formula(paste("Surv(time_t2, censor) ~", paste(c(valid_covariates, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurv::flexsurvreg(form_t2, data = Data_model, weights = iptw, dist = "llogis"), error = function(e) NULL)
  if (!is.null(fit_t2)) fit_t2$my_formula <- form_t2

  shiny::validate(shiny::need(!is.null(fit_t1) && !is.null(fit_t2), "Model fitting failed. The dataset might be too small or sparse."))

  set.seed(2024)
  idx_pseudo <- sample(seq_len(nrow(Data_model)), size = 5000, replace = TRUE, prob = Data_model$iptw)
  pseudo_base <- Data_model[idx_pseudo, ]

  simulated_data <- list()
  cap_days <- 365.25 * 10

  for (g in c("1", "2")) {
    pseudo_g <- pseudo_base
    pseudo_g$Group <- factor(g, levels = levels(Data_model$Group))

    sim_T1 <- gen_sim_times_robust(fit_t1, pseudo_g, dist_type = "weibull")

    pseudo_g$logT1_scale <- log(pmax(sim_T1 / 365.25, 1e-6))
    ns_sim <- tryCatch(predict(ns_obj, pseudo_g$logT1_scale), error = function(e) matrix(0, nrow(pseudo_g), attr(ns_obj, "df")))
    colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
    for(j in seq_len(ncol(ns_sim))) pseudo_g[[colnames(ns_sim)[j]]] <- ns_sim[, j]

    sim_T2 <- gen_sim_times_robust(fit_t2, pseudo_g, dist_type = "llogis")
    sim_OS <- sim_T1 + sim_T2

    final_censor <- ifelse(sim_OS > cap_days, 0, 1)
    final_os <- pmin(sim_OS, cap_days)
    final_t1 <- pmin(sim_T1, cap_days - 0.1)

    simulated_data[[length(simulated_data) + 1]] <- data.frame(
      C.CAT調査結果.基本項目.ハッシュID = pseudo_g$C.CAT調査結果.基本項目.ハッシュID,
      time_pre = final_t1,
      time_all = final_os,
      censor = final_censor,
      Group = g
    )
  }

  Data_survival_simulated <- do.call(rbind, simulated_data) %>% dplyr::filter(is.finite(time_all))
  shiny::validate(shiny::need(nrow(Data_survival_simulated) > 0, "Simulation yielded no data."))

  survival_compare_and_plot_CTx(
    data = Data_survival_simulated,
    time_var1 = "time_pre",
    time_var2 = "time_all",
    status_var = "censor",
    group_var = "Group",
    plot_title = "Standardized OS (Ver 2.3.2: Calibrated G-computation)",
    adjustment = FALSE,
    color_var_surv_CTx_1 = "diagnosis",
    weights_var = NULL
  )
})

# =========================================================================
# 第2部：フォレストプロット (forest_plot_multivariate)
# =========================================================================
output$forest_plot_multivariate = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control)

  shiny::validate(shiny::need(requireNamespace("patchwork", quietly = TRUE), "Please install the 'patchwork' package."))
  shiny::validate(shiny::need(requireNamespace("splines", quietly = TRUE), "Please install the 'splines' package."))

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  selected_genes <- input$gene_survival_interactive_1_P_1_control_forest

  ref_surv_list <- list()
  age_groups <- c("40未満", "40代", "50代", "60代", "70代", "80以上")

  if (!is.null(input$survival_data_source)) {
    if (input$survival_data_source == "registry") {
      if(exists("Data_age_survival_5_year") && input$registry_cancer_type %in% names(Data_age_survival_5_year)) {
        cancer_data <- Data_age_survival_5_year[[input$registry_cancer_type]]
        fallback_surv <- cancer_data[["全年齢"]]
        for (ag in age_groups) {
          if (!is.null(cancer_data[[ag]]) && length(cancer_data[[ag]]) == 5) {
            ref_surv_list[[ag]] <- as.numeric(cancer_data[[ag]])
          } else if (!is.null(fallback_surv) && length(fallback_surv) == 5) {
            ref_surv_list[[ag]] <- as.numeric(fallback_surv)
          }
        }
      }
    } else if (input$survival_data_source == "manual_all" || input$survival_data_source == "manual_age") {
      for (ag in age_groups) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
    }
  }

  if (length(ref_surv_list) == 0) {
    for (ag in age_groups) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
  } else {
    for (ag in age_groups) { if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- ref_surv_list[[1]] }
  }

  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = as.numeric(症例.基本情報.年齢),
      age_class = dplyr::case_when(
        age_num < 40 ~ "40未満", age_num < 50 ~ "40代", age_num < 60 ~ "50代",
        age_num < 70 ~ "60代", age_num < 80 ~ "70代", age_num >= 80 ~ "80以上",
        TRUE ~ "全年齢"
      )
    )

  if (length(ref_surv_list) > 0) {
    Data_survival_interactive <- calculate_calibrated_iptw(Data_survival_interactive, ref_surv_list)
  } else {
    Data_survival_interactive$iptw <- 1.0
  }
  Data_survival_interactive$iptw <- ifelse(is.na(Data_survival_interactive$iptw) | Data_survival_interactive$iptw <= 0, 1.0, Data_survival_interactive$iptw)

  Data_forest <- Data_survival_interactive %>%
    dplyr::filter(time_pre > 0, time_all > time_pre) %>%
    dplyr::mutate(
      time_t2 = pmax(time_all - time_pre, 0.1),
      time_pre_event = 1,
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.character(Cancers)
    )

  shiny::validate(shiny::need(nrow(Data_forest) > 50, "Not enough valid cases for multivariate analysis."))
  total_cohort_N <- nrow(Data_forest)

  if(length(levels(Data_forest$Sex)) > 1) {
    top_sex <- names(sort(table(Data_forest$Sex), decreasing = TRUE))[1]
    Data_forest$Sex <- relevel(Data_forest$Sex, ref = top_sex)
  }

  cancer_counts <- table(Data_forest$Cancers)
  rare_cancers <- names(cancer_counts)[cancer_counts < 50]
  Data_forest$Cancers[Data_forest$Cancers %in% rare_cancers] <- "Others"
  Data_forest$Cancers <- as.factor(Data_forest$Cancers)
  top_cancer <- names(sort(table(Data_forest$Cancers), decreasing = TRUE))[1]
  Data_forest$Cancers <- relevel(Data_forest$Cancers, ref = top_cancer)

  valid_gene_covariates <- c()
  if (length(selected_genes) > 0) {
    for (gene in selected_genes) {
      mutated_IDs <- (Data_MAF_target %>% dplyr::filter(Hugo_Symbol == gene))$Tumor_Sample_Barcode
      col_name <- paste0("Gene_", make.names(gene))
      Data_forest[[col_name]] <- factor(ifelse(Data_forest$C.CAT調査結果.基本項目.ハッシュID %in% mutated_IDs, "Mutated", "WT"), levels = c("WT", "Mutated"))
      if (length(unique(Data_forest[[col_name]])) > 1) valid_gene_covariates <- c(valid_gene_covariates, col_name)
    }
  }

  base_covariates <- c("age_num", "Sex", "Cancers")
  all_covariates <- c()
  for(cov in base_covariates) {
    if(length(unique(na.omit(Data_forest[[cov]]))) > 1) all_covariates <- c(all_covariates, cov)
  }
  all_covariates <- c(all_covariates, valid_gene_covariates)
  shiny::validate(shiny::need(length(all_covariates) > 0, "No valid covariates found."))

  library(splines)
  Data_forest$logT1_scale <- log(pmax(Data_forest$time_pre / 365.25, 1e-6))
  ns_obj <- ns(Data_forest$logT1_scale, df = 3)
  ns_mat <- as.matrix(ns_obj)
  colnames(ns_mat) <- paste0("ns", seq_len(ncol(ns_mat)))
  for(j in seq_len(ncol(ns_mat))) Data_forest[[colnames(ns_mat)[j]]] <- ns_mat[, j]

  req(requireNamespace("flexsurv", quietly = TRUE))

  # 【修正】モデルの式を安全に保持させる
  form_t1 <- as.formula(paste("Surv(time_pre, time_pre_event) ~", paste(all_covariates, collapse = " + ")))
  fit_t1 <- tryCatch(flexsurv::flexsurvreg(form_t1, data = Data_forest, weights = iptw, dist = "weibull"), error = function(e) NULL)
  if (!is.null(fit_t1)) fit_t1$my_formula <- form_t1

  form_t2 <- as.formula(paste("Surv(time_t2, censor) ~", paste(c(all_covariates, colnames(ns_mat)), collapse = " + ")))
  fit_t2 <- tryCatch(flexsurv::flexsurvreg(form_t2, data = Data_forest, weights = iptw, dist = "llogis"), error = function(e) NULL)
  if (!is.null(fit_t2)) fit_t2$my_formula <- form_t2

  shiny::validate(shiny::need(!is.null(fit_t1) && !is.null(fit_t2), "T1/T2 modeling failed. Matrix may be singular."))

  set.seed(2024)
  M_pseudo <- 10000
  idx_pseudo <- sample(seq_len(nrow(Data_forest)), size = M_pseudo, replace = TRUE, prob = Data_forest$iptw)
  pseudo_forest <- Data_forest[idx_pseudo, ]

  pseudo_forest$sim_T1 <- gen_sim_times_robust(fit_t1, pseudo_forest, dist_type = "weibull")

  pseudo_forest$logT1_sim_scale <- log(pmax(pseudo_forest$sim_T1 / 365.25, 1e-6))
  ns_sim <- tryCatch(predict(ns_obj, pseudo_forest$logT1_sim_scale), error = function(e) matrix(0, nrow(pseudo_forest), attr(ns_obj, "df")))
  colnames(ns_sim) <- paste0("ns", seq_len(ncol(ns_sim)))
  for(j in seq_len(ncol(ns_sim))) pseudo_forest[[colnames(ns_sim)[j]]] <- ns_sim[, j]

  pseudo_forest$sim_T2 <- gen_sim_times_robust(fit_t2, pseudo_forest, dist_type = "llogis")
  pseudo_forest$sim_OS <- pseudo_forest$sim_T1 + pseudo_forest$sim_T2
  pseudo_forest$sim_Event <- 1

  form_final <- as.formula(paste("Surv(sim_OS, sim_Event) ~", paste(all_covariates, collapse = " + ")))
  fit_forest <- tryCatch(flexsurv::flexsurvreg(form_final, data = pseudo_forest, dist = "llogis"), error = function(e) NULL)

  shiny::validate(shiny::need(!is.null(fit_forest), "Final synchronized multivariate model failed to converge."))

  res <- fit_forest$res
  cov_names <- rownames(res)[!rownames(res) %in% c("shape", "scale")]

  plot_data <- data.frame(
    Variable_Raw = cov_names,
    Estimate = exp(res[cov_names, "est"]),
    Lower = exp(res[cov_names, "L95%"]),
    Upper = exp(res[cov_names, "U95%"]),
    Category = NA, Label = NA, Positive_N = NA, Total_N = total_cohort_N, Text_CI = NA
  )

  for (i in 1:nrow(plot_data)) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^age_num", var)) {
      plot_data$Estimate[i] <- exp(res[var, "est"] * 10)
      plot_data$Lower[i] <- exp(res[var, "L95%"] * 10)
      plot_data$Upper[i] <- exp(res[var, "U95%"] * 10)
    }
  }

  plot_data$Text_CI <- sprintf("%.2f (%.2f-%.2f)", plot_data$Estimate, plot_data$Lower, plot_data$Upper)

  for (i in 1:nrow(plot_data)) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^age_num", var)) {
      plot_data$Label[i] <- "Age (per +10 years)"
      plot_data$Category[i] <- "Demographics"
      plot_data$Positive_N[i] <- NA
    } else if (grepl("^Sex", var)) {
      level_name <- sub("^Sex", "", var)
      plot_data$Label[i] <- paste0("Sex: ", level_name, " (vs ", levels(Data_forest$Sex)[1], ")")
      plot_data$Category[i] <- "Demographics"
      plot_data$Positive_N[i] <- sum(Data_forest$Sex == level_name, na.rm = TRUE)
    } else if (grepl("^Cancers", var)) {
      cancer_type <- sub("^Cancers", "", var)
      plot_data$Label[i] <- paste0("Histology: ", cancer_type, " (vs ", levels(Data_forest$Cancers)[1], ")")
      plot_data$Category[i] <- "Histology"
      plot_data$Positive_N[i] <- sum(Data_forest$Cancers == cancer_type, na.rm = TRUE)
    } else if (grepl("^Gene_", var)) {
      gene_name <- sub("^Gene_", "", var)
      gene_name <- sub("Mutated$", "", gene_name)
      col_name <- paste0("Gene_", gene_name)
      plot_data$Label[i] <- paste0("Gene: ", gene_name, " (Mutated vs WT)")
      plot_data$Category[i] <- "Genomics"
      plot_data$Positive_N[i] <- sum(Data_forest[[col_name]] == "Mutated", na.rm = TRUE)
    } else {
      plot_data$Label[i] <- var
      plot_data$Category[i] <- "Others"
    }
  }

  plot_data$Category <- factor(plot_data$Category, levels = c("Genomics", "Demographics", "Histology", "Others"))
  plot_data <- plot_data[order(plot_data$Category, plot_data$Estimate), ]
  plot_data$Label <- factor(plot_data$Label, levels = plot_data$Label)
  plot_data$Text_PosN <- ifelse(is.na(plot_data$Positive_N), "-", as.character(plot_data$Positive_N))

  library(patchwork)

  p_left <- ggplot(plot_data, aes(x = Estimate, y = Label, color = Category)) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.3, linewidth = 0.8) +
    geom_point(size = 4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7f8c8d", linewidth = 1) +
    scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_color_manual(values = c("Genomics" = "#e74c3c", "Demographics" = "#2980b9", "Histology" = "#27ae60", "Others" = "#8e44ad")) +
    theme_minimal(base_size = 14) +
    labs(x = "Time Ratio (Log Scale, 95% CI)", y = "") +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.position = "none"
    )

  p_right <- ggplot(plot_data, aes(y = Label)) +
    geom_text(aes(x = 0, label = Total_N), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 1, label = Text_PosN), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 2, label = Text_CI), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    scale_x_continuous(limits = c(-0.5, 2.5), breaks = c(0, 1, 2), labels = c("Total N\n(Real)", "Positive N\n(Real)", "TR\n(95% CI)")) +
    theme_minimal(base_size = 14) +
    labs(x = "", y = "") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(face = "bold", color = "#2c3e50"),
      axis.ticks = element_blank()
    )

  final_plot <- p_left + p_right + plot_layout(widths = c(2, 1.2)) +
    plot_annotation(
      title = "Independent Prognostic Impact (Standardized AFT Model)",
      subtitle = "Computed via Calibrated G-computation (Ver 2.3.2) correcting Immortal Time & Dependent Truncation\nTR > 1: Prolonged Survival | TR < 1: Shortened Survival",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "#34495e")
      )
    )

  return(final_plot)

}, height = function() {
  Data_whole <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  if (is.null(Data_whole)) return(400)
  n_genes <- length(input$gene_survival_interactive_1_P_1_control_forest)
  cancer_counts <- table(Data_whole$Cancers)
  n_cancers <- sum(cancer_counts >= 50)
  if (n_cancers <= 1) n_cancers <- 1
  estimated_items <- 1 + 1 + (n_cancers - 1) + n_genes
  calculated_height <- max(450, estimated_items * 35 + 200)
  return(calculated_height)
})
