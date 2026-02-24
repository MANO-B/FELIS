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


# output$figure_survival_CTx_interactive_1_control = renderPlot({
#   req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
#       OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control,
#       OUTPUT_DATA$figure_surv_CTx_Data_drug_control)
#
#   Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
#   Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
#   Data_drug = OUTPUT_DATA$figure_surv_CTx_Data_drug_control
#   # ========================================================
#   # IPTW calculation for left truncation bias adjustment
#   # ========================================================
#   ref_surv_list <- list()
#   age_groups <- c("40未満", "40代", "50代", "60代", "70代", "80以上")
#
#   if (!is.null(input$survival_data_source)) {
#     if (input$survival_data_source == "registry") {
#       # Load from JSON data
#       if(exists("Data_age_survival_5_year") && input$registry_cancer_type %in% names(Data_age_survival_5_year)) {
#
#         cancer_data <- Data_age_survival_5_year[[input$registry_cancer_type]]
#
#         # Extract the overall survival rates as fallback
#         fallback_surv <- cancer_data[["全年齢"]]
#
#         for (ag in age_groups) {
#           # Check if age-specific data exists and has exactly 5 years of data
#           if (!is.null(cancer_data[[ag]]) && length(cancer_data[[ag]]) == 5) {
#             ref_surv_list[[ag]] <- as.numeric(cancer_data[[ag]])
#           } else if (!is.null(fallback_surv) && length(fallback_surv) == 5) {
#             # Fallback to "全年齢" (overall) survival rates if age-specific is missing
#             ref_surv_list[[ag]] <- as.numeric(fallback_surv)
#           }
#         }
#       }
#     } else if (input$survival_data_source == "manual_all") {
#       # Apply overall inputs to all age groups
#       vals <- c(input$surv_all_1y, input$surv_all_2y, input$surv_all_3y, input$surv_all_4y, input$surv_all_5y)
#       for (ag in age_groups) ref_surv_list[[ag]] <- vals
#     } else if (input$survival_data_source == "manual_age") {
#       # Apply detailed inputs by age group
#       ref_surv_list[["40未満"]] <- c(input$surv_u40_1y, input$surv_u40_2y, input$surv_u40_3y, input$surv_u40_4y, input$surv_u40_5y)
#       ref_surv_list[["40代"]] <- c(input$surv_40s_1y, input$surv_40s_2y, input$surv_40s_3y, input$surv_40s_4y, input$surv_40s_5y)
#       ref_surv_list[["50代"]] <- c(input$surv_50s_1y, input$surv_50s_2y, input$surv_50s_3y, input$surv_50s_4y, input$surv_50s_5y)
#       ref_surv_list[["60代"]] <- c(input$surv_60s_1y, input$surv_60s_2y, input$surv_60s_3y, input$surv_60s_4y, input$surv_60s_5y)
#       ref_surv_list[["70代"]] <- c(input$surv_70s_1y, input$surv_70s_2y, input$surv_70s_3y, input$surv_70s_4y, input$surv_70s_5y)
#       ref_surv_list[["80以上"]] <- c(input$surv_80s_1y, input$surv_80s_2y, input$surv_80s_3y, input$surv_80s_4y, input$surv_80s_5y)
#     }
#   }
#
#   # If list is successfully constructed
#   if (length(ref_surv_list) > 0) {
#
#     # Ultimate fallback just in case some rare cancers have NO overall data either
#     for (ag in age_groups) {
#       if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- c(10, 5, 3, 2, 1) # Fallback to a poor prognosis
#     }
#
#     Data_survival_interactive <- calculate_iptw_age(Data_survival_interactive, ref_surv_list, time_var = "time_pre", age_var = "症例.基本情報.年齢")
#   } else {
#     Data_survival_interactive$iptw <- 1.0 # Fallback
#   }
#
#   # ID抽出関数を定義
#   extract_group_ids <- function(group_num) {
#     # 初期IDセット
#     IDs <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)
#
#     # 動的に入力名を構築
#     input_prefix <- paste0("gene_survival_interactive_", group_num, "_")
#
#     # 遺伝子フィルタ（P_1: 必須変異1）
#     p1_input <- input[[paste0(input_prefix, "P_1_control")]]
#     if(!all(is.null(p1_input))) {
#       IDs <- intersect(IDs, (Data_MAF_target %>%
#                                dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
#     }
#
#     # 遺伝子フィルタ（P_2: 必須変異2）
#     p2_input <- input[[paste0(input_prefix, "P_2_control")]]
#     if(!all(is.null(p2_input))) {
#       IDs <- intersect(IDs, (Data_MAF_target %>%
#                                dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
#     }
#
#     # 遺伝子除外（W: 除外変異）
#     w_input <- input[[paste0(input_prefix, "W_control")]]
#     if(!all(is.null(w_input))) {
#       IDs <- setdiff(IDs, (Data_MAF_target %>%
#                              dplyr::filter(Hugo_Symbol %in% w_input))$Tumor_Sample_Barcode)
#     }
#
#     # 臨床データフィルタ
#     clinical_filters <- list(
#       A_control = "YoungOld",
#       S_control = "症例.基本情報.性別.名称.",
#       H_control = "Cancers"
#     )
#
#     for(filter_key in names(clinical_filters)) {
#       filter_input <- input[[paste0(input_prefix, filter_key)]]
#       if(!all(is.null(filter_input))) {
#         column_name <- clinical_filters[[filter_key]]
#         filter_expr <- paste0(column_name, " %in% filter_input")
#         IDs <- intersect(IDs, (Data_survival_interactive %>%
#                                  dplyr::filter(!!parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
#       }
#     }
#
#     # 薬剤フィルタ（D）
#     d_input <- input[[paste0(input_prefix, "D_control")]]
#     if(!all(is.null(d_input))) {
#       IDs <- intersect(IDs, (Data_drug %>%
#                                dplyr::filter(Drug %in% d_input))$ID)
#     }
#
#     return(IDs)
#   }
#
#   # Group1とGroup2のIDを取得
#   ID_1 <- extract_group_ids(1)
#   ID_2 <- extract_group_ids(2)
#
#   # データ作成
#   Data_survival_1 <- Data_survival_interactive %>%
#     dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>%
#     dplyr::mutate(Group = 1)
#
#   Data_survival_2 <- Data_survival_interactive %>%
#     dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>%
#     dplyr::mutate(Group = 2)
#
#   Data_survival <- rbind(Data_survival_1, Data_survival_2)
#
#   # サバイバル解析実行
#   survival_compare_and_plot_CTx(
#     data = Data_survival,
#     time_var1 = "time_pre",
#     time_var2 = "time_all",
#     status_var = "censor",
#     group_var = "Group",
#     plot_title = "Survival analisys based on cohort data",
#     adjustment = FALSE,
#     color_var_surv_CTx_1 = "diagnosis",
#     weights_var = "iptw"
#   )
# })

output$figure_survival_CTx_interactive_1_control = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control,
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control)

  # =========================================================================
  # Force Reactivity: Evaluate inputs explicitly so Shiny redraws on changes
  # =========================================================================
  lapply(1:2, function(i) {
    prefix <- paste0("gene_survival_interactive_", i, "_")
    list(
      input[[paste0(prefix, "P_1_control")]],
      input[[paste0(prefix, "P_2_control")]],
      input[[paste0(prefix, "W_control")]],
      input[[paste0(prefix, "A_control")]],
      input[[paste0(prefix, "S_control")]],
      input[[paste0(prefix, "H_control")]],
      input[[paste0(prefix, "D_control")]]
    )
  })

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  Data_drug = OUTPUT_DATA$figure_surv_CTx_Data_drug_control

  # ========================================================
  # Prepare Reference Survival Data (Macro Data)
  # ========================================================
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

  # Fallback generation for empty lists
  if (length(ref_surv_list) == 0) {
    for (ag in age_groups) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
  } else {
    for (ag in age_groups) {
      if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }

  # Pre-calculate age class for simulation sampling
  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = as.numeric(症例.基本情報.年齢),
      age_class = dplyr::case_when(
        age_num < 40 ~ "40未満",
        age_num < 50 ~ "40代",
        age_num < 60 ~ "50代",
        age_num < 70 ~ "60代",
        age_num < 80 ~ "70代",
        age_num >= 80 ~ "80以上",
        TRUE ~ "全年齢"
      )
    )

  # Define ID extraction function
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

    return(IDs)
  }

  # Get IDs for Group 1 and Group 2
  ID_1 <- extract_group_ids(1)
  ID_2 <- extract_group_ids(2)

  Data_survival_1 <- Data_survival_interactive %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>% dplyr::mutate(Group = "1")
  Data_survival_2 <- Data_survival_interactive %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>% dplyr::mutate(Group = "2")

  Data_survival <- rbind(Data_survival_1, Data_survival_2)
  req(nrow(Data_survival) > 0)

  # =====================================================================
  # Data Pre-processing for Anchored Models
  # =====================================================================
  # Create a baseline reference group (EP_treat == 0)
  Data_EP0 <- Data_survival_interactive %>%
    dplyr::filter(EP_treat == 0) %>%
    dplyr::mutate(Group = "EP0")

  # Combine EP0, Group 1, and Group 2 for modeling
  Data_model <- rbind(Data_EP0, Data_survival_1, Data_survival_2) %>%
    dplyr::mutate(time_t2 = time_all - time_pre) %>%
    dplyr::filter(time_pre > 0, time_t2 > 0)

  shiny::validate(shiny::need(nrow(Data_model) > 0, "No valid data to fit the model."))

  # Set EP0 as the reference baseline level
  Data_model$Group <- factor(Data_model$Group, levels = c("EP0", "1", "2"))

  shiny::validate(shiny::need(nrow(Data_model %>% dplyr::filter(Group == "EP0")) >= 5,
                              "Not enough EP_treat == 0 patients to establish the reference baseline."))

  # =====================================================================
  # Fit Conditional T2 Model to Bypass Dependent Truncation (Kendall tau != 0)
  # =====================================================================
  # By conditioning on time_pre, we neutralize the dependent left truncation bias.
  # We use Log-logistic distribution to extract the unbiased Acceleration Factor (AF).
  req(requireNamespace("flexsurv", quietly = TRUE))

  fit_cgp_t2 <- tryCatch({
    flexsurv::flexsurvreg(Surv(time_t2, censor) ~ time_pre + Group, data = Data_model, dist = "llogis")
  }, error = function(e) { NULL })

  shiny::validate(shiny::need(!is.null(fit_cgp_t2), "Failed to fit conditional T2 model. Data might be too sparse."))

  # =====================================================================
  # Extract Baseline Parameters from Macro Data (Log-logistic OS)
  # =====================================================================
  macro_surv <- ref_surv_list[["全年齢"]]
  if(is.null(macro_surv)) macro_surv <- ref_surv_list[[1]]

  S_t <- macro_surv[1:5] / 100
  S_t <- pmax(pmin(S_t, 0.999), 0.001)
  t_points <- 1:5

  # Linearize Log-logistic: log(1/S(t) - 1) = p * log(lambda) + p * log(t)
  y_llogis <- log(1/S_t - 1)
  x_llogis <- log(t_points)
  fit_macro <- lm(y_llogis ~ x_llogis)

  # Extract shape (p) and scale (alpha = 1/lambda)
  shape_macro <- coef(fit_macro)[2]
  log_lambda_macro <- coef(fit_macro)[1] / shape_macro
  scale_macro_years <- exp(-log_lambda_macro) # alpha = 1 / lambda
  scale_macro_days <- scale_macro_years * 365.25

  # =====================================================================
  # Direct OS Simulation via Log-logistic Distribution
  # =====================================================================
  n_sim_per_group <- 5000
  simulated_data <- list()

  # Extract log(Acceleration Factor) for Group 1 and 2 compared to EP0
  coef_g1 <- ifelse("Group1" %in% rownames(fit_cgp_t2$res), fit_cgp_t2$res["Group1", "est"], 0)
  coef_g2 <- ifelse("Group2" %in% rownames(fit_cgp_t2$res), fit_cgp_t2$res["Group2", "est"], 0)

  for (g in c("1", "2")) {
    # Check if the requested group actually had data; skip if empty
    if (g == "1" && nrow(Data_survival_1) == 0) next
    if (g == "2" && nrow(Data_survival_2) == 0) next

    # Bounded uniform to avoid Inf/NaN
    u <- runif(n_sim_per_group, min = 0.001, max = 0.999)

    # Apply the unbiased Acceleration Factor to the Macro OS scale
    af_g <- ifelse(g == "1", coef_g1, coef_g2)
    scale_g <- scale_macro_days * exp(af_g)

    # Inverse Transform Sampling for Log-logistic OS
    # Formula: t = scale * ((1 - u) / u)^(1 / shape)
    sim_os_days <- scale_g * ((1 - u) / u)^(1 / shape_macro)
    sim_os_days <- pmin(sim_os_days, 365.25 * 20) # Cap at 20 years

    simulated_data[[length(simulated_data) + 1]] <- data.frame(
      C.CAT調査結果.基本項目.ハッシュID = paste0("SIM_", g, "_", 1:n_sim_per_group),
      time_pre = 0, # Unused, set to 0
      time_all = sim_os_days,
      censor = 1, # All simulated baseline patients die eventually
      Group = g
    )
  }

  Data_survival_simulated <- do.call(rbind, simulated_data)
  shiny::validate(shiny::need(!is.null(Data_survival_simulated) && nrow(Data_survival_simulated) > 0,
                              "Simulation yielded no data."))

  # =====================================================================
  # Plot the Unbiased Simulated Cohort Data
  # =====================================================================
  survival_compare_and_plot_CTx(
    data = Data_survival_simulated,
    time_var1 = "time_pre", # Not used for adjustment=FALSE
    time_var2 = "time_all",
    status_var = "censor",
    group_var = "Group",
    plot_title = "Anchored OS Simulation (Log-logistic, AF from conditional T2)",
    adjustment = FALSE,
    color_var_surv_CTx_1 = "diagnosis",
    weights_var = NULL
  )
})
