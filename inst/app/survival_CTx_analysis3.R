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

  # ========================================================
  # Calculate IPTW (Inverse Probability of Treatment Weighting)
  # ========================================================
  # This step projects the CGP cohort to match the general population survival baseline
  if (length(ref_surv_list) > 0) {
    # NOTE: Assuming 'calculate_iptw_age' or 'calculate_iptw_loglogistic_age' is defined
    Data_survival_interactive <- calculate_iptw_age(Data_survival_interactive, ref_surv_list, time_var = "time_pre", age_var = "症例.基本情報.年齢")
  } else {
    Data_survival_interactive$iptw <- 1.0
  }

  # Sanitize IPTW to ensure it is a pure numeric vector
  Data_survival_interactive$iptw <- as.numeric(unlist(Data_survival_interactive$iptw))
  Data_survival_interactive$iptw <- ifelse(is.na(Data_survival_interactive$iptw) | Data_survival_interactive$iptw <= 0, 1.0, Data_survival_interactive$iptw)

  # Pre-calculate age class for safety
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

  # ========================================================
  # Extract Group IDs
  # ========================================================
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

  # =====================================================================
  # Data Pre-processing for Anchored Models
  # =====================================================================
  # Create a baseline reference group (EP_treat == 0)
  Data_EP0 <- Data_survival_interactive %>%
    dplyr::filter(EP_treat == 0) %>%
    dplyr::mutate(Group = "EP0")

  # Combine EP0, Group 1, and Group 2 for OS modeling
  Data_model <- rbind(Data_EP0, Data_survival_1, Data_survival_2) %>%
    dplyr::mutate(time_t2 = time_all - time_pre) %>%
    dplyr::filter(time_pre > 0, time_all > time_pre)

  shiny::validate(shiny::need(nrow(Data_model) > 0, "No valid data to fit the model."))

  # Set EP0 as the reference baseline level
  Data_model$Group <- factor(Data_model$Group, levels = c("EP0", "1", "2"))

  shiny::validate(shiny::need(nrow(Data_model %>% dplyr::filter(Group == "EP0")) >= 5,
                              "Not enough EP_treat == 0 patients to establish the reference baseline."))

  # =====================================================================
  # Fit Multivariate Left-Truncated OS Model with IPTW (Doubly Robust)
  # =====================================================================
  req(requireNamespace("flexsurv", quietly = TRUE))

  # 1. Define potential confounding factors to include in the multivariate model
  # Ensure the column names exactly match your dataset
  potential_covariates <- c("age_num", "`症例.基本情報.性別.名称.`", "Cancers")

  # 2. Dynamically build the formula by checking if columns exist and have sufficient variance
  # (Prevents model crash if a variable has only 1 level in the filtered data)
  valid_covariates <- c("Group")
  for (cov in potential_covariates) {
    clean_cov <- gsub("`", "", cov) # Remove backticks for checking
    if (clean_cov %in% colnames(Data_model)) {
      # Only add if there are at least 2 distinct values (to avoid singularity)
      if (length(unique(na.omit(Data_model[[clean_cov]]))) > 1) {
        valid_covariates <- c(valid_covariates, cov)
      }
    }
  }

  # Construct the multivariate formula: Surv(...) ~ Group + age_num + Cancers...
  formula_str <- paste("Surv(time_pre, time_all, censor) ~", paste(valid_covariates, collapse = " + "))
  os_formula <- as.formula(formula_str)
  dist_choice <- "llogis"
  # 3. Fit the Doubly Robust AFT Model
  fit_os <- tryCatch({
    flexsurv::flexsurvreg(os_formula,
                          data = Data_model,
                          weights = iptw,
                          dist = dist_choice)
  }, error = function(e) {
    # Fallback to unweighted if the multi-dimensional matrix is singular
    tryCatch({
      flexsurv::flexsurvreg(os_formula, data = Data_model, dist = dist_choice)
    }, error = function(e2) { NULL })
  })

  # Fail-safe: Fallback to Log-logistic if Weibull fails in multivariate space
  if (is.null(fit_os) && dist_choice == "weibull") {
    dist_choice <- "llogis"
    fit_os <- tryCatch({
      flexsurv::flexsurvreg(os_formula, data = Data_model, weights = iptw, dist = "llogis")
    }, error = function(e) { NULL })
  }

  shiny::validate(shiny::need(!is.null(fit_os), paste("Failed to fit the multivariate left-truncated OS model. Too many rare cancer types might cause a singular matrix.")))

  # 4. Extract the MULTIVARIATE-ADJUSTED pure Acceleration Factors
  # These coefficients now represent the isolated effect of the group, independent of age/histology
  af_g1 <- ifelse("Group1" %in% rownames(fit_os$res), fit_os$res["Group1", "est"], 0)
  af_g2 <- ifelse("Group2" %in% rownames(fit_os$res), fit_os$res["Group2", "est"], 0)
  OUTPUT_DATA$fit_os <- fit_os

  # =====================================================================
  # Extract Age-Stratified Baseline Parameters from Macro Data (Registry)
  # =====================================================================
  # [FIX] Instead of "全年齢" (All Ages), we calculate parameters for EACH age group
  macro_models <- list()
  for (ag in names(ref_surv_list)) {
    S_t <- ref_surv_list[[ag]][1:5] / 100
    S_t <- pmax(pmin(S_t, 0.999), 0.001)

    y_llogis <- log(1/S_t - 1)
    x_llogis <- log(1:5)
    fit_macro <- lm(y_llogis ~ x_llogis)

    shape <- coef(fit_macro)[2]
    scale_days <- exp(-coef(fit_macro)[1] / shape) * 365.25

    macro_models[[ag]] <- list(shape = shape, scale = scale_days)
  }

  # =====================================================================
  # Simulation: Absolute OS Nearest Neighbor Matching (Age-Stratified)
  # =====================================================================
  simulated_data <- list()

  # [修正箇所] 先生のご指摘通り、不完全な比率の適用を防ぐため、
  # 比率抽出のドナープールは「死亡確認例（censor == 1）」のみに限定します。
  emp_ep0 <- Data_EP0 %>%
    dplyr::filter(time_all > 0, time_pre > 0, time_all > time_pre, censor == 1) %>% # censor == 1 を追加
    dplyr::mutate(
      actual_t2 = time_all - time_pre,
      t2_ratio = actual_t2 / time_all
    ) %>%
    dplyr::arrange(time_all)

  # EP0群の死亡例が極端に少ない場合のフォールバック（全体から死亡例を探す）
  if(nrow(emp_ep0) < 5) {
    emp_ep0 <- Data_model %>%
      dplyr::filter(time_all > 0, time_pre > 0, time_all > time_pre, censor == 1) %>% # censor == 1 を追加
      dplyr::mutate(
        actual_t2 = time_all - time_pre,
        t2_ratio = actual_t2 / time_all
      ) %>%
      dplyr::arrange(time_all)
  }

  # 万が一、データ全体でも死亡例が5例未満という異常事態のための最終安全装置
  if(nrow(emp_ep0) < 5) {
    emp_ep0 <- Data_model %>%
      dplyr::filter(time_all > 0, time_pre > 0, time_all > time_pre) %>% # ここだけは打ち切りも含める
      dplyr::mutate(
        actual_t2 = time_all - time_pre,
        t2_ratio = actual_t2 / time_all
      ) %>%
      dplyr::arrange(time_all)
  }

  actual_os_array <- emp_ep0$time_all
  actual_t2_ratio_array <- emp_ep0$t2_ratio

  for (g in c("1", "2")) {
    group_data <- if(g == "1") Data_survival_1 else Data_survival_2
    n_sim <- nrow(group_data)

    if (n_sim == 0) next

    # 1. Generate Age-Stratified Baseline OS (Deterministic Quantile Sampling)
    # [FIX] Replace 'runif' with evenly spaced quantiles to eliminate simulation variance
    sim_os_macro <- numeric(n_sim)

    unique_ages <- unique(group_data$age_class)
    for (ag in unique_ages) {
      idx <- which(group_data$age_class == ag)
      n_ag <- length(idx)

      if (n_ag == 0) next

      # Generate evenly spaced probabilities: (0.5/n, 1.5/n, ..., (n-0.5)/n)
      u_ag <- (1:n_ag - 0.5) / n_ag
      # Bound strictly to avoid Inf at absolute extremes if n_ag is very small
      u_ag <- pmax(pmin(u_ag, 0.999), 0.001)

      ag_model <- ag
      if (is.null(macro_models[[ag_model]])) ag_model <- "全年齢"
      if (is.null(macro_models[[ag_model]])) ag_model <- names(macro_models)[1] # Fallback

      shape_m <- macro_models[[ag_model]]$shape
      scale_m <- macro_models[[ag_model]]$scale

      # Calculate deterministic expected OS values directly from the CDF
      sim_os_macro[idx] <- scale_m * ((1 - u_ag) / u_ag)^(1 / shape_m)
    }

    # 2. Nearest Neighbor Matching based on ABSOLUTE OS length
    pos <- findInterval(sim_os_macro, actual_os_array, all.inside = TRUE)
    diff1 <- abs(sim_os_macro - actual_os_array[pos])
    diff2 <- abs(sim_os_macro - actual_os_array[pos + 1])
    idx <- ifelse(diff1 <= diff2, pos, pos + 1)

    matched_t2_ratio <- actual_t2_ratio_array[idx]
    sim_t2_base <- sim_os_macro * matched_t2_ratio

    # 3. Apply Mutation/Treatment Effects (AF) to the ENTIRE OS
    af <- ifelse(g == "1", af_g1, af_g2)
    sim_os_treated <- sim_os_macro * exp(af)

    # 4. Reconstruct timelines expanding proportionally
    sim_t2_treated <- sim_t2_base * exp(af)
    sim_t1_treated <- sim_os_treated - sim_t2_treated

    # 10-Year Administrative Censoring
    cap_days <- 365.25 * 10
    final_censor <- ifelse(sim_os_treated > cap_days, 0, 1)
    final_os <- pmin(sim_os_treated, cap_days)
    final_t1 <- pmin(sim_t1_treated, cap_days - 0.1)

    simulated_data[[length(simulated_data) + 1]] <- data.frame(
      C.CAT調査結果.基本項目.ハッシュID = group_data$C.CAT調査結果.基本項目.ハッシュID,
      time_pre = final_t1,
      time_all = final_os,
      censor = final_censor,
      Group = g
    )
  }

  Data_survival_simulated <- do.call(rbind, simulated_data) %>%
    dplyr::filter(is.finite(time_all))

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
    plot_title = "Unbiased OS Simulation (IPTW + Left-Truncated Llogis + Rank Match)",
    adjustment = FALSE,
    color_var_surv_CTx_1 = "diagnosis",
    weights_var = NULL
  )
})

output$forest_plot_whole_cohort <- renderPlot({
  # 必要なデータと入力変数が揃っているか確認
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control)

  Data_whole <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control

  # =====================================================================
  # 1. データの前処理（全集団）
  # =====================================================================
  # 左側切断の条件（time_pre > 0, time_all > time_pre）を満たす全患者を抽出
  Data_forest <- Data_whole %>%
    dplyr::filter(time_pre > 0, time_all > time_pre) %>%
    dplyr::mutate(
      age_num = as.numeric(症例.基本情報.年齢),
      # 推奨治療の有無をわかりやすいラベルに変換（リファレンスを"No"にするためfactor化）
      Treatment = factor(ifelse(EP_treat == 1, "Received", "Not_Received"), levels = c("Not_Received", "Received")),
      Sex = factor(`症例.基本情報.性別.名称.`)
    )

  shiny::validate(shiny::need(nrow(Data_forest) > 50, "Not enough valid cases for multivariate analysis."))

  # =====================================================================
  # 2. 注目する変異（G1条件）のフラグ立て
  # =====================================================================
  # G1の抽出条件（P1, P2, Wなど）に合致した患者IDのリスト（以前のextract_group_ids関数の結果など）を取得
  # ※ ここでは 'ID_1' がスコープ内にあるか、再計算して取得する想定です
  if (exists("ID_1")) {
    Data_forest$Target_Mutation <- factor(
      ifelse(Data_forest$C.CAT調査結果.基本項目.ハッシュID %in% ID_1, "Mutated", "Wild-type"),
      levels = c("Wild-type", "Mutated") # Wild-typeをリファレンス（基準）にする
    )
  } else {
    # ID_1が取得できない場合の安全装置
    Data_forest$Target_Mutation <- factor("Wild-type", levels = c("Wild-type", "Mutated"))
  }

  # =====================================================================
  # 3. 多変量モデルの構築（IPTW + 左側切断 + 対数ロジスティックAFT）
  # =====================================================================
  # IPTWが計算されていなければ 1.0 とする
  if (!"iptw" %in% colnames(Data_forest)) Data_forest$iptw <- 1.0
  Data_forest$iptw <- as.numeric(unlist(Data_forest$iptw))
  Data_forest$iptw <- ifelse(is.na(Data_forest$iptw) | Data_forest$iptw <= 0, 1.0, Data_forest$iptw)

  req(requireNamespace("flexsurv", quietly = TRUE))

  # 組織型（Cancers）を含めると変数が多すぎてSingular matrixになる場合は外す
  # formula_str <- "Surv(time_pre, time_all, censor) ~ Target_Mutation + Treatment + age_num + Sex + Cancers"
  formula_str <- "Surv(time_pre, time_all, censor) ~ Target_Mutation + Treatment + age_num + Sex"

  fit_forest <- tryCatch({
    flexsurv::flexsurvreg(as.formula(formula_str),
                          data = Data_forest,
                          weights = iptw,
                          dist = "llogis")
  }, error = function(e) { NULL })

  shiny::validate(shiny::need(!is.null(fit_forest), "Multivariate model failed to converge. Please check covariate variance."))

  # =====================================================================
  # 4. フォレストプロット用データの抽出と整形
  # =====================================================================
  res <- fit_forest$res
  # 'shape' と 'scale' はベースラインハザードのパラメータなので除外
  cov_names <- rownames(res)[!rownames(res) %in% c("shape", "scale")]

  plot_data <- data.frame(
    Variable = cov_names,
    # AFTモデルの係数を指数変換して「時間比（Time Ratio: 加速係数）」に戻す
    Estimate = exp(res[cov_names, "est"]),
    Lower = exp(res[cov_names, "L95%"]),
    Upper = exp(res[cov_names, "U95%"])
  )

  # 変数名をグラフ上で読みやすくリネーム（適宜変更してください）
  plot_data$Variable <- gsub("Target_MutationMutated", "Target Mutation (Mutated vs WT)", plot_data$Variable)
  plot_data$Variable <- gsub("TreatmentReceived", "Recommended Treatment (Recv vs Not)", plot_data$Variable)
  plot_data$Variable <- gsub("age_num", "Age (per +1 year)", plot_data$Variable)
  plot_data$Variable <- gsub("Sex男性", "Sex (Male vs Female)", plot_data$Variable)

  # =====================================================================
  # 5. ggplot2 による描画
  # =====================================================================

  ggplot(plot_data, aes(x = Estimate, y = reorder(Variable, Estimate))) +
    # エラーバーとポイントの描画
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#2c3e50", size = 0.8) +
    geom_point(size = 4, color = "#e74c3c") +
    # TR = 1.0 (影響なし) の垂直線
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7f8c8d", size = 1) +
    # AFTの時間比（倍率）なので、対数スケールにすることで対称的な視覚表現になる
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5)) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Multivariate Analysis of Overall Survival (Entire Cohort)",
      subtitle = "Time Ratio (TR) > 1: Prolonged Survival | TR < 1: Shortened Survival",
      x = "Time Ratio (Log Scale, 95% CI)",
      y = ""
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(margin = margin(t = 10))
    )
})


output$forest_plot_multivariate = renderPlot({
  # =========================================================================
  # Require necessary data and packages
  # =========================================================================
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control)

  shiny::validate(shiny::need(requireNamespace("patchwork", quietly = TRUE),
                              "Please install the 'patchwork' package. Run: install.packages('patchwork')"))

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  selected_genes <- input$gene_survival_interactive_1_P_1_control_forest

  # =========================================================================
  # Prepare Reference Survival Data for IPTW
  # =========================================================================
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
    for (ag in age_groups) {
      if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }

  # Calculate IPTW
  if (length(ref_surv_list) > 0) {
    Data_survival_interactive <- calculate_iptw_age(Data_survival_interactive, ref_surv_list, time_var = "time_pre", age_var = "症例.基本情報.年齢")
  } else {
    Data_survival_interactive$iptw <- 1.0
  }

  Data_survival_interactive$iptw <- as.numeric(unlist(Data_survival_interactive$iptw))
  Data_survival_interactive$iptw <- ifelse(is.na(Data_survival_interactive$iptw) | Data_survival_interactive$iptw <= 0, 1.0, Data_survival_interactive$iptw)

  # =========================================================================
  # Prepare Data for Multivariate Modeling (Robust Factor Handling)
  # =========================================================================
  Data_forest <- Data_survival_interactive %>%
    dplyr::filter(time_pre > 0, time_all > time_pre) %>%
    dplyr::mutate(
      age_num = as.numeric(症例.基本情報.年齢),
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.character(Cancers)
    )

  shiny::validate(shiny::need(nrow(Data_forest) > 50, "Not enough valid cases for multivariate analysis."))
  total_cohort_N <- nrow(Data_forest)

  # Set Sex reference
  if(length(levels(Data_forest$Sex)) > 1) {
    top_sex <- names(sort(table(Data_forest$Sex), decreasing = TRUE))[1]
    Data_forest$Sex <- relevel(Data_forest$Sex, ref = top_sex)
  }

  # [MODIFIED] Lump rare cancers (< 50 cases) into "Others"
  cancer_counts <- table(Data_forest$Cancers)
  rare_cancers <- names(cancer_counts)[cancer_counts < 50]
  Data_forest$Cancers[Data_forest$Cancers %in% rare_cancers] <- "Others"
  Data_forest$Cancers <- as.factor(Data_forest$Cancers)

  # Ensure the most frequent cancer type is the reference level
  top_cancer <- names(sort(table(Data_forest$Cancers), decreasing = TRUE))[1]
  Data_forest$Cancers <- relevel(Data_forest$Cancers, ref = top_cancer)

  # =========================================================================
  # Dynamically add selected genes
  # =========================================================================
  valid_gene_covariates <- c()
  if (length(selected_genes) > 0) {
    for (gene in selected_genes) {
      mutated_IDs <- (Data_MAF_target %>% dplyr::filter(Hugo_Symbol == gene))$Tumor_Sample_Barcode
      col_name <- paste0("Gene_", make.names(gene))

      Data_forest[[col_name]] <- factor(
        ifelse(Data_forest$C.CAT調査結果.基本項目.ハッシュID %in% mutated_IDs, "Mutated", "WT"),
        levels = c("WT", "Mutated")
      )

      if (length(unique(Data_forest[[col_name]])) > 1) {
        valid_gene_covariates <- c(valid_gene_covariates, col_name)
      }
    }
  }

  # =========================================================================
  # Build and Fit Model
  # =========================================================================
  base_covariates <- c("age_num", "Sex", "Cancers")
  all_covariates <- c()

  for(cov in base_covariates) {
    if(length(unique(na.omit(Data_forest[[cov]]))) > 1) {
      all_covariates <- c(all_covariates, cov)
    }
  }
  all_covariates <- c(all_covariates, valid_gene_covariates)

  shiny::validate(shiny::need(length(all_covariates) > 0, "No valid covariates found."))
  formula_str <- paste("Surv(time_pre, time_all, censor) ~", paste(all_covariates, collapse = " + "))

  req(requireNamespace("flexsurv", quietly = TRUE))
  fit_forest <- tryCatch({
    flexsurv::flexsurvreg(as.formula(formula_str), data = Data_forest, weights = iptw, dist = "llogis")
  }, error = function(e) {
    tryCatch({ flexsurv::flexsurvreg(as.formula(formula_str), data = Data_forest, dist = "llogis") }, error = function(e2) { NULL })
  })

  shiny::validate(shiny::need(!is.null(fit_forest), "Multivariate model failed to converge."))

  # =====================================================================
  # Extract and Format Data for Plot and Table
  # =====================================================================
  res <- fit_forest$res
  cov_names <- rownames(res)[!rownames(res) %in% c("shape", "scale")]

  plot_data <- data.frame(
    Variable_Raw = cov_names,
    Estimate = exp(res[cov_names, "est"]),
    Lower = exp(res[cov_names, "L95%"]),
    Upper = exp(res[cov_names, "U95%"]),
    Category = NA, Label = NA, Positive_N = NA, Total_N = total_cohort_N, Text_CI = NA
  )

  # [MODIFIED] Adjust Estimate and CI for Age (per +10 years)
  for (i in 1:nrow(plot_data)) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^age_num", var)) {
      # Multiply log-coefficient by 10, then exponentiate to get TR per 10 years
      plot_data$Estimate[i] <- exp(res[var, "est"] * 10)
      plot_data$Lower[i] <- exp(res[var, "L95%"] * 10)
      plot_data$Upper[i] <- exp(res[var, "U95%"] * 10)
    }
  }

  # Format CI text after age correction
  plot_data$Text_CI <- sprintf("%.2f (%.2f-%.2f)", plot_data$Estimate, plot_data$Lower, plot_data$Upper)

  for (i in 1:nrow(plot_data)) {
    var <- plot_data$Variable_Raw[i]

    if (grepl("^age_num", var)) {
      plot_data$Label[i] <- "Age (per +10 years)" # Updated label
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

  # =====================================================================
  # Build Two-Panel Plot (Forest Plot + Text Table) using Patchwork
  # =====================================================================
  library(ggplot2)
  library(patchwork)

  # Left: Forest Plot
  p_left <- ggplot(plot_data, aes(x = Estimate, y = Label, color = Category)) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.3, size = 0.8) +
    geom_point(size = 4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7f8c8d", size = 1) +
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

  # Right: Text Table
  p_right <- ggplot(plot_data, aes(y = Label)) +
    geom_text(aes(x = 0, label = Total_N), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 1, label = Text_PosN), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 2, label = Text_CI), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    scale_x_continuous(limits = c(-0.5, 2.5), breaks = c(0, 1, 2), labels = c("Total N", "Positive N", "TR (95% CI)")) +
    theme_minimal(base_size = 14) +
    labs(x = "", y = "") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(face = "bold", color = "#2c3e50"),
      axis.ticks = element_blank()
    )

  # Merge with patchwork
  # [MODIFIED] Added N >= 50 condition to the subtitle
  final_plot <- p_left + p_right + plot_layout(widths = c(2, 1.2)) +
    plot_annotation(
      title = "Independent Prognostic Impact (Multivariate AFT Model)",
      subtitle = "Model adjusted for Age, Sex, Histology (N >= 50), and concurrent Mutations (IPTW applied)\nTR > 1: Prolonged Survival | TR < 1: Shortened Survival",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "#34495e")
      )
    )

  return(final_plot)

}, height = function() {
  # =====================================================================
  # Dynamically calculate plot height
  # =====================================================================
  Data_whole <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  if (is.null(Data_whole)) return(400)

  n_genes <- length(input$gene_survival_interactive_1_P_1_control_forest)

  # [MODIFIED] Accurately count histology levels (>= 50 cases)
  cancer_counts <- table(Data_whole$Cancers)
  n_cancers <- sum(cancer_counts >= 50)
  if (n_cancers <= 1) n_cancers <- 1

  # Age(1) + Sex(1) + Histologies(n_cancers - 1 ref) + Genes(n_genes)
  estimated_items <- 1 + 1 + (n_cancers - 1) + n_genes

  calculated_height <- max(450, estimated_items * 35 + 200)

  return(calculated_height)
})

# =========================================================================
# Server Logic for Simulation Study (Fixed Version)
# =========================================================================

# Helper function for censoring scale estimation (Fixed interval for 'Days' scale)
find_censoring_scale <- function(T2, target_rate, pattern) {
  obj_func <- function(scale_param) {
    if (pattern == "indep") {
      C2 <- rexp(length(T2), rate = 1 / scale_param)
    } else if (pattern == "early") {
      C2 <- rweibull(length(T2), shape = 0.5, scale = scale_param)
    } else if (pattern == "late") {
      C2 <- rweibull(length(T2), shape = 3.0, scale = scale_param)
    } else if (pattern == "ushape") {
      u <- runif(length(T2))
      C2 <- ifelse(u < 0.5,
                   rweibull(length(T2), shape = 0.5, scale = scale_param),
                   rweibull(length(T2), shape = 3.0, scale = scale_param * 2))
    }
    censored <- sum(C2 < T2) / length(T2)
    return(censored - target_rate)
  }

  # [FIX] Expanded the interval to accommodate 'Days' (e.g., up to 1,000,000 days)
  optimal <- tryCatch({
    uniroot(obj_func, interval = c(1, 1000000))$root
  }, error = function(e) { 365.25 * 5 }) # Fallback to 5 years

  return(optimal)
}

observeEvent(input$run_sim, {
  req(input$sim_n, input$sim_true_af, input$sim_cens_rate, input$sim_mut_freq)

  showNotification("Running multivariate simulation...", type = "message", duration = 3)

  N_target <- input$sim_n
  True_AF_X <- input$sim_true_af
  Mut_Freq <- input$sim_mut_freq / 100
  t1_pat <- input$sim_t1_pattern
  cens_pat <- input$sim_cens_pattern
  cens_rate <- input$sim_cens_rate / 100

  # =========================================================================
  # 1. Generate Macro Population with Covariates
  # =========================================================================
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
    ifelse(Histology == "READ", True_AF_Hist_READ,
           ifelse(Histology == "COADREAD", True_AF_Hist_COADREAD, 1.0))

  shape_base <- 1.5
  scale_base <- 3.0 * 365.25

  u <- runif(N_macro)
  # [FIX] Bound 'u' to prevent infinite survival times
  u <- pmax(pmin(u, 0.999), 0.001)

  T_base <- scale_base * ((1 - u) / u)^(1 / shape_base)
  T_true <- T_base * AF_bg * ifelse(X == 1, True_AF_X, 1.0)

  # =========================================================================
  # 2. Generate T1 (Time to CGP) and Sample Cohort
  # =========================================================================
  if (t1_pat == "indep") {
    T1 <- runif(N_macro, 30, 365.25 * 3)
  } else if (t1_pat == "real") {
    T1 <- T_true * runif(N_macro, 0.3, 0.8)
  } else if (t1_pat == "rev") {
    ratio <- pmax(0.1, 1 - exp(-T_true / 365.25)) * runif(N_macro, 0.2, 0.9)
    T1 <- T_true * ratio
  }

  cgp_indices <- which(T_true > T1)

  prob_select <- ifelse(Age[cgp_indices] < 60, 0.8, 0.4)
  selected_indices <- sample(cgp_indices, size = min(N_target, length(cgp_indices)), prob = prob_select, replace = FALSE)

  Data_cgp <- data.frame(
    ID = 1:length(selected_indices),
    Age = Age[selected_indices],
    Sex = factor(Sex[selected_indices], levels = c("Male", "Female")),
    Histology = factor(Histology[selected_indices], levels = c("COAD", "READ", "COADREAD")),
    X = factor(ifelse(X[selected_indices] == 1, "Mutated", "WT"), levels = c("WT", "Mutated")),
    T_true = T_true[selected_indices],
    T1 = T1[selected_indices]
  )
  Data_cgp$T2_true <- Data_cgp$T_true - Data_cgp$T1

  # =========================================================================
  # 3. Apply Informative Censoring (C2)
  # =========================================================================
  scale_opt <- find_censoring_scale(Data_cgp$T2_true, cens_rate, cens_pat)

  if (cens_pat == "indep") { C2 <- rexp(nrow(Data_cgp), rate = 1 / scale_opt)
  } else if (cens_pat == "early") { C2 <- rweibull(nrow(Data_cgp), shape = 0.5, scale = scale_opt)
  } else if (cens_pat == "late") { C2 <- rweibull(nrow(Data_cgp), shape = 3.0, scale = scale_opt)
  } else if (cens_pat == "ushape") {
    u_c <- runif(nrow(Data_cgp))
    C2 <- ifelse(u_c < 0.5, rweibull(nrow(Data_cgp), shape = 0.5, scale = scale_opt), rweibull(nrow(Data_cgp), shape = 3.0, scale = scale_opt * 2))
  }

  Data_cgp$T2_obs <- pmin(Data_cgp$T2_true, C2)

  # [FIX] Renamed to 'Event' to clarify 1 = Death, 0 = Censored
  Data_cgp$Event <- ifelse(Data_cgp$T2_true <= C2, 1, 0)
  Data_cgp$T_obs <- Data_cgp$T1 + Data_cgp$T2_obs

  # Check if events exist to prevent model explosion
  shiny::validate(shiny::need(sum(Data_cgp$Event) > 20, "Simulation generated too few events (deaths). The model cannot be fitted."))

  # =========================================================================
  # 4. Simplified IPTW based on Age
  # =========================================================================
  Age_class_macro <- ifelse(Age < 60, "Young", "Old")
  pt_young <- sum(T_true[Age_class_macro == "Young"])
  pt_old <- sum(T_true[Age_class_macro == "Old"])

  Data_cgp$Age_class <- ifelse(Data_cgp$Age < 60, "Young", "Old")
  weight_young <- pt_young / sum(Data_cgp$Age_class == "Young")
  weight_old <- pt_old / sum(Data_cgp$Age_class == "Old")

  Data_cgp$iptw <- ifelse(Data_cgp$Age_class == "Young", weight_young, weight_old)
  Data_cgp$iptw <- Data_cgp$iptw / mean(Data_cgp$iptw)

  # =========================================================================
  # 5. Fit Proposed Multivariate Left-Truncated AFT Model
  # =========================================================================
  library(flexsurv)
  # [FIX] using 'Event' instead of 'Censor'
  fit_prop <- tryCatch({
    flexsurvreg(Surv(T1, T_obs, Event) ~ X + Age + Sex + Histology, data = Data_cgp, weights = iptw, dist = "llogis")
  }, error = function(e) { NULL })

  shiny::validate(shiny::need(!is.null(fit_prop), "The survival model failed to converge. Please try running the simulation again."))

  # =========================================================================
  # 6. Render Multivariate Results Table
  # =========================================================================
  output$sim_result_table <- renderTable({
    req(fit_prop)
    res <- fit_prop$res

    est_X <- exp(res["XMutated", "est"])
    est_Age <- exp(res["Age", "est"] * 10)
    est_Sex <- exp(res["SexFemale", "est"])
    est_READ <- exp(res["HistologyREAD", "est"])
    est_COADREAD <- exp(res["HistologyCOADREAD", "est"])

    ci_X <- sprintf("%.2f - %.2f", exp(res["XMutated", "L95%"]), exp(res["XMutated", "U95%"]))
    ci_Age <- sprintf("%.2f - %.2f", exp(res["Age", "L95%"]*10), exp(res["Age", "U95%"]*10))
    ci_Sex <- sprintf("%.2f - %.2f", exp(res["SexFemale", "L95%"]), exp(res["SexFemale", "U95%"]))
    ci_READ <- sprintf("%.2f - %.2f", exp(res["HistologyREAD", "L95%"]), exp(res["HistologyREAD", "U95%"]))
    ci_COADREAD <- sprintf("%.2f - %.2f", exp(res["HistologyCOADREAD", "L95%"]), exp(res["HistologyCOADREAD", "U95%"]))

    res_df <- data.frame(
      Covariate = c("Target Gene (Mutated vs WT)",
                    "Age (per +10 years)",
                    "Sex (Female vs Male)",
                    "Histology: READ (vs COAD)",
                    "Histology: COADREAD (vs COAD)"),
      True_AF = c(sprintf("%.2f", True_AF_X),
                  sprintf("%.2f", True_AF_Age_10yr),
                  sprintf("%.2f", True_AF_Sex_Female),
                  sprintf("%.2f", True_AF_Hist_READ),
                  sprintf("%.2f", True_AF_Hist_COADREAD)),
      Estimated_AF = c(sprintf("%.2f", est_X),
                       sprintf("%.2f", est_Age),
                       sprintf("%.2f", est_Sex),
                       sprintf("%.2f", est_READ),
                       sprintf("%.2f", est_COADREAD)),
      `95%_CI` = c(ci_X, ci_Age, ci_Sex, ci_READ, ci_COADREAD),
      Bias = c(sprintf("%+.3f", est_X - True_AF_X),
               sprintf("%+.3f", est_Age - True_AF_Age_10yr),
               sprintf("%+.3f", est_Sex - True_AF_Sex_Female),
               sprintf("%+.3f", est_READ - True_AF_Hist_READ),
               sprintf("%+.3f", est_COADREAD - True_AF_Hist_COADREAD))
    )
    return(res_df)
  }, bordered = TRUE, striped = TRUE, hover = TRUE, width = "100%", align = "c")

  # =========================================================================
  # 7. Render Reconstructed Survival Plot
  # =========================================================================
  output$sim_survival_plot <- renderPlot({
    req(fit_prop)
    library(ggplot2)

    shape_val <- fit_prop$res["shape", "est"]
    log_scale_base <- fit_prop$res["scale", "est"] + (fit_prop$res["Age", "est"] * 60)
    scale_wt <- exp(log_scale_base)

    af_X <- exp(fit_prop$res["XMutated", "est"])
    scale_mut <- scale_wt * af_X

    t_seq <- seq(0, 365.25 * 5, length.out = 100)

    S_wt <- 1 / (1 + (t_seq / scale_wt)^shape_val)
    S_mut <- 1 / (1 + (t_seq / scale_mut)^shape_val)

    plot_df <- data.frame(
      Time_Years = rep(t_seq / 365.25, 2),
      Survival = c(S_wt, S_mut),
      Group = rep(c("Wild-type", "Mutated (Target Gene)"), each = length(t_seq))
    )

    ggplot(plot_df, aes(x = Time_Years, y = Survival, color = Group)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("Wild-type" = "#2980b9", "Mutated (Target Gene)" = "#e74c3c")) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste0("Reconstructed Survival for Target Gene (True AF = ", True_AF_X, ")"),
        subtitle = "Curves represent a 60-year-old Male patient with COAD (Reference profile)",
        x = "Time from Diagnosis (Years)",
        y = "Overall Survival Probability"
      ) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()
      )
  })
})
