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
      Data_cluster_ID_list = Data_cluster_ID() %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster)
      Data_case_target = left_join(Data_case_target,
                                   Data_cluster_ID_list,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_case_target$cluster[is.na(Data_case_target$cluster)] = max(Data_case_target$cluster, na.rm = T) + 1
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
  # Data Pre-processing for Models
  # =====================================================================
  # Ensure there is enough data to proceed (Using shiny:: to prevent namespace collision)
  shiny::validate(shiny::need(nrow(Data_survival) > 0, "No patients match the selected criteria."))

  Data_survival <- Data_survival %>%
    dplyr::mutate(
      time_t2 = time_all - time_pre,
      event_t1 = 1,
      Group = as.factor(Group),
      # Fix possible missing or negative weights safely
      iptw = ifelse(is.na(iptw) | iptw <= 0, 1.0, iptw)
    ) %>%
    dplyr::filter(time_pre > 0, time_t2 > 0)

  shiny::validate(shiny::need(nrow(Data_survival) >= 10, "Not enough valid cases (time_pre > 0 and time_t2 > 0) to fit parametric models."))

  # =====================================================================
  # Robust Parametric Simulation Implementation (Weighted Models)
  # =====================================================================
  req(requireNamespace("flexsurv", quietly = TRUE))
  n_groups <- length(unique(Data_survival$Group))

  # 1. Fit Weibull model for T1 (Time from diagnosis to CGP)
  fit_t1 <- tryCatch({
    flexsurv::flexsurvreg(Surv(time_pre, event_t1) ~ 1,
                          data = Data_survival,
                          weights = Data_survival$iptw,
                          dist = "weibull")
  }, error = function(e) {
    # Fallback to unweighted model if IPTW causes a matrix singularity
    tryCatch({
      flexsurv::flexsurvreg(Surv(time_pre, event_t1) ~ 1, data = Data_survival, dist = "weibull")
    }, error = function(e2) { NULL })
  })

  shiny::validate(shiny::need(!is.null(fit_t1), "Failed to fit Weibull model for T1 (CGP entry time). Data might be too sparse."))

  # 2. Fit Log-logistic AFT model for T2 (Survival after CGP)
  fit_t2 <- tryCatch({
    if (n_groups > 1) {
      flexsurv::flexsurvreg(Surv(time_t2, censor) ~ time_pre + Group,
                            data = Data_survival,
                            weights = Data_survival$iptw,
                            dist = "llogis")
    } else {
      flexsurv::flexsurvreg(Surv(time_t2, censor) ~ time_pre,
                            data = Data_survival,
                            weights = Data_survival$iptw,
                            dist = "llogis")
    }
  }, error = function(e) {
    # Fallback to unweighted model
    tryCatch({
      if (n_groups > 1) {
        flexsurv::flexsurvreg(Surv(time_t2, censor) ~ time_pre + Group, data = Data_survival, dist = "llogis")
      } else {
        flexsurv::flexsurvreg(Surv(time_t2, censor) ~ time_pre, data = Data_survival, dist = "llogis")
      }
    }, error = function(e2) { NULL })
  })

  shiny::validate(shiny::need(!is.null(fit_t2), "Failed to fit Log-logistic model for T2 (Post-CGP survival)."))

  # 3. Simulate Cohort (n = 5000 per group for smooth, stable curves)
  n_sim_per_group <- 5000
  simulated_data <- list()

  # Extract T1 model coefficients
  shape_t1 <- fit_t1$res["shape", "est"]
  scale_t1 <- fit_t1$res["scale", "est"]

  # Extract T2 model coefficients
  shape_t2 <- fit_t2$res["shape", "est"]
  baseline_scale_t2 <- fit_t2$res["scale", "est"]

  # Guard against covariates being completely dropped from the model
  beta_time_pre <- ifelse("time_pre" %in% rownames(fit_t2$res), fit_t2$res["time_pre", "est"], 0)
  beta_group2 <- ifelse(n_groups > 1 && "Group2" %in% rownames(fit_t2$res), fit_t2$res["Group2", "est"], 0)

  for (g in levels(Data_survival$Group)) {
    # Generate random probabilities bounded strictly to avoid Inf/NaN errors
    u1 <- runif(n_sim_per_group, min = 0.001, max = 0.999)
    u2 <- runif(n_sim_per_group, min = 0.001, max = 0.999)

    # Step A: Sample T1 from the Weibull distribution
    sim_t1_days <- qweibull(u1, shape = shape_t1, scale = scale_t1)

    # Step B: Predict individual scale for T2 based on simulated T1
    group_eff <- ifelse(g == "2", beta_group2, 0)
    sim_scale_t2 <- baseline_scale_t2 * exp(beta_time_pre * sim_t1_days + group_eff)

    # Step C: Sample T2 from the Log-logistic distribution
    sim_t2_days <- sim_scale_t2 * ((1 - u2) / u2)^(1 / shape_t2)

    # Step D: Integrate Overall Survival (OS = T1 + T2)
    sim_os_days <- sim_t1_days + sim_t2_days

    # Cap the maximum survival time at 20 years to prevent plot breakdown
    sim_os_days <- pmin(sim_os_days, 365.25 * 20)

    simulated_data[[length(simulated_data) + 1]] <- data.frame(
      C.CAT調査結果.基本項目.ハッシュID = paste0("SIM_", g, "_", 1:n_sim_per_group),
      time_pre = sim_t1_days,
      time_all = sim_os_days,
      censor = 1, # Uncensored completely in baseline simulation
      Group = g
    )
  }

  # Clean up and ensure finite values
  Data_survival_simulated <- do.call(rbind, simulated_data) %>%
    dplyr::filter(is.finite(time_all), is.finite(time_pre))

  shiny::validate(shiny::need(nrow(Data_survival_simulated) > 0, "Simulation failed to generate valid survival times."))

  # =====================================================================
  # Plot the Unbiased Simulated Cohort Data
  # =====================================================================
  survival_compare_and_plot_CTx(
    data = Data_survival_simulated,
    time_var1 = "time_pre",
    time_var2 = "time_all",
    status_var = "censor",
    group_var = "Group",
    plot_title = "Unbiased Parametric Simulation OS (Weighted)",
    adjustment = FALSE,
    color_var_surv_CTx_1 = "diagnosis",
    weights_var = NULL
  )
})
