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
      OUTPUT_DATA$figure_surv_CTx_Data_drug = Data_drug

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
        adjustment = FALSE
        OUTPUT_DATA$figure_surv_CTx_adjustment = adjustment
        if(is.null(input$color_var_surv_CTx_3)){
          adjustment = TRUE
        } else if(input$color_var_surv_CTx_3 == "Yes"){
          adjustment = TRUE
        } else {
          adjustment = FALSE
        }
        OUTPUT_DATA$figure_surv_CTx_adjustment = adjustment

        # ========================================================
        # IPTW calculation for left truncation bias adjustment
        # ========================================================
        ref_surv_list <- list()
        age_groups <- c("40未満", "40代", "50代", "60代", "70代", "80以上")

        if (!is.null(input$survival_data_source)) {
          if (input$survival_data_source == "registry") {
            # Load from JSON data
            if(exists("Data_age_survival_5_year") && input$registry_cancer_type %in% names(Data_age_survival_5_year)) {

              cancer_data <- Data_age_survival_5_year[[input$registry_cancer_type]]

              # Extract the overall survival rates as fallback
              fallback_surv <- cancer_data[["全年齢"]]

              for (ag in age_groups) {
                # Check if age-specific data exists and has exactly 5 years of data
                if (!is.null(cancer_data[[ag]]) && length(cancer_data[[ag]]) == 5) {
                  ref_surv_list[[ag]] <- as.numeric(cancer_data[[ag]])
                } else if (!is.null(fallback_surv) && length(fallback_surv) == 5) {
                  # Fallback to "全年齢" (overall) survival rates if age-specific is missing
                  ref_surv_list[[ag]] <- as.numeric(fallback_surv)
                }
              }
            }
          } else if (input$survival_data_source == "manual_all") {
            # Apply overall inputs to all age groups
            vals <- c(input$surv_all_1y, input$surv_all_2y, input$surv_all_3y, input$surv_all_4y, input$surv_all_5y)
            for (ag in age_groups) ref_surv_list[[ag]] <- vals
          } else if (input$survival_data_source == "manual_age") {
            # Apply detailed inputs by age group
            ref_surv_list[["40未満"]] <- c(input$surv_u40_1y, input$surv_u40_2y, input$surv_u40_3y, input$surv_u40_4y, input$surv_u40_5y)
            ref_surv_list[["40代"]] <- c(input$surv_40s_1y, input$surv_40s_2y, input$surv_40s_3y, input$surv_40s_4y, input$surv_40s_5y)
            ref_surv_list[["50代"]] <- c(input$surv_50s_1y, input$surv_50s_2y, input$surv_50s_3y, input$surv_50s_4y, input$surv_50s_5y)
            ref_surv_list[["60代"]] <- c(input$surv_60s_1y, input$surv_60s_2y, input$surv_60s_3y, input$surv_60s_4y, input$surv_60s_5y)
            ref_surv_list[["70代"]] <- c(input$surv_70s_1y, input$surv_70s_2y, input$surv_70s_3y, input$surv_70s_4y, input$surv_70s_5y)
            ref_surv_list[["80以上"]] <- c(input$surv_80s_1y, input$surv_80s_2y, input$surv_80s_3y, input$surv_80s_4y, input$surv_80s_5y)
          }
        }

        # If list is successfully constructed
        if (length(ref_surv_list) > 0) {

          # Ultimate fallback just in case some rare cancers have NO overall data either
          for (ag in age_groups) {
            if (is.null(ref_surv_list[[ag]])) ref_surv_list[[ag]] <- c(10, 5, 3, 2, 1) # Fallback to a poor prognosis
          }

          Data_case_target <- calculate_iptw_age(Data_case_target, ref_surv_list, time_var = "time_pre", age_var = "症例.基本情報.年齢")
        } else {
          Data_case_target$iptw <- 1.0 # Fallback
        }
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
        OUTPUT_DATA$figure_surv_CTx_Data_MAF_target = Data_MAF_target
        incProgress(1 / 13)

        Data_survival_interactive = Data_survival
        OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control = Data_survival_interactive
      }
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}


output$figure_survival_CTx_1_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive,
      input$color_var_surv_CTx_2,
      input$color_var_surv_CTx_1)
  # 変数の初期化
  plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive
  group_var <- ""
  plot_title <- ""

  # ユーザーの選択に応じて変数を設定
  if(input$color_var_surv_CTx_2 == "entire"){
    independence_check = with(plot_data %>% sample_n(min(length(plot_data$C.CAT調査結果.基本項目.ハッシュID), 40000)), cKendall(
      trun = time_pre,
      obs = time_all,
      delta = censor,
      method = "MB",
      trans = "linear"))
    group_var <- "1"
    plot_title <- paste0(nrow(plot_data), " patients, cKendall tau= ",
                         format_p(independence_check$PE,digits=2),
                         ", SE= ", format_p(independence_check$SE,digits=2),
                         ", p= ", format_p(independence_check$p.value,digits=3))

  } else if(input$color_var_surv_CTx_2 == "treat_option"){
    group_var <- "EP_option"
    plot_title <- "Expert panel recommended treatment option existed or not"

  } else if(input$color_var_surv_CTx_2 == "treat_option_pos"){
    plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>% dplyr::filter(EP_option == 1)
    group_var <- "treat_group_2"
    plot_title <- "Only patients with expert panel recommended treatment option"

  } else if(input$color_var_surv_CTx_2 == "treat_group"){
    group_var <- "treat_group"
    plot_title <- "EP option and treatment"

  } else if(input$color_var_surv_CTx_2 == "treat_group_2"){
    group_var <- "treat_group_2"
    plot_title <- "EP option and treatment"

  } else if(input$color_var_surv_CTx_2 == "treat_group_3"){
    group_var <- "treat_group_3"
    plot_title <- "EP option and treatment"

  } else if(input$color_var_surv_CTx_2 == "treat_group_4"){
    plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>%
      dplyr::mutate(CTx_lines_before_CGP = case_when(
        CTx_lines_before_CGP %in% c("0", "1", "2", "3") ~ CTx_lines_before_CGP,
        TRUE ~ "4~"
      ))
    group_var <- "CTx_lines_before_CGP"
    plot_title <- "CTx lines before CGP"

  } else if(input$color_var_surv_CTx_2 == "treat_group_5"){
    group_var <- "pre_CGP_best_RECIST"
    plot_title <- "Best CTx effect before CGP"

  } else if(input$color_var_surv_CTx_2 == "treat_group_6"){
    plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>% dplyr::filter(treat_group_2 != "No treatment done")
    group_var <- "treat_group_2"
    plot_title <- "Recommended vs not-recommended treatment"

  } else if(input$color_var_surv_CTx_2 == "treat_group_7"){
    plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>% dplyr::mutate(PS = `症例.背景情報.ECOG.PS.名称.`)
    group_var <- "PS"
    plot_title <- "ECOG Performance status"

  } else if(input$color_var_surv_CTx_2 == "treat_group_8"){
    plot_data <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>%
      dplyr::filter(症例.背景情報.ECOG.PS.名称. != "Unknown") %>%
      dplyr::mutate(PS = case_when(
        症例.背景情報.ECOG.PS.名称. == "0" ~ "0",
        症例.背景情報.ECOG.PS.名称. == "1" ~ "1",
        TRUE ~ "2_4"
      ))
    group_var <- "PS"
    plot_title <- "ECOG Performance status (unknown patients excluded)"

  } else if(input$color_var_surv_CTx_2 == "treat_group_9"){
    Data_survival_tmp <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>%
      dplyr::filter(Cancers %in% names(sort(table(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive$Cancers), decreasing = TRUE))[1:min(7, length(unique(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive$Cancers)))])
    Data_survival_tmp2 <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive %>%
      dplyr::mutate(Cancers = " ALL")
    plot_data <- rbind(Data_survival_tmp, Data_survival_tmp2)
    group_var <- "Cancers"
    plot_title <- "Diagnosis"
  }

  # 最後に一度だけ関数を呼び出す
  survival_compare_and_plot_CTx(
    data = plot_data,
    time_var1 = "time_pre",
    time_var2 = "time_all",
    status_var = "censor",
    group_var = group_var,
    plot_title = plot_title,
    adjustment = OUTPUT_DATA$figure_surv_CTx_adjustment,
    color_var_surv_CTx_1 = input$color_var_surv_CTx_1
  )
})

output$figure_survival_CTx_2_data_2_raw = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive, "FELIS downloaded raw data in Custom survival analysis tab")
})

output$figure_survival_CTx_interactive_1 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target,
      OUTPUT_DATA$figure_surv_CTx_Data_drug,
      input$color_var_surv_CTx_1)

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target
  Data_drug = OUTPUT_DATA$figure_surv_CTx_Data_drug

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

    # 遺伝子フィルタ（P_1: 必須変異1）
    p1_input <- input[[paste0(input_prefix, "P_1")]]
    if(!all(is.null(p1_input))) {
      IDs <- intersect(IDs, (Data_MAF_target %>%
                               dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
    }

    # 遺伝子フィルタ（P_2: 必須変異2）
    p2_input <- input[[paste0(input_prefix, "P_2")]]
    if(!all(is.null(p2_input))) {
      IDs <- intersect(IDs, (Data_MAF_target %>%
                               dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
    }

    # 遺伝子除外（W: 除外変異）
    w_input <- input[[paste0(input_prefix, "W")]]
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
      H = "Cancers"
    )

    for(filter_key in names(clinical_filters)) {
      filter_input <- input[[paste0(input_prefix, filter_key)]]
      if(!all(is.null(filter_input))) {
        column_name <- clinical_filters[[filter_key]]
        filter_expr <- paste0(column_name, " %in% filter_input")
        IDs <- intersect(IDs, (Data_survival_interactive %>%
                                 dplyr::filter(!!parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
      }
    }

    # 薬剤フィルタ（D）
    d_input <- input[[paste0(input_prefix, "D")]]
    if(!all(is.null(d_input))) {
      IDs <- intersect(IDs, (Data_drug %>%
                               dplyr::filter(Drug %in% d_input))$ID)
    }

    return(IDs)
  }

  # Group1とGroup2のIDを取得
  ID_1 <- extract_group_ids(1)
  ID_2 <- extract_group_ids(2)

  # データ作成
  Data_survival_1 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>%
    dplyr::mutate(Group = 1)

  Data_survival_2 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>%
    dplyr::mutate(Group = 2)

  Data_survival <- rbind(Data_survival_1, Data_survival_2)

  # サバイバル解析実行
  survival_compare_and_plot_CTx(
    data = Data_survival,
    time_var1 = "time_pre",
    time_var2 = "time_all",
    status_var = "censor",
    group_var = "Group",
    plot_title = "Survival analisys based on variants and drugs",
    adjustment = OUTPUT_DATA$figure_surv_CTx_adjustment,
    color_var_surv_CTx_1 = input$color_var_surv_CTx_1
  )
})
