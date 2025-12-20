survival_CTx_analysis2_logic <- function() {
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

        if(is.null(input$color_var_surv_CTx_1)){
          Data_case_target$time_pre = Data_case_target$time_palliative_enroll
          Data_case_target$time_all = Data_case_target$time_palliative_final
          Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)
        } else if(input$color_var_surv_CTx_1 == "diagnosis"){
          Data_case_target$time_pre = Data_case_target$time_diagnosis_enroll
          Data_case_target$time_all = Data_case_target$time_diagnosis_final
          Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)
        } else if(input$color_var_surv_CTx_1 == "1L_CTx"){
          Data_case_target$time_pre = Data_case_target$time_palliative_enroll
          Data_case_target$time_all = Data_case_target$time_palliative_final
          Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)
        } else if(input$color_var_surv_CTx_1 == "2L_CTx"){
          Data_case_target$time_pre = Data_case_target$time_2L_enroll
          Data_case_target$time_all = Data_case_target$time_2L_final
          Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)
        } else {
          Data_case_target$time_pre = rep(0, nrow(Data_case_target))
          Data_case_target$time_all = Data_case_target$time_enroll_final
        }
        if(is.null(input$color_var_surv_CTx_3)){
          adjustment = TRUE
        } else if(input$color_var_surv_CTx_3 == "Yes"){
          adjustment = TRUE
        } else {
          adjustment = FALSE
        }
        OUTPUT_DATA$figure_surv_CTx_adjustment = adjustment
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
        OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive = Data_survival_interactive

        candidate_genes = sort(unique(c(Data_MAF_target$Hugo_Symbol,
                                        paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                                        paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                                        paste0(input$special_gene, "_NOS"))))
        candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
        candidate_genes = candidate_genes[!is.na(candidate_genes)]
        Top_gene = unique(names(sort(table(Data_MAF_target$Hugo_Symbol), decreasing = T)))
        Top_gene = Top_gene[Top_gene %in% candidate_genes]
        candidate_drugs = sort(unique(c(Data_drug$Drug)))
        OUTPUT_DATA$figure_surv_interactive_Top_gene = Top_gene
        OUTPUT_DATA$figure_surv_interactive_candidate_genes = candidate_genes
        OUTPUT_DATA$figure_surv_interactive_candidate_drugs = candidate_drugs[!is.na(candidate_drugs)]
        OUTPUT_DATA$figure_surv_interactive_candidate_lines = sort(unique(Data_survival_interactive$CTx_lines_before_CGP))
        OUTPUT_DATA$figure_surv_interactive_candidate_RECIST = sort(unique(Data_survival_interactive$pre_CGP_best_RECIST))
        OUTPUT_DATA$figure_surv_interactive_candidate_Age = sort(unique(Data_survival_interactive$YoungOld))
        OUTPUT_DATA$figure_surv_interactive_candidate_Sex = sort(unique(Data_survival_interactive$症例.基本情報.性別.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_Histology = sort(unique(Data_survival_interactive$Cancers))
        OUTPUT_DATA$figure_surv_interactive_candidate_cluster = sort(unique(Data_survival_interactive$cluster))
        OUTPUT_DATA$figure_surv_interactive_candidate_PS = sort(unique(Data_survival_interactive$症例.背景情報.ECOG.PS.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_meta = c('Lymph_met','Brain_met','Lung_met','Bone_met','Liver_met')


        oncogenic_genes = Data_MAF_target %>%
          dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
          dplyr::distinct() %>%
          dplyr::count(Hugo_Symbol) %>%
          dplyr::arrange(-n)
        colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
        oncogenic_genes = rbind(oncogenic_genes %>% dplyr::filter(gene_mutation %in% input$gene),
                                oncogenic_genes %>% dplyr::filter(!gene_mutation %in% input$gene))
        oncogenic_genes = oncogenic_genes[1:min(length(oncogenic_genes$all_patients), input$gene_no),]
        OUTPUT_DATA$figure_surv_CTx_oncogenic_genes = oncogenic_genes

        incProgress(1 / 13)

        Data_survival = Data_survival %>%
          dplyr::mutate(Cancer_name_event = case_when(
            censor == 1 ~ Cancers,
            TRUE ~ "Others"
          ))
        Cancer_table = sort(table(Data_survival$Cancer_name_event),decreasing = T)
        Cancer_names = names(Cancer_table[Cancer_table>=3])
        Cancer_names = Cancer_names[Cancer_names != "Others"]
        Data_survival = Data_survival %>%
          dplyr::mutate(Cancer_name = case_when(
            Cancers %in% Cancer_names ~ Cancers,
            TRUE ~ "Others"
          ))
        Cancer_names = c(Cancer_names, "Others")
        Data_survival$Cancer_name = as.factor(Data_survival$Cancer_name)
        Data_survival$Cancer_name = relevel(Data_survival$Cancer_name,
                                            ref=Cancer_names[[1]])
        Data_survival$cluster = as.character(Data_survival$cluster)
        Data_survival = Data_survival %>%
          dplyr::mutate(cluster_name_event = case_when(
            censor == 1 ~ cluster,
            TRUE ~ "Others"
          ))
        cluster_table = sort(table(Data_survival$cluster_name_event),decreasing = T)
        cluster_names = names(cluster_table[cluster_table>=3])
        cluster_names = cluster_names[cluster_names != "Others"]
        Data_survival = Data_survival %>%
          dplyr::mutate(cluster = case_when(
            cluster %in% cluster_names ~ cluster,
            TRUE ~ "Others"
          ))
        cluster_names = c(cluster_names, "Others")
        Data_survival$cluster = as.factor(Data_survival$cluster)
        Data_survival$cluster = relevel(Data_survival$cluster,
                                        ref=cluster_names[[1]])
        OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer = Data_survival %>% dplyr::filter(Cancer_name != "Others")
        OUTPUT_DATA$figure_surv_CTx_Data_survival_cluster = Data_survival %>% dplyr::filter(cluster != "Others")
        OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer_cluster = OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer %>% dplyr::filter(cluster != "Others")
      }
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}

output$figure_survival_CTx_4_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer)
  ggforest(coxph(Surv(time = time_palliative_enroll,
                      time2 = time_palliative_final, censor)~Cancer_name,
                 data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer), data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer)
})
output$figure_survival_CTx_5_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_cluster)
  ggforest(coxph(Surv(time = time_palliative_enroll,
                      time2 = time_palliative_final, censor)~cluster,
                 data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cluster), data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cluster)
})
output$figure_survival_CTx_5_1_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer_cluster)
  ggforest(coxph(Surv(time = time_palliative_enroll,
                      time2 = time_palliative_final, censor)~cluster+Cancer_name,
                 data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer_cluster), data=OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer_cluster)
})

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
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive, "FELIS downloaded raw data in Survival and clinical information tab")
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
      H = "Cancers",
      C = "cluster"
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

    # 転移部位フィルタ（M: 複数選択可能、OR条件）
    m_input <- input[[paste0(input_prefix, "M")]]
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

output$figure_survival_CTx_3_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target,
      OUTPUT_DATA$figure_surv_CTx_oncogenic_genes,
      input$color_var_surv_CTx_1)
  oncogenic_genes = OUTPUT_DATA$figure_surv_CTx_oncogenic_genes
  Data_survival = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target
  adjustment = OUTPUT_DATA$figure_surv_CTx_adjustment
  # analysis for common oncogenic mutations
  hs = list()

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

  k = 1
  withProgress(message = "Analyzing common mutated genes", {
    for(i in 1:nrow(oncogenic_genes)){
      incProgress(1 / nrow(oncogenic_genes))
      ID_mutation = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% oncogenic_genes$gene_mutation[i]))$Tumor_Sample_Barcode
      Data_survival = Data_survival %>% dplyr::mutate(
        gene_mut = case_when(
          C.CAT調査結果.基本項目.ハッシュID %in% ID_mutation ~ 1,
          TRUE ~ 0
        )
      )
      if(adjustment){
        traditional_fit = survfit(Surv(time = time_pre,
                                       time2 = time_all, censor) ~ gene_mut,
                                  data = Data_survival)
      } else {
        traditional_fit = survfit(Surv(time = time_all, censor) ~ gene_mut,
                                  data = Data_survival)
      }
      if(length(Data_survival[,1]) != traditional_fit[[1]][[1]]){
        if(adjustment){
          diff_0 = coxph(Surv(time = time_pre,
                              time2 = time_all, censor)~gene_mut,
                         data=Data_survival)
          Xlab = paste0("Time from ", input$color_var_surv_CTx_1, ", risk-set adjusted (months)")
        } else {
          diff_0 = coxph(Surv(time = time_all, censor)~gene_mut,
                         data=Data_survival)
          Xlab = paste0("Time from ", input$color_var_surv_CTx_1, ", bias not adj (months)")
        }
        mut_gene = paste(oncogenic_genes$gene_mutation[i],
                         ", top ", i, " gene", sep="")
        tmp = data.frame(summary(traditional_fit)$table)
        gene_table$positive_median[i] = round(tmp$median[2] / 365.25 * 12, digits = 1)
        gene_table$positive_LL[i] = round(tmp$X0.95LCL[2] / 365.25 * 12, digits = 1)
        gene_table$positive_UL[i] = round(tmp$X0.95UCL[2] / 365.25 * 12, digits = 1)
        gene_table$negative_median[i] = round(tmp$median[1] / 365.25 * 12, digits = 1)
        gene_table$negative_LL[i] = round(tmp$X0.95LCL[1] / 365.25 * 12, digits = 1)
        gene_table$negative_UL[i] = round(tmp$X0.95UCL[1] / 365.25 * 12, digits = 1)
        data_tmp = tidy(diff_0, exponentiate=TRUE, conf.int=TRUE)
        gene_table$diff_median[i] = round(data_tmp[[2]], digits = 2)
        gene_table$diff_LL[i] = round(data_tmp[[6]], digits = 2)
        gene_table$diff_UL[i] = round(data_tmp[[7]], digits = 2)

        hs[[k]] = surv_curv_CTx(traditional_fit, Data_survival,
                                paste0(oncogenic_genes$gene_mutation[i], " mutation and overall survival"),
                                c("Not mut", paste0(oncogenic_genes$gene_mutation[i], " mut")), diff_0, Xlab)
        k = k + 1
      } else{
        tmp = summary(traditional_fit)$table
        legends = paste0(format_p(tmp[[7]] / 365.25 * 12, digits = 1), " (", format_p(tmp[[8]] / 365.25 * 12, digits = 1), "-", format_p(tmp[[9]] / 365.25 * 12, digits = 1),")")
        gene_table$negative_median[i] = round(tmp[[7]] / 365.25 * 12, digits = 1)
        gene_table$negative_LL[i] = round(tmp[[8]] / 365.25 * 12, digits = 1)
        gene_table$negative_UL[i] = round(tmp[[9]] / 365.25 * 12, digits = 1)
      }
    }
  })
  Gene_arrange = gene_table$Gene
  gene_table$Gene = factor(gene_table$Gene, levels = Gene_arrange)
  # -Inf, Inf を置き換える新しい列を作る
  min_val <- min(gene_table$diff_LL[gene_table$diff_LL != 0], na.rm = TRUE)
  max_val <- max(gene_table$diff_UL[is.finite(gene_table$diff_UL)], na.rm = TRUE)
  gene_table$diff_LL_plot <- ifelse(gene_table$diff_LL != 0, gene_table$diff_LL, min_val)
  gene_table$diff_UL_plot <- ifelse(is.finite(gene_table$diff_UL), gene_table$diff_UL, max_val)
  gene_table$diff_median_plot <- ifelse(gene_table$diff_median >= gene_table$diff_LL_plot &
                                          gene_table$diff_median <= gene_table$diff_UL_plot, gene_table$diff_median, NA)
  withProgress(message = "Drawing survival curves...", {
    p_mid <-
      gene_table |>
      ggplot(aes(y = fct_rev(Gene))) +
      theme_classic() +  scale_x_log10() +
      geom_point(aes(x=diff_median_plot), shape=15, size=3) +
      geom_linerange(aes(xmin=diff_LL_plot, xmax=diff_UL_plot))  +
      geom_vline(xintercept = 1, linetype="dashed") +
      labs(x="Hazard ratio", y="Genes with oncogenic alteration") +
      coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                      xlim=c( min_val ^ 1.1,
                              max_val ^ 1.1)) +
      annotate("text",
               x = sqrt(max_val ^ 1.1),
               y = length(Gene_arrange) + 1,
               label = "mut=short survival") +
      annotate("text",
               x = sqrt(min_val ^ 1.1),
               y = length(Gene_arrange) + 1,
               label = "mut=long survival") +
      theme(legend.position = "none",
            axis.line.y = element_blank(),
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
          estimate_lab_1 = "Survival, mut (+)",
          estimate_lab_2 = "Survival, mut (-)",
          estimate_lab_3 = "Hazard ratio",
          patients = "Patients"
        ), gene_table)
    gene_table$Gene = factor(gene_table$Gene, levels = Gene_arrange)
    OUTPUT_DATA$figure_surv_CTx_gene_table = gene_table

    p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
    p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                 fontface = ifelse(gene_table$Gene == "Gene", "bold", "bold.italic"))
    p_left <- p_left + geom_text(aes(x = 1.8, label = estimate_lab_1), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_1 == "Survival, mut (+)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 3.3, label = estimate_lab_2), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_2 == "Survival, mut (-)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 4.8, label = estimate_lab_3), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_3 == "Hazard ratio", "bold", "plain"))
    p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 7.5))

    # right side of plot - pvalues
    p_right <- gene_table  |> ggplot() +
      geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
      theme_void()

    for(kk in k:32){
      hs[[kk]] = ggsurvplot_empty()
      incProgress(1 / 32)
    }

    # final plot arrangement
    layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 7.5),
                patchwork::area(t = 1, l = 7.5, b = 30, r = 12.5),
                patchwork::area(t = 0, l = 13.5, b = 30, r = 14))
    OUTPUT_DATA$figure_surv_CTx_forest_plot = p_left + p_mid + p_right + plot_layout(design = layout)
  })
  arrange_ggsurvplots(hs,print=TRUE,ncol=4,nrow=8,surv.plot.height = 0.8,risk.table.height = 0.2)
})

output$figure_survival_CTx_2_2 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_forest_plot)
  withProgress(message = "Drawing forest plot...", {
    OUTPUT_DATA$figure_surv_CTx_forest_plot
  })
})
output$figure_survival_CTx_2_data_2 = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$figure_surv_CTx_gene_table)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_CTx_gene_table, "FELIS downloaded raw data in Frequent variants and survival tab")
})
