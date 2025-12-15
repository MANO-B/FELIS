drug_table_logic <- function() {
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
                      YoungOld,
                      Lymph_met,
                      Brain_met,
                      Lung_met,
                      Bone_met,
                      Liver_met,
                      Other_met,
                      EP_option,
                      EP_treat,
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
                      time_2L_enroll,
                      PD_L1
        )
      if(!"LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
        Data_case_target = Data_case_target %>% dplyr::select(
          -PD_L1
        )
      }
      
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
      
      
      Data_MAF_target = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
      if(input$patho == "Only pathogenic muts"){
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Evidence_level == "F")
      } else {
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(!Hugo_Symbol %in% c("TMB", "MSI") | Evidence_level == "F")
        
      }
      Data_MAF_target_tmp = Data_MAF_target  %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        Hugo_Symbol,
                        amino.acid.change,
                        .keep_all = TRUE)
      Data_MAF_target = Data_MAF_target  %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        Hugo_Symbol,
                        .keep_all = TRUE)
      if(nrow(Data_case_target)>0){
        Data_drug = Data_drug_raw_rename()
        Data_drug = Data_drug %>%
          dplyr::filter(ID %in%
                          Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
        Data_drug$AE = ifelse(Data_drug$Adverse_effect == "AE (-)", 0, 1)
        Data_drug$治療ライン = as.character(Data_drug$CTx_line)
        Drug_summary =  Data_drug %>%
          dplyr::select(
            ID,
            Drug,
            CTx_line,
            最良総合効果,
            終了理由,
            Adverse_effect,
            Overall_AE_grage,
            all_of(names(classification_rules)),
            TTD,
            censor,
            ToT,
            ToT_censor,
            TTF,
            TTF_censor,
            TTAE
          ) %>%
          dplyr::distinct()
        colnames(Drug_summary)=
          c("ID", "Drugs", "CTx line", "RECIST", "Reason to finish treatment",
            "Any adverse effect (G3-G5)",
            "Max AE grade",
            names(classification_rules),
            "Time to death", "Time to death, finished or censored",
            "Time on treatment", "ToT, finished or censored",
            "Time to treatment failure", "TTF, finished or censored", "Time to adverse effect")
        Drug_summary$`Reason to finish treatment`[is.na(Drug_summary$`Reason to finish treatment`)] = "Unknown"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "その他理由で中止"] = "Other reason"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "不明"] = "Unknown"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "副作用等で中止"] = "Side effect"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "本人希望により中止"] = "Patient will"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "無効中止"] = "Not effective"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "死亡中止"] = "Death"
        Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "計画通り終了"] = "As planned"
        Drug_summary$`Time to death, finished or censored` = as.character(Drug_summary$`Time to death, finished or censored`)
        Drug_summary$`Time to death, finished or censored`[is.na(Drug_summary$`Time to death, finished or censored`)] = "No TTD Data"
        Drug_summary$`Time to death, finished or censored`[Drug_summary$`Time to death, finished or censored` == "1"] = "Finished"
        Drug_summary$`Time to death, finished or censored`[Drug_summary$`Time to death, finished or censored` == "0"] = "Censored"
        Drug_summary$`ToT, finished or censored` = as.character(Drug_summary$`ToT, finished or censored`)
        Drug_summary$`ToT, finished or censored`[is.na(Drug_summary$`ToT, finished or censored`)] = "No ToT Data"
        Drug_summary$`ToT, finished or censored`[Drug_summary$`ToT, finished or censored` == "1"] = "Finished"
        Drug_summary$`ToT, finished or censored`[Drug_summary$`ToT, finished or censored` == "0"] = "Censored"
        Drug_summary$`TTF, finished or censored` = as.character(Drug_summary$`TTF, finished or censored`)
        Drug_summary$`TTF, finished or censored`[is.na(Drug_summary$`TTF, finished or censored`)] = "No TTF Data"
        Drug_summary$`TTF, finished or censored`[Drug_summary$`TTF, finished or censored` == "1"] = "Finished"
        Drug_summary$`TTF, finished or censored`[Drug_summary$`TTF, finished or censored` == "0"] = "Censored"
        if(!"LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
          Data_case_age_sex = Data_case_target %>%
            dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                          症例.基本情報.年齢,
                          YoungOld,
                          症例.基本情報.性別.名称.,
                          症例.背景情報.喫煙歴有無.名称.,
                          症例.背景情報.アルコール多飲有無.名称.,
                          症例.基本情報.がん種.OncoTree..名称.
            ) %>% distinct()
          colnames(Data_case_age_sex) = c(
            "ID",
            "Age",
            paste0("Age <=", input$mid_age, " or older"),
            "Sex",
            "Smoking history",
            "Heavy alcohol",
            "Diagnosis"
          )
        } else{
          Data_case_age_sex = Data_case_target %>%
            dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                          症例.基本情報.年齢,
                          YoungOld,
                          症例.基本情報.性別.名称.,
                          症例.背景情報.喫煙歴有無.名称.,
                          症例.背景情報.アルコール多飲有無.名称.,
                          PD_L1,
                          症例.基本情報.がん種.OncoTree..名称.
            ) %>% distinct()
          colnames(Data_case_age_sex) = c(
            "ID",
            "Age",
            paste0("Age <=", input$mid_age, " or older"),
            "Sex",
            "Smoking history",
            "Heavy alcohol",
            "PD-L1 (lung only)",
            "Diagnosis"
          )
        }
        Data_cluster_ID_list = Data_cluster_ID() %>%
          dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster)
        Data_case_age_sex = left_join(Data_case_age_sex,
                                     Data_cluster_ID_list,
                                     by = c("ID" = "C.CAT調査結果.基本項目.ハッシュID"))
        Data_case_age_sex$cluster = as.factor(Data_case_age_sex$cluster)
        
        Drug_summary = Drug_summary %>%
          dplyr::left_join(Data_case_age_sex, by = "ID")
        
        remove_cols <- names(Drug_summary) %in% names(classification_rules) & 
          sapply(Drug_summary, is.logical)
        Drug_summary <- Drug_summary[, !remove_cols, with = FALSE] %>%
          dplyr::select(where(~ any(!is.na(.))))
        # Drug_summaryの各Drugの行数をカウントし、閾値未満を"Others"に変換
        Drug_summary <- Drug_summary %>%
          add_count(Drugs, name = "drug_count") %>%
          mutate(Drugs = ifelse(drug_count < ifelse(!is.null(input$minimum_courses), input$minimum_courses, 10), "Others", Drugs)) %>%
          select(-drug_count)
        OUTPUT_DATA$drug_table_Drug_summary = Drug_summary
        OUTPUT_DATA$drug_table_Data_case_target = Data_case_target
        Drug_summary_renamed = Drug_summary
        if(!is.null(input$gene_group_1)){
          if(!is.null(input$gene_group_2)){
            Data_MAF_target_tmp = Data_MAF_target_tmp %>%
              dplyr::filter(Tumor_Sample_Barcode %in%
                              Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
              dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)
            
            ID_special_gene_mutation_1 = (Data_MAF_target_tmp %>%
                                            dplyr::filter(Hugo_Symbol %in% input$gene_group_1 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_1, "_"),collapse ="|")) &
                                                            amino.acid.change %in% input$special_gene_mutation_1 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_1, "_"),collapse ="|")) &
                                                            amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
            ID_special_gene_mutation_2 = (Data_MAF_target_tmp %>%
                                            dplyr::filter(Hugo_Symbol %in% input$gene_group_2 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_2, "_"),collapse ="|")) &
                                                            amino.acid.change %in% input$special_gene_mutation_1 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_2, "_"),collapse ="|")) &
                                                            amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
            Drug_summary_line_renamed = Drug_summary_renamed %>% dplyr::mutate(
              separation_value = case_when(
                ID %in% ID_special_gene_mutation_1 &
                  ID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut, ", paste0(paste(input$gene_group_2, collapse="/"), " mut")),
                ID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut, ", paste0(paste(input$gene_group_2, collapse="/"), " WT")),
                ID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$gene_group_1, collapse="/"), " WT, ", paste0(paste(input$gene_group_2, collapse="/"), " mut")),
                TRUE ~ paste0("No mut in ",  paste(input$gene_group_1, collapse="/"), "/", paste(input$gene_group_2, collapse="/"))
              )
            )
          } else {
            Data_MAF_target_tmp = Data_MAF_target_tmp %>%
              dplyr::filter(Tumor_Sample_Barcode %in%
                              Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
              dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)
            
            ID_special_gene_mutation_1 = (Data_MAF_target_tmp %>%
                                            dplyr::filter(Hugo_Symbol %in% input$gene_group_1 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_1, "_"),collapse ="|")) &
                                                            amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
            Drug_summary_line_renamed = Drug_summary_renamed %>% dplyr::mutate(
              separation_value = case_when(
                ID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut"),
                TRUE ~ paste0("No mut in ", paste(input$gene_group_1, collapse="/"))
              )
            )
          }
        } else if(!input$special_gene == ""){
          Data_MAF_target_tmp = Data_MAF_target_tmp %>%
            dplyr::filter(Tumor_Sample_Barcode %in%
                            Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
            dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)
          
          ID_special_gene = (Data_MAF_target_tmp %>%
                               dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_"))))$Tumor_Sample_Barcode
          ID_special_gene_mutation_1 = (Data_MAF_target_tmp %>%
                                          dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                          amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
          ID_special_gene_mutation_2 = (Data_MAF_target_tmp %>%
                                          dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                          amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
          Drug_summary_line_renamed = Drug_summary_renamed %>% dplyr::mutate(
            separation_value = case_when(
              ID %in% ID_special_gene_mutation_1 ~ paste0(input$special_gene_mutation_1_name, " in ", input$special_gene),
              ID %in% ID_special_gene_mutation_2 ~ paste0(input$special_gene_mutation_2_name, " in ", input$special_gene),
              ID %in% ID_special_gene ~ paste0("Other mut in ", input$special_gene),
              TRUE ~ paste0("No mut in ", input$special_gene)
            )
          )
        } else {
          Drug_summary_line_renamed = Drug_summary_renamed %>% dplyr::mutate(
            separation_value = "No mutation pattern")
        }
        # if(!is.null(input$target_line)){
        #   Drug_summary_line_renamed = Drug_summary_line_renamed %>% dplyr::filter(`CTx line` %in% input$target_line)
        # }
        incProgress(1 / 3)
        
        Drug_summary_diagnosis_renamed = Drug_summary_line_renamed
        ID_diagnosis_list = Data_case_target %>% dplyr::select(C.CAT調査結果.基本項目.ハッシュID,症例.基本情報.がん種.OncoTree.) %>% dplyr::distinct()
        Drug_summary_diagnosis_renamed$diagnosis = unlist(lapply(list(Drug_summary_diagnosis_renamed$ID), function(x) {
          as.vector(ID_diagnosis_list$症例.基本情報.がん種.OncoTree.[match(x, ID_diagnosis_list$C.CAT調査結果.基本項目.ハッシュID)])}))
        mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
        if(length(mut_gene_)==0){
          mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
          ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
          Drug_summary_mutation = Drug_summary_diagnosis_renamed %>% dplyr::mutate(
            mutation = case_when(
              ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
              TRUE ~ paste0(mut_gene_, " mut(-)")
            )
          )
        } else{
          if(is.null(input$gene_group_2)){
            ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
            Drug_summary_mutation = Drug_summary_diagnosis_renamed %>% dplyr::mutate(
              mutation = case_when(
                ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
                TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
              )
            )
          } else{
            ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
            ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
            Drug_summary_mutation = Drug_summary_diagnosis_renamed %>% dplyr::mutate(
              mutation = case_when(
                ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
              )
            )
          }
        }
        OUTPUT_DATA$drug_table_Drug_summary_line_renamed = Drug_summary_line_renamed
        OUTPUT_DATA$drug_table_Drug_summary_diagnosis_renamed = Drug_summary_diagnosis_renamed
        OUTPUT_DATA$drug_table_Drug_summary_mutation = Drug_summary_mutation
      }
      incProgress(1 / 3)
    })
  })
}

# メイン関数：集計テーブル作成
create_summary_table <- function(data, group_var = NULL, caption, is_detailed = TRUE) {
  
  # data.table形式の場合はdata.frameに変換
  if ("data.table" %in% class(data)) {
    data <- as.data.frame(data)
  }
  
  # 集計表の作成
  if (!is_detailed) {
    # 簡略表示バージョン
    tbl_output =
      tbl_summary(
        data,
        by = group_var,
        statistic = list(
          all_continuous() ~ "{mean} ({sd})",
          all_categorical() ~ "{n}"
        ),
        type = list(
          all_continuous() ~ "continuous2",
          all_dichotomous() ~ "categorical"
        ),
        digits = list(
          all_continuous() ~ list(mean = 1, sd = 1),
          all_categorical() ~ 0
        ),
        missing = "ifany"
      ) %>%
      modify_caption(caption) %>%
      add_n() %>%
      bold_labels()
    if (!is.null(group_var)) {
      tbl_output <- tbl_output %>%
        add_p(
          test = list(
            all_continuous() ~ "kruskal.test",
            all_categorical() ~ "chisq.test.no.correct",
            all_dichotomous() ~ "fisher.test"
          ),
          pvalue_fun = ~ style_pvalue(., digits = 3)
        ) %>%
        bold_p(t = 0.05) %>%           # p < 0.05を太字
        add_significance_stars()
    }
    tbl_output = tbl_output %>% as_gt()
  } else {
    # 詳細表示バージョン
    tbl_output =
      tbl_summary(
        data,
        by = group_var,
        statistic = list(
          all_continuous() ~ c("{N_nonmiss}",
                               "{mean} ({sd})",
                               "{median} ({p25}, {p75})", 
                               "{min}, {max}"),
          all_categorical() ~ "{n} ({p}%)"
        ),
        type = list(
          all_continuous() ~ "continuous2",
          all_dichotomous() ~ "categorical"
        ),
        digits = list(
          all_continuous() ~ list(N_nonmiss = 0, mean = 1, sd = 1, median = 1,
                                  p25 = 1, p75 = 1, min = 1, max = 1),
          all_categorical() ~ c(0, 1)
        ),
        missing = "ifany"
      ) %>%
      modify_table_body(
        ~ .x %>%
          mutate(across(starts_with("stat_"),
                        ~ str_replace_all(.x, "\\(0\\.0%\\)", "(<0.1%)")))
      ) %>%
      modify_caption(caption) %>%
      add_n() %>%
      bold_labels()
    if (!is.null(group_var)) {
      tbl_output <- tbl_output %>%
        add_p(
          test = list(
            all_continuous() ~ "kruskal.test",
            all_categorical() ~ "chisq.test.no.correct",
            all_dichotomous() ~ "fisher.test"
          ),
          pvalue_fun = ~ style_pvalue(., digits = 3)
        ) %>%
        bold_p(t = 0.05) %>%           # p < 0.05を太字
        add_significance_stars()
    }
    tbl_output = tbl_output %>% as_gt()
  }
  return(tbl_output)
}

# データ準備関数
prepare_table_data <- function(base_data, table_var, target_line, mid_age = NULL) {
  # データをdata.frameに変換（data.table競合を回避）
  if ("data.table" %in% class(base_data)) {
    base_data <- as.data.frame(base_data)
  }
  
  # 基本フィルタリング
  filtered_data <- base_data %>%
    dplyr::filter(`CTx line` %in% target_line) %>%
    dplyr::select(-ID)
  
  # 特殊なデータセット処理
  if (table_var %in% c("mutations", "diagnosis", "gene")) {
    # これらの場合は異なるデータセットを使用
    return(get_special_dataset(table_var))
  }
  
  return(filtered_data)
}

# 特殊データセット取得関数
get_special_dataset <- function(table_var) {
  result <- switch(table_var,
                   "mutations" = OUTPUT_DATA$drug_table_Drug_summary_line_renamed %>% dplyr::select(-ID),
                   "diagnosis" = OUTPUT_DATA$drug_table_Drug_summary_diagnosis_renamed %>% dplyr::select(-ID, -separation_value),
                   "gene" = OUTPUT_DATA$drug_table_Drug_summary_mutation %>% dplyr::select(-ID, -separation_value, -diagnosis)
  )
  
  # data.table形式の場合はdata.frameに変換
  if ("data.table" %in% class(result)) {
    result <- as.data.frame(result)
  }
  
  return(result)
}

# グループ変数とキャプション設定関数
get_group_and_caption <- function(table_var, mid_age = NULL, target_line = NULL) {
  switch(table_var,
         "All" = list(
           group_var = NULL,
           caption = "**Palliative CTx with treatment duration information**"
         ),
         "Adverse_effect" = list(
           group_var = "Any adverse effect (G3-G5)",
           caption = "**Palliative CTx with treatment duration information, by adverse effect**"
         ),
         "lines" = list(
           group_var = "CTx line",
           caption = "**Palliative CTx with treatment duration information, by line**"
         ),
         "RECIST" = list(
           group_var = "RECIST",
           caption = "**Palliative CTx with treatment duration information, by treatment response**"
         ),
         "age" = list(
           group_var = paste0("Age <=", mid_age, " or older"),
           caption = paste0("**Palliative CTx with treatment duration information, younger:age <=", mid_age, "**")
         ),
         "sex" = list(
           group_var = "Sex",
           caption = "**Palliative CTx with treatment duration information, by sex**"
         ),
         "PD-L1" = list(
           group_var = "PD-L1 (lung only)",
           caption = "**Palliative CTx with treatment duration information, by PD-L1 status**"
         ),
         "mutations" = list(
           group_var = "separation_value",
           caption = "**Palliative CTx, selected lines, by mutation pattern**"
         ),
         "Diagnosis" = list(
           group_var = "Diagnosis",
           caption = "**Palliative CTx, selected lines, by diagnosis**"
         ),
         "cluster" = list(
           group_var = "cluster",
           caption = "**Palliative CTx, selected lines, by mutation-based cluster**"
         ),
         "gene" = list(
           group_var = "mutation",
           caption = "**Palliative CTx, selected lines, by mutation**"
         ),
         # デフォルト
         list(group_var = NULL, caption = "**Summary Table**")
  )
}


# メインのテーブル出力処理
output$table_drug_all_1 <- render_gt({
  # 必要な入力値チェック
  req(input$target_line, input$table_var_drug_1, input$drug_table_layout)
  
  # 表示モード確認（詳細/簡略）
  is_detailed = input$drug_table_layout != "No"

  # データ準備
  table_data <- prepare_table_data(
    base_data = OUTPUT_DATA$drug_table_Drug_summary,
    table_var = input$table_var_drug_1,
    target_line = input$target_line,
    mid_age = input$mid_age
  )
  
  # グループ変数とキャプション設定
  config <- get_group_and_caption(
    table_var = input$table_var_drug_1,
    mid_age = input$mid_age,
    target_line = input$target_line
  )
  
  # テーブル作成
  create_summary_table(
    data = table_data,
    group_var = config$group_var,
    caption = config$caption,
    is_detailed = is_detailed
  )
})

output$select_table_var_drug_1 = renderUI({
  choices = c("All" = "All",
              "Treatment line" = "lines",
              "Treatment effect (RECIST)" = "RECIST",
              "Mutation-based cluster" = "cluster",
              "Diagnosis" = "Diagnosis",
              "Age" = "age",
              "Sex" = "sex",
              "Adverse effect" = "Adverse_effect",
              "Mutation in a gene set" = "gene",
              "PD-L1" = "PD-L1"
  )
  if(!"LUNG" %in% OUTPUT_DATA$drug_table_Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
    choices = choices[choices != "PD-L1"]
  }
  radioButtons(
    inputId = "table_var_drug_1",
    label = "Grouped by",
    choices = choices,
    selected = "All")
})
