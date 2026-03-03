summary_base_logic <- function() {
  analysis_env <- new.env()
  clear_reactive_data(OUTPUT_DATA)
  gc()
  with(analysis_env, {
    withProgress(message = sample(nietzsche)[1], {
      Data_case_target = Data_case()
      if(input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed", "Only cases without mutations in the gene set are analyzed")){
        Data_case_target = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                          Data_report()$Tumor_Sample_Barcode)
      }
      incProgress(1 / 8)
      Data_summary = Data_case_target %>%
        dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID,.keep_all=TRUE) %>%
        dplyr::select(
          C.CAT調査結果.基本項目.ハッシュID,
          症例.基本情報.がん種.OncoTree.,
          症例.基本情報.がん種.OncoTree..名称.,
          症例.基本情報.がん種.OncoTree.LEVEL1.,
          症例.基本情報.性別.名称.,
          症例.基本情報.年齢,
          症例.背景情報.診断日,
          症例.管理情報.登録日,
          症例.背景情報.家族歴有無.名称.,
          症例.背景情報.喫煙歴有無.名称.,
          症例.背景情報.喫煙年数,
          症例.背景情報.アルコール多飲有無.名称.,
          症例.背景情報.ECOG.PS.名称.,
          症例.背景情報.重複がん有無.異なる臓器..名称.,
          症例.背景情報.多発がん有無.同一臓器..名称.,
          症例.検体情報.パネル.名称.,
          # 症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称.,
          # 症例.EP後レジメン情報.提示された治療薬を投与した.名称.,
          症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称.,
          症例.がん種情報.肺.EGFR.名称.,
          症例.がん種情報.肺.ALK融合.名称.,
          症例.がん種情報.肺.ROS1.名称.,
          症例.がん種情報.肺.BRAF.V600..名称.,
          症例.がん種情報.肺.MET遺伝子エクソン14スキッピング変異.名称.,
          症例.がん種情報.肺.KRAS.G12C遺伝子変異.名称.,
          症例.がん種情報.肺.RET融合遺伝子.名称.,
          PD_L1,
          症例.がん種情報.乳.HER2.IHC..名称.,
          症例.がん種情報.乳.HER2.FISH..名称.,
          症例.がん種情報.乳.ER.名称.,
          症例.がん種情報.乳.PgR.名称.,
          症例.がん種情報.固形がん.マイクロサテライト不安定性.名称.,
          症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称.,
          症例.がん種情報.食道.胃.腸.HER2.名称.,
          症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称.,
          症例.がん種情報.登録時転移有無.名称.,
          症例.背景情報.初回治療前のステージ分類.名称.,
          Lymph_met,
          Brain_met,
          Lung_met,
          Bone_met,
          Liver_met,
          Other_met,
          Not_brain_bone_liver_met,
          CTx_lines_before_CGP,
          EP_treat,
          EP_option,
          YoungOld)
      incProgress(1 / 8)
      if(!"BOWEL" %in% Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1. &
         !"STOMACH" %in% Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1.){
        Data_summary = Data_summary %>% dplyr::select(
          -症例.がん種情報.食道.胃.腸.HER2.名称.,
          -症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称.
        )
      }
      if(!"LUNG" %in% Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1.){
        Data_summary = Data_summary %>% dplyr::select(
          -症例.がん種情報.肺.EGFR.名称.,
          -症例.がん種情報.肺.ALK融合.名称.,
          -症例.がん種情報.肺.ROS1.名称.,
          -症例.がん種情報.肺.BRAF.V600..名称.,
          -症例.がん種情報.肺.MET遺伝子エクソン14スキッピング変異.名称.,
          -症例.がん種情報.肺.KRAS.G12C遺伝子変異.名称.,
          -症例.がん種情報.肺.RET融合遺伝子.名称.,
          -PD_L1
        )
      }
      if(!"BREAST" %in% Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1.){
        Data_summary = Data_summary %>% dplyr::select(
          -症例.がん種情報.乳.HER2.IHC..名称.,
          -症例.がん種情報.乳.HER2.FISH..名称.,
          -症例.がん種情報.乳.ER.名称.,
          -症例.がん種情報.乳.PgR.名称.
        )
      }
      if(length(unique(Data_summary$症例.基本情報.がん種.OncoTree.)) >
         ifelse(is.null(input$table_summary_no), 20,
                input$table_summary_no)){
        Data_summary$症例.基本情報.がん種.OncoTree. = Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1.
        Data_summary$症例.基本情報.がん種.OncoTree..名称. = Data_summary$症例.基本情報.がん種.OncoTree.LEVEL1.
      }

      Data_summary = Data_summary %>% dplyr::select(
        -症例.基本情報.がん種.OncoTree.LEVEL1.,
      )
      Data_Best_Evidence_Level = Data_report() %>%
        dplyr::filter(
          !str_detect(Hugo_Symbol, ",") &
            Hugo_Symbol != "" &
            Evidence_level %in% c("A","B","C","D","E") &
            Variant_Classification != "expression"
        ) %>%
        dplyr::arrange(Evidence_level) %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        .keep_all = TRUE) %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_summary$C.CAT調査結果.基本項目.ハッシュID) %>%
        dplyr::select(Tumor_Sample_Barcode, Evidence_level)
      colnames(Data_Best_Evidence_Level) = c("C.CAT調査結果.基本項目.ハッシュID", "Best Evidence Level")
      Data_summary = left_join(Data_summary, Data_Best_Evidence_Level,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_summary$`Best Evidence Level`[is.na(Data_summary$`Best Evidence Level`)] = "None"

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
      incProgress(1 / 8)

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
      incProgress(1 / 8)
      Data_MAF_target = Data_MAF_target %>%
        dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change, .keep_all = T) %>%
        dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)
      incProgress(1 / 8)
      if(!is.null(input$table_summary_1_gene_group_1)){
        if(!is.null(input$table_summary_1_gene_group_2)){
          ID_special_gene_mutation_1 = (Data_MAF_target %>%
                                          dplyr::filter(Hugo_Symbol %in% input$table_summary_1_gene_group_1 | str_detect(Hugo_Symbol, paste(paste0(input$table_summary_1_gene_group_1, "_"),collapse ="|")) &
                                                          amino.acid.change %in% input$special_gene_mutation_1 | str_detect(Hugo_Symbol, paste(paste0(input$table_summary_1_gene_group_1, "_"),collapse ="|")) &
                                                          amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
          ID_special_gene_mutation_2 = (Data_MAF_target %>%
                                          dplyr::filter(Hugo_Symbol %in% input$table_summary_1_gene_group_2 | str_detect(Hugo_Symbol, paste(paste0(input$table_summary_1_gene_group_2, "_"),collapse ="|")) &
                                                          amino.acid.change %in% input$special_gene_mutation_1 | str_detect(Hugo_Symbol, paste(paste0(input$table_summary_1_gene_group_2, "_"),collapse ="|")) &
                                                          amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
          Data_summary = Data_summary %>% dplyr::mutate(
            separation_value = case_when(
              C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_1 &
                C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$table_summary_1_gene_group_1, collapse="/"), " mut, ", paste0(paste(input$table_summary_1_gene_group_2, collapse="/"), " mut")),
              C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$table_summary_1_gene_group_1, collapse="/"), " mut, ", paste0(paste(input$table_summary_1_gene_group_2, collapse="/"), " WT")),
              C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$table_summary_1_gene_group_1, collapse="/"), " WT, ", paste0(paste(input$table_summary_1_gene_group_2, collapse="/"), " mut")),
              TRUE ~ paste0("No mut in ",  paste(input$table_summary_1_gene_group_1, collapse="/"), "/", paste(input$table_summary_1_gene_group_2, collapse="/"))
            )
          )
        } else {
          ID_special_gene_mutation_1 = (Data_MAF_target %>%
                                          dplyr::filter(Hugo_Symbol %in% input$table_summary_1_gene_group_1 | str_detect(Hugo_Symbol, paste(paste0(input$table_summary_1_gene_group_1, "_"),collapse ="|")) &
                                                          amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
          Data_summary = Data_summary %>% dplyr::mutate(
            separation_value = case_when(
              C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$table_summary_1_gene_group_1, collapse="/"), " mut"),
              TRUE ~ paste0("No mut in ", paste(input$table_summary_1_gene_group_1, collapse="/"))
            )
          )
        }
      } else if(!input$special_gene == ""){
        ID_special_gene = (Data_MAF_target %>%
                             dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_"))))$Tumor_Sample_Barcode
        ID_special_gene_mutation_1 = (Data_MAF_target %>%
                                        dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                        amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
        ID_special_gene_mutation_2 = (Data_MAF_target %>%
                                        dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                        amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
        Data_summary = Data_summary %>% dplyr::mutate(
          separation_value = case_when(
            C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_1 ~ paste0(input$special_gene_mutation_1_name, " in ", input$special_gene),
            C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene_mutation_2 ~ paste0(input$special_gene_mutation_2_name, " in ", input$special_gene),
            C.CAT調査結果.基本項目.ハッシュID %in% ID_special_gene ~ paste0("Other mut in ", input$special_gene),
            TRUE ~ paste0("No mut in ", input$special_gene)
          )
        )
      } else {
        Data_summary = Data_summary %>% dplyr::mutate(
          separation_value = "No mutation pattern")
      }
      incProgress(1 / 8)
      Data_summary_ID = Data_summary$C.CAT調査結果.基本項目.ハッシュID
      Data_summary = Data_summary %>% dplyr::select(
        -C.CAT調査結果.基本項目.ハッシュID)

      Data_summary$Lymph_met = as.character(Data_summary$Lymph_met)
      Data_summary$Brain_met = as.character(Data_summary$Brain_met)
      Data_summary$Lung_met = as.character(Data_summary$Lung_met)
      Data_summary$Bone_met = as.character(Data_summary$Bone_met)
      Data_summary$Liver_met = as.character(Data_summary$Liver_met)
      Data_summary$Other_met = as.character(Data_summary$Other_met)
      Data_summary$Not_brain_bone_liver_met = as.character(Data_summary$Not_brain_bone_liver_met)

      Data_summary$症例.基本情報.年齢.診断時 =
        Data_summary$症例.基本情報.年齢 -
        ceiling(as.integer(as.Date(Data_summary$症例.管理情報.登録日) - as.Date(Data_summary$症例.背景情報.診断日))/365.25)
      Data_summary$症例.基本情報.年齢.診断時[Data_summary$症例.基本情報.年齢.診断時 < 0] = NA_integer_
      Data_summary = Data_summary %>% dplyr::select(
        -症例.背景情報.診断日,
        -症例.管理情報.登録日
      )

      Data_summary = Data_summary %>%
        dplyr::mutate(
          YoungOld = case_when(
            YoungOld == "Younger" ~ paste0("Younger than ", input$mid_age + 1),
            TRUE ~ paste0(input$mid_age + 1, " and older")
          ),
          EP_option = as.character(EP_option),
          EP_treat = as.character(EP_treat)
        )
      translation_map <- list(
        "Family cancer history" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
        "Double cancer in a different organ" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
        "Multiple tumor nodules in the organ" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
        "Treatment recommended by expert panel" = c("1" = "Yes", "0" = "No"),
        "Indicated treatment administered" = c("1" = "Yes", "0" = "No"),
        "Any distant metastasis" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
        "EGFR muts in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "ALK fusion in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "ROS1 fusion in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "BRAF V600E in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "MET ex14 skipping in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "KRAS G12C in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "RET fusion in lung" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "HER2 IHC in breast" = c("陽性（3+）" = "Positive", "陰性" = "Negative", "陰性（1+）" = "Negative", "境界域（2+）" = "Intermediate", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "HER2 FISH in breast" = c("陽性" = "Positive", "陰性" = "Negative", "equivocal" = "Equivocal", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "ER IHC in breast" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "PgR IHC in breast" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "MSI PCR" = c("陽性" = "Positive", "陰性" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "MMR IHC" = c("dMMR(欠損)" = "dMMR", "pMMR(正常)" = "pMMR", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "HER2 IHC in stomach/bowel" = c("陽性（3+）" = "Positive", "境界域（2+）" = "Intermediate", "陰性" = "Negative", "陰性（1+）" = "Negative", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "HER2 FISH in stomach/bowel" = c("陽性" = "Positive", "陰性" = "Negative", "equivocal" = "Equivocal", "不明or未検査" = "Unknown", "判定不能" = "Undetermined"),
        "Reason for no indicated CTx" = c(
          "EP提示の治験等に参加できなかった" = "Could not participate in EP-proposed clinical trials",
          "不明" = "Unknown",
          "主治医の主に臨床的な判断" = "Clinical judgment by the attending physician",
          "患者が化学療法を希望しなかった" = "Patient declined chemotherapy",
          "患者が治験等を希望したが、適格・除外基準や登録期間外のため参加できなかった" = "Patient wished to join a trial but was ineligible or outside registration period",
          "患者の全身状態不良により化学療法ができなかった" = "Chemotherapy not possible due to poor general condition",
          "患者の経済的事情により化学療法ができなかった" = "Chemotherapy not possible due to financial reasons",
          "患者側の希望または事情" = "Patient's preference or personal reasons",
          "提示された治療薬以外の化学療法を行った" = "Received chemotherapy other than the proposed regimen",
          "死亡" = "Patient deceased"
        )
      )
      # 列名マッピング辞書
      column_map <- c(
        "CTx_lines_before_CGP" = "CTx lines before CGP",
        "症例.背景情報.初回治療前のステージ分類.名称." = "Stage at diagnosis",
        "症例.基本情報.がん種.OncoTree." = "Diagnosis",
        "症例.基本情報.がん種.OncoTree..名称." = "Diagnosis (OncoTree)",
        "症例.基本情報.性別.名称." = "Sex",
        "症例.基本情報.年齢" = "Age at CGP",
        "症例.基本情報.年齢.診断時" = "Age at diagnosis",
        "症例.背景情報.家族歴有無.名称." = "Family cancer history",
        "症例.背景情報.喫煙歴有無.名称." = "Smoking history",
        "症例.背景情報.喫煙年数" = "Years of smoking",
        "症例.背景情報.アルコール多飲有無.名称." = "Alcoholic history",
        "症例.背景情報.ECOG.PS.名称." = "Performance status",
        "症例.背景情報.重複がん有無.異なる臓器..名称." = "Double cancer in a different organ",
        "症例.背景情報.多発がん有無.同一臓器..名称." = "Multiple tumor nodules in the organ",
        "症例.検体情報.パネル.名称." = "Cancer gene panel",
        "EP_option" = "Treatment recommended by expert panel",
        "EP_treat" = "Indicated treatment administered",
        "症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称." = "Reason for no indicated CTx",
        "症例.がん種情報.登録時転移有無.名称." = "Any distant metastasis",
        "Lymph_met" = "Lymphatic metastasis",
        "Brain_met" = "Brain metastasis",
        "Lung_met" = "Lung metastasis",
        "Bone_met" = "Bone metastasis",
        "Liver_met" = "Liver metastasis",
        "Other_met" = "Other organ metastasis",
        "Not_brain_bone_liver_met" = "Metastasis other than brain/bone/liver",
        "症例.がん種情報.肺.EGFR.名称." = "EGFR muts in lung",
        "症例.がん種情報.肺.ALK融合.名称." = "ALK fusion in lung",
        "症例.がん種情報.肺.ROS1.名称." = "ROS1 fusion in lung",
        "症例.がん種情報.肺.BRAF.V600..名称." = "BRAF V600E in lung",
        "症例.がん種情報.肺.MET遺伝子エクソン14スキッピング変異.名称." = "MET ex14 skipping in lung",
        "症例.がん種情報.肺.KRAS.G12C遺伝子変異.名称." = "KRAS G12C in lung",
        "症例.がん種情報.肺.RET融合遺伝子.名称." = "RET fusion in lung",
        "PD_L1" = "PD-L1 IHC in lung",
        "症例.がん種情報.乳.HER2.IHC..名称." = "HER2 IHC in breast",
        "症例.がん種情報.乳.HER2.FISH..名称." = "HER2 FISH in breast",
        "症例.がん種情報.乳.ER.名称." = "ER IHC in breast",
        "症例.がん種情報.乳.PgR.名称." = "PgR IHC in breast",
        "症例.がん種情報.固形がん.マイクロサテライト不安定性.名称." = "MSI PCR",
        "症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称." = "MMR IHC",
        "症例.がん種情報.食道.胃.腸.HER2.名称." = "HER2 IHC in stomach/bowel",
        "症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称." = "HER2 FISH in stomach/bowel"
      )
      # data.tableの場合
      setDT(Data_summary)
      setnames(Data_summary, names(column_map), column_map, skip_absent = TRUE)
      # YoungOld列の特別処理
      if ("YoungOld" %in% names(Data_summary)) {
        setnames(Data_summary, "YoungOld", paste0(input$mid_age + 1, " and older"))
      }
      # 対象列を特定
      target_cols <- intersect(names(translation_map), names(Data_summary))
      # 列変換
      for (col in target_cols) {
        Data_summary[is.na(get(col)), (col) := "Unknown"]
        Data_summary[, (col) := translation_map[[col]][get(col)]]
        Data_summary[is.na(get(col)), (col) := "Unknown"]
      }
      OUTPUT_DATA$table_summary_1_Data_summary = Data_summary
      incProgress(1 / 8)
    })
  })
  rm(analysis_env)
  gc()
}

output$table_summary_1 <- render_gt({
  req(OUTPUT_DATA$table_summary_1_Data_summary,
      input$table_summary,
      input$table_summary_layout)
  withProgress(message = "Drawing table...", {
    # データの前処理
    data <- OUTPUT_DATA$table_summary_1_Data_summary
    # 表示形式に応じた除外カラムとgrouping変数の設定
    data <- switch(input$table_summary,
                   "all" = data %>% select(-Diagnosis),
                   "gene" = data %>% select(-Diagnosis),
                   "panel" = data %>% select(-`Diagnosis (OncoTree)`),
                   "diagnosis" = data %>% select(-`Diagnosis (OncoTree)`)
    )

    group_var <- switch(input$table_summary,
                        "all" = NULL,
                        "gene" = "separation_value",
                        "panel" = "Cancer gene panel",
                        "diagnosis" = "Diagnosis"
    )

    caption <- switch(input$table_summary,
                      "all" = "**Patient Characteristics** (N = {N})",
                      "gene" = "**Patient Characteristics** (N = {N})",
                      "panel" = "**Patient Characteristics, by panel** (N = {N})",
                      "diagnosis" = "**Patient Characteristics, by diagnosis** (N = {N})"
    )
    # 集計表の作成
    if (input$table_summary_layout == "No") {
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
        tbl_output %>% as_gt()
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
      tbl_output %>% as_gt()
    }
  })
})
