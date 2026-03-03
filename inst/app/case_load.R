Data_case_raw_pre =  reactive({
  withProgress(message = "Clinical csv files loading.", {
    req(input$clinical_files)
    Encode = guess_encoding(input$clinical_files[[1, 'datapath']])[[1]][[1]]
    if(Encode == "Shift_JIS"){
      Encode = "CP932"
    }
    clin_tmp <- read_csv(col_names = TRUE,
                         file = input$clinical_files$datapath,
                         locale = locale(encoding=Encode),
                         num_threads=max(1, parallel::detectCores() - 1, na.rm = TRUE),
                         progress=FALSE,
                         show_col_types=FALSE,
                         name_repair=make.names,
                         col_select = c(
                           C.CAT調査結果.基本項目.ハッシュID,
                           症例.基本情報.年齢,
                           症例.背景情報.病理診断名,
                           症例.背景情報.臨床診断名,
                           症例.検体情報.病理診断名,
                           症例.基本情報.性別.名称.,
                           症例.背景情報.初回治療前のステージ分類.名称.,
                           症例.基本情報.がん種.OncoTree.,
                           症例.基本情報.がん種.OncoTree..名称.,
                           症例.基本情報.がん種.OncoTree.LEVEL1.,
                           症例.検体情報.パネル.名称.,
                           症例.背景情報.ECOG.PS.名称.,
                           症例.背景情報.家族歴有無.名称.,
                           症例.背景情報.喫煙歴有無.名称.,
                           症例.背景情報.喫煙年数,
                           #症例.背景情報.喫煙１日本数,
                           症例.背景情報.アルコール多飲有無.名称.,
                           症例.背景情報.重複がん有無.異なる臓器..名称.,
                           症例.背景情報.多発がん有無.同一臓器..名称.,
                           症例.がん種情報.登録時転移有無.名称.,
                           症例.がん種情報.登録時転移部位.名称.,
                           症例.検体情報.腫瘍細胞含有割合,
                           症例.検体情報.検体採取部位.名称.,
                           症例.検体情報.検体採取日.腫瘍組織.,
                           症例.EP前レジメン情報.治療ライン.名称.,
                           症例.EP前レジメン情報.実施目的.名称.,
                           症例.EP前レジメン情報.化学療法レジメン名称,
                           症例.EP前レジメン情報.薬剤名.YJ一般名.EN.,
                           症例.EP前レジメン情報.投与開始日,
                           症例.EP前レジメン情報.投与終了日,
                           症例.EP前レジメン情報.終了理由.名称.,
                           症例.EP前レジメン情報.レジメン継続区分.名称.,
                           症例.EP前レジメン情報.最良総合効果.名称.,
                           症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称.,
                           症例.EP後レジメン情報.提示された治療薬を投与した.名称.,
                           症例.EP後レジメン情報.治療ライン,
                           症例.EP後レジメン情報.治療ライン.名称.,
                           症例.EP後レジメン情報.化学療法レジメン名称,
                           症例.EP後レジメン情報.薬剤名.YJ一般名.EN.,
                           症例.EP後レジメン情報.投与開始日,
                           症例.EP後レジメン情報.投与終了日,
                           症例.EP後レジメン情報.エキスパートパネル開催日,
                           症例.EP後レジメン情報.終了理由.名称.,
                           症例.EP後レジメン情報.レジメン継続区分.名称.,
                           症例.EP後レジメン情報.最良総合効果.名称.,
                           症例.転帰情報.転帰.名称.,
                           症例.背景情報.診断日,
                           症例.管理情報.登録日,
                           症例.転帰情報.最終生存確認日,
                           症例.転帰情報.死亡日,
                           症例.転帰情報.死因.名称.,
                           症例.がん種情報.固形がん.マイクロサテライト不安定性.名称.,
                           症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称.,
                           症例.がん種情報.肺.EGFR.名称.,
                           症例.がん種情報.肺.ALK融合.名称.,
                           症例.がん種情報.肺.ROS1.名称.,
                           症例.がん種情報.肺.BRAF.V600..名称.,
                           症例.がん種情報.肺.MET遺伝子エクソン14スキッピング変異.名称.,
                           症例.がん種情報.肺.KRAS.G12C遺伝子変異.名称.,
                           症例.がん種情報.肺.RET融合遺伝子.名称.,
                           症例.がん種情報.肺.PD.L1.IHC..名称.,
                           症例.がん種情報.乳.HER2.IHC..名称.,
                           症例.がん種情報.乳.HER2.FISH..名称.,
                           症例.がん種情報.乳.ER.名称.,
                           症例.がん種情報.乳.PgR.名称.,
                           症例.がん種情報.食道.胃.腸.HER2.名称.,
                           症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称.,
                           症例.背景情報.発症年齢,
                           症例.がん種情報.肝.HBsAg.名称.,
                           症例.がん種情報.肝.HBs抗体.名称.,
                           症例.がん種情報.肝.HCV抗体.名称.,
                           症例.EP前レジメン情報.増悪確認日,
                           症例.EP前レジメン情報.Grade3以上有害事象の有無.名称.,
                           症例.EP前レジメン情報.中止に至った有害事象名.英語.,
                           症例.EP前レジメン情報.最悪Grade.名称.,
                           症例.EP前副作用情報.発症.覚知.日付,
                           症例.EP前副作用情報.CTCAEv5.0名称英語,
                           症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称.,
                           症例.EP前レジメン情報.最悪Grade.名称.,
                           症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称.,
                           症例.EP後レジメン情報.提示された治療薬とは異なる薬剤の投与有無.名称.,
                           症例.EP後レジメン情報.増悪確認日,
                           症例.EP後レジメン情報.Grade3以上有害事象の有無.名称.,
                           症例.EP後レジメン情報.中止に至った有害事象名.英語.,
                           症例.EP後レジメン情報.最悪Grade.名称.,
                           症例.EP後副作用情報.発症.覚知.日付,
                           症例.EP後副作用情報.CTCAEv5.0名称英語,
                           症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称.
                         ),
                         col_types = cols(
                           C.CAT調査結果.基本項目.ハッシュID = col_character(),
                           症例.基本情報.年齢 = col_integer(),
                           症例.背景情報.病理診断名 = col_character(),
                           症例.背景情報.臨床診断名 = col_character(),
                           症例.検体情報.病理診断名 = col_character(),
                           症例.基本情報.性別.名称. = col_character(),
                           症例.背景情報.初回治療前のステージ分類.名称. = col_character(),
                           症例.基本情報.がん種.OncoTree. = col_character(),
                           症例.基本情報.がん種.OncoTree..名称. = col_character(),
                           症例.基本情報.がん種.OncoTree.LEVEL1. = col_character(),
                           症例.検体情報.パネル.名称. = col_character(),
                           症例.背景情報.ECOG.PS.名称. = col_character(),
                           症例.背景情報.家族歴有無.名称. = col_character(),
                           症例.背景情報.喫煙歴有無.名称. = col_character(),
                           症例.背景情報.喫煙年数 = col_integer(),
                           #症例.背景情報.喫煙１日本数 = col_integer(),
                           症例.背景情報.アルコール多飲有無.名称. = col_character(),
                           症例.背景情報.重複がん有無.異なる臓器..名称. = col_character(),
                           症例.背景情報.多発がん有無.同一臓器..名称. = col_character(),
                           症例.がん種情報.登録時転移有無.名称. = col_character(),
                           症例.がん種情報.登録時転移部位.名称. = col_character(),
                           症例.検体情報.腫瘍細胞含有割合 = col_integer(),
                           症例.検体情報.検体採取部位.名称. = col_character(),
                           症例.検体情報.検体採取日.腫瘍組織. = col_character(),
                           症例.EP前レジメン情報.治療ライン.名称. = col_character(),
                           症例.EP前レジメン情報.実施目的.名称. = col_character(),
                           症例.EP前レジメン情報.化学療法レジメン名称 = col_character(),
                           症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = col_character(),
                           症例.EP前レジメン情報.投与開始日 = col_character(),
                           症例.EP前レジメン情報.投与終了日 = col_character(),
                           症例.EP前レジメン情報.終了理由.名称. = col_character(),
                           症例.EP前レジメン情報.レジメン継続区分.名称. = col_character(),
                           症例.EP前レジメン情報.最良総合効果.名称. = col_character(),
                           症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称. = col_character(),
                           症例.EP後レジメン情報.提示された治療薬を投与した.名称. = col_character(),
                           症例.EP後レジメン情報.治療ライン = col_integer(),
                           症例.EP後レジメン情報.治療ライン.名称. = col_character(),
                           症例.EP後レジメン情報.化学療法レジメン名称 = col_character(),
                           症例.EP後レジメン情報.薬剤名.YJ一般名.EN. = col_character(),
                           症例.EP後レジメン情報.投与開始日 = col_character(),
                           症例.EP後レジメン情報.投与終了日 = col_character(),
                           症例.EP後レジメン情報.エキスパートパネル開催日 = col_character(),
                           症例.EP後レジメン情報.終了理由.名称. = col_character(),
                           症例.EP後レジメン情報.レジメン継続区分.名称. = col_character(),
                           症例.EP後レジメン情報.最良総合効果.名称. = col_character(),
                           症例.転帰情報.転帰.名称. = col_character(),
                           症例.背景情報.診断日 = col_character(),
                           症例.管理情報.登録日 = col_character(),
                           症例.転帰情報.最終生存確認日 = col_character(),
                           症例.転帰情報.死亡日 = col_character(),
                           症例.転帰情報.死因.名称. = col_character(),
                           症例.がん種情報.固形がん.マイクロサテライト不安定性.名称. = col_character(),
                           症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称. = col_character(),
                           症例.がん種情報.肺.EGFR.名称. = col_character(),
                           症例.がん種情報.肺.ALK融合.名称. = col_character(),
                           症例.がん種情報.肺.ROS1.名称. = col_character(),
                           症例.がん種情報.肺.BRAF.V600..名称. = col_character(),
                           症例.がん種情報.肺.MET遺伝子エクソン14スキッピング変異.名称. = col_character(),
                           症例.がん種情報.肺.KRAS.G12C遺伝子変異.名称. = col_character(),
                           症例.がん種情報.肺.RET融合遺伝子.名称. = col_character(),
                           症例.がん種情報.肺.PD.L1.IHC..名称. = col_character(),
                           症例.がん種情報.乳.HER2.IHC..名称. = col_character(),
                           症例.がん種情報.乳.HER2.FISH..名称. = col_character(),
                           症例.がん種情報.乳.ER.名称. = col_character(),
                           症例.がん種情報.乳.PgR.名称. = col_character(),
                           症例.がん種情報.食道.胃.腸.HER2.名称. = col_character(),
                           症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称. = col_character(),
                           症例.背景情報.発症年齢 = col_integer(),
                           症例.がん種情報.肝.HBsAg.名称. = col_character(),
                           症例.がん種情報.肝.HBs抗体.名称. = col_character(),
                           症例.がん種情報.肝.HCV抗体.名称. = col_character(),
                           症例.EP前レジメン情報.増悪確認日 = col_character(),
                           症例.EP前レジメン情報.Grade3以上有害事象の有無.名称. = col_character(),
                           症例.EP前レジメン情報.中止に至った有害事象名.英語. = col_character(),
                           症例.EP前レジメン情報.最悪Grade.名称. = col_character(),
                           症例.EP前副作用情報.発症.覚知.日付 = col_character(),
                           症例.EP前副作用情報.CTCAEv5.0名称英語 = col_character(),
                           症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称. = col_character(),
                           症例.EP前レジメン情報.最悪Grade.名称. = col_character(),
                           症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称. = col_character(),
                           症例.EP後レジメン情報.提示された治療薬とは異なる薬剤の投与有無.名称. = col_character(),
                           症例.EP後レジメン情報.増悪確認日 = col_character(),
                           症例.EP後レジメン情報.Grade3以上有害事象の有無.名称. = col_character(),
                           症例.EP後レジメン情報.中止に至った有害事象名.英語. = col_character(),
                           症例.EP後レジメン情報.最悪Grade.名称. = col_character(),
                           症例.EP後副作用情報.発症.覚知.日付 = col_character(),
                           症例.EP後副作用情報.CTCAEv5.0名称英語 = col_character(),
                           症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称. = col_character(),
                           .default = col_guess()
                         )) %>%
      dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
    incProgress(1 / 4)
    clin_tmp = data.table(clin_tmp)
    Lymph_met_ID = unique((clin_tmp %>% dplyr::filter(
      症例.がん種情報.登録時転移部位.名称. %in% c("リンパ節", "リンパ節/リンパ管")))$C.CAT調査結果.基本項目.ハッシュID)
    Brain_met_ID = unique((clin_tmp %>% dplyr::filter(
      症例.がん種情報.登録時転移部位.名称. %in% c("脳", "中枢神経系")))$C.CAT調査結果.基本項目.ハッシュID)
    Lung_met_ID = unique((clin_tmp %>% dplyr::filter(
      症例.がん種情報.登録時転移部位.名称. %in% c("肺")))$C.CAT調査結果.基本項目.ハッシュID)
    Bone_met_ID = unique((clin_tmp %>% dplyr::filter(
      症例.がん種情報.登録時転移部位.名称. %in% c("骨髄", "骨")))$C.CAT調査結果.基本項目.ハッシュID)
    Liver_met_ID = unique((clin_tmp %>% dplyr::filter(
      症例.がん種情報.登録時転移部位.名称. %in% c("肝")))$C.CAT調査結果.基本項目.ハッシュID)
    Other_met_ID = unique((clin_tmp %>% dplyr::filter(
      !症例.がん種情報.登録時転移部位.名称. %in%
        c("リンパ節", "リンパ節/リンパ管", "脳", "中枢神経系", "肺", "骨髄", "骨", "肝", "")))$C.CAT調査結果.基本項目.ハッシュID)
    Not_brain_bone_liver_met_ID = unique((clin_tmp %>% dplyr::filter(
      !症例.がん種情報.登録時転移部位.名称. %in%
        c("脳", "中枢神経系", "骨髄", "骨", "肝", "")))$C.CAT調査結果.基本項目.ハッシュID)
    EP_option_ID = unique((clin_tmp %>% dplyr::filter(
      症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称. == "はい"))$C.CAT調査結果.基本項目.ハッシュID)
    EP_treat_ID = unique((clin_tmp %>% dplyr::filter(
      症例.EP後レジメン情報.提示された治療薬を投与した.名称. == "投与した"))$C.CAT調査結果.基本項目.ハッシュID)
    EP_treat_data = clin_tmp %>%
      dplyr::filter(症例.EP後レジメン情報.提示された治療薬を投与した.名称. == "投与した") %>%
      dplyr::select(C.CAT調査結果.基本項目.ハッシュID, 症例.EP後レジメン情報.投与開始日) %>%
      dplyr::arrange(症例.EP後レジメン情報.投与開始日) %>%
      dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID,.keep_all = T)
    colnames(EP_treat_data) = c("C.CAT調査結果.基本項目.ハッシュID", "EP_treat_date")
    # Define column names
    id_col     <- "C.CAT調査結果.基本項目.ハッシュID"
    reason_col <- "症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称."
    # Broadcast the first non-NA value to all rows within each group
    clin_tmp[, (reason_col) := {
      # Get the vector for the current group
      x <- .SD[[reason_col]]
      # Pick the first non-NA value
      v <- x[!is.na(x)][1L]
      # If v is NA (all values in group are NA), return x.
      # Otherwise, return v (it will be recycled to all rows in the group).
      if (is.na(v)) x else v
    }, by = id_col, .SDcols = reason_col]
    clin_tmp = clin_tmp %>%
      left_join(EP_treat_data, by="C.CAT調査結果.基本項目.ハッシュID")
    clin_tmp = clin_tmp %>%
      dplyr::mutate(
        症例.背景情報.ECOG.PS.名称. = as.character(症例.背景情報.ECOG.PS.名称.),
        症例.基本情報.がん種.OncoTree.LEVEL1. = case_when(
          症例.基本情報.がん種.OncoTree.LEVEL1. == "WHO_BRAIN" ~ "BRAIN",
          TRUE ~ 症例.基本情報.がん種.OncoTree.LEVEL1.
        ),
        症例.基本情報.がん種.OncoTree..名称. =
          str_split(str_replace(症例.基本情報.がん種.OncoTree..名称.,
                                "\\)", "\\)__"), "__", simplify = TRUE)[,1],
        症例.基本情報.がん種.OncoTree..名称. = case_when(
          str_detect(症例.基本情報.がん種.OncoTree..名称., paste0(" \\(", 症例.基本情報.がん種.OncoTree.LEVEL1., "\\)")) ~ 症例.基本情報.がん種.OncoTree.LEVEL1.,
          TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
        ),
        症例.基本情報.がん種.OncoTree. = str_replace(str_replace(症例.基本情報.がん種.OncoTree..名称., "\\)", ""), ".*\\(",""),
        症例.基本情報.がん種.OncoTree. = case_when(
          is.na(症例.基本情報.がん種.OncoTree.) ~ 症例.基本情報.がん種.OncoTree.LEVEL1.,
          TRUE ~ 症例.基本情報.がん種.OncoTree.
        ),症例.がん種情報.登録時転移有無.名称. = case_when(
          is.na(症例.がん種情報.登録時転移有無.名称.) ~ "不明",
          症例.がん種情報.登録時転移有無.名称. == "" ~ "不明",
          TRUE ~ 症例.がん種情報.登録時転移有無.名称.
        ),
        症例.背景情報.家族歴有無.名称. = case_when(
          is.na(症例.背景情報.家族歴有無.名称.) ~ "不明",
          症例.背景情報.家族歴有無.名称. == "" ~ "不明",
          TRUE ~ 症例.背景情報.家族歴有無.名称.
        ),
        症例.背景情報.ECOG.PS.名称. = case_when(
          is.na(症例.背景情報.ECOG.PS.名称.) ~ "不明",
          症例.背景情報.ECOG.PS.名称. == "" ~ "不明",
          TRUE ~ 症例.背景情報.ECOG.PS.名称.
        ),
        症例.基本情報.性別.名称. = case_when(
          is.na(症例.基本情報.性別.名称.) ~ "未入力・不明",
          症例.基本情報.性別.名称. == "" ~ "未入力・不明",
          TRUE ~ 症例.基本情報.性別.名称.
        ),
        症例.EP前副作用情報.CTCAEv5.0名称英語 = case_when(
          症例.EP前副作用情報.CTCAEv5.0名称英語 == "12085330" ~ "Neutrophil count decreased",
          症例.EP前副作用情報.CTCAEv5.0名称英語 == "Neutrophil　count decreased" ~ "Neutrophil count decreased",
          TRUE ~ 症例.EP前副作用情報.CTCAEv5.0名称英語
        ),
        症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称. = case_when(
          症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称. == "不明" ~ "Grade  Unknown",
          TRUE ~ 症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称.
        ),
        症例.EP後副作用情報.CTCAEv5.0名称英語 = case_when(
          症例.EP後副作用情報.CTCAEv5.0名称英語 == "12085330" ~ "Neutrophil count decreased",
          症例.EP後副作用情報.CTCAEv5.0名称英語 == "Neutrophil　count decreased" ~ "Neutrophil count decreased",
          TRUE ~ 症例.EP後副作用情報.CTCAEv5.0名称英語
        ),
        症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称. = case_when(
          症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称. == "不明" ~ "Grade  Unknown",
          TRUE ~ 症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称.
        ),
        across(
          c(
            EP_treat_date,
            症例.EP前レジメン情報.投与開始日,
            症例.EP後レジメン情報.投与開始日,
            症例.EP前レジメン情報.投与終了日,
            症例.EP後レジメン情報.投与終了日,
            症例.検体情報.検体採取日.腫瘍組織.,
            症例.EP後レジメン情報.エキスパートパネル開催日,
            症例.EP前レジメン情報.増悪確認日,
            症例.EP前副作用情報.発症.覚知.日付,
            症例.EP後レジメン情報.増悪確認日,
            症例.EP後副作用情報.発症.覚知.日付
          ),
          convert_date
        ),
        Lymph_met = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% Lymph_met_ID,
          "Yes", "No"
        ),
        Brain_met = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% Brain_met_ID,
          "Yes", "No"
        ),
        Lung_met = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% Lung_met_ID,
          "Yes", "No"
        ),
        Bone_met = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% Bone_met_ID,
          "Yes", "No"
        ),
        Liver_met = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% Liver_met_ID,
          "Yes", "No"
        ),
        Other_met = if_else(
          !C.CAT調査結果.基本項目.ハッシュID %in% Other_met_ID,
          "Yes", "No"
        ),
        Not_brain_bone_liver_met = if_else(
          !C.CAT調査結果.基本項目.ハッシュID %in% Not_brain_bone_liver_met_ID,
          "Yes", "No"
        ),
        HER2_IHC = case_when(
          症例.がん種情報.乳.HER2.IHC..名称. == "陽性（3+）" |
            症例.がん種情報.食道.胃.腸.HER2.名称. == "陽性（3+）" |
            症例.がん種情報.乳.HER2.IHC..名称. == "境界域（2+）" & 症例.がん種情報.乳.HER2.FISH..名称. == "陽性" |
            症例.がん種情報.食道.胃.腸.HER2.名称. == "境界域（2+）" & 症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称. == "陽性" ~ "Positive",
          症例.がん種情報.乳.HER2.IHC..名称. %in% c("陰性",  "陰性（1+）") |
            症例.がん種情報.食道.胃.腸.HER2.名称. %in% c("陰性",  "陰性（1+）") |
            症例.がん種情報.乳.HER2.IHC..名称. == "境界域（2+）" & 症例.がん種情報.乳.HER2.FISH..名称. %in% c("equivocal", "陰性") |
            症例.がん種情報.食道.胃.腸.HER2.名称. == "境界域（2+）" & 症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称. %in% c("equivocal", "陰性") ~ "Negative",
          TRUE ~ "Unknown"
        ),
        MSI_PCR = case_when(
          症例.がん種情報.固形がん.マイクロサテライト不安定性.名称. == "陽性" ~ "Positive",
          症例.がん種情報.固形がん.マイクロサテライト不安定性.名称. == "陰性" ~ "Negative",
          TRUE ~ "Unknown"
        ),
        MMR_IHC = case_when(
          症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称. == "dMMR(欠損)" ~ "dMMR",
          症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称. == "pMMR(正常)" ~ "pMMR",
          TRUE ~ "Unknown"
        ),
        PD_L1 = case_when(
          症例.がん種情報.肺.PD.L1.IHC..名称. == "陰性" ~ "Negative",
          症例.がん種情報.肺.PD.L1.IHC..名称. == "陽性" ~ "Positive",
          症例.がん種情報.肺.PD.L1.IHC..名称. == "不明or未検査" ~ "Unknown",
          症例.がん種情報.肺.PD.L1.IHC..名称. == "判定不能" ~ "Undetermined",
          TRUE ~ "Unknown"
        ),
        EP_option = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% EP_option_ID,
          1, 0
        ),
        EP_treat = if_else(
          C.CAT調査結果.基本項目.ハッシュID %in% EP_treat_ID,
          1, 0
        ),
        症例.背景情報.アルコール多飲有無.名称. = case_when(
          is.na(症例.背景情報.アルコール多飲有無.名称.) ~ "Unknown",
          症例.背景情報.アルコール多飲有無.名称. == "あり" ~ "Yes",
          症例.背景情報.アルコール多飲有無.名称. == "なし" ~ "No",
          症例.背景情報.アルコール多飲有無.名称. == "不明" ~ "Unknown",
          症例.背景情報.アルコール多飲有無.名称. == "" ~ "Unknown",
          TRUE ~ "Unknown"
        ),
        症例.背景情報.重複がん有無.異なる臓器..名称. = case_when(
          症例.背景情報.重複がん有無.異なる臓器..名称. == "" ~ "不明",
          TRUE ~ 症例.背景情報.重複がん有無.異なる臓器..名称.
        ),
        症例.背景情報.多発がん有無.同一臓器..名称. = case_when(
          症例.背景情報.多発がん有無.同一臓器..名称. == "" ~ "不明",
          TRUE ~ 症例.背景情報.多発がん有無.同一臓器..名称.
        ),
        症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称. = case_when(
          is.na(症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称.) ~ "いいえ",
          症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称. == "" ~ "いいえ",
          TRUE ~ 症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称.
        ),
        症例.EP後レジメン情報.提示された治療薬を投与した.名称. = case_when(
          is.na(症例.EP後レジメン情報.提示された治療薬を投与した.名称.) ~ "投与しなかった",
          症例.EP後レジメン情報.提示された治療薬を投与した.名称. == "不明" ~ "投与しなかった",
          症例.EP後レジメン情報.提示された治療薬を投与した.名称. == "" ~ "投与しなかった",
          TRUE ~ 症例.EP後レジメン情報.提示された治療薬を投与した.名称.
        ),
        症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称. = case_when(
          is.na(症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称.) ~ "不明",
          症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称. == "その他・不明" ~ "不明",
          症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称. == "その他" ~ "不明",
          症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称. == "" ~ "不明",
          TRUE ~ 症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称.
        ),
        症例.背景情報.喫煙歴有無.名称. = case_when(
          is.na(症例.背景情報.喫煙歴有無.名称.) ~ "不明",
          症例.背景情報.喫煙歴有無.名称. == "" ~ "不明",
          TRUE ~ 症例.背景情報.喫煙歴有無.名称.
        ),
        across(
          c(
            症例.検体情報.検体採取部位.名称.,
            症例.背景情報.ECOG.PS.名称.,
            症例.背景情報.重複がん有無.異なる臓器..名称.,
            症例.背景情報.多発がん有無.同一臓器..名称.
          ),
          ~ replace_na(.x, "不明")
        ),
        across(
          c(
            症例.背景情報.喫煙歴有無.名称.,
            症例.検体情報.パネル.名称.,
            症例.がん種情報.登録時転移部位.名称.,
            症例.転帰情報.死亡日,
            症例.転帰情報.最終生存確認日,
            症例.管理情報.登録日,
            EP_treat_date,
            症例.EP前レジメン情報.投与開始日,
            症例.EP前レジメン情報.投与終了日,
            症例.EP後レジメン情報.投与開始日,
            症例.EP後レジメン情報.投与終了日,
            症例.EP後レジメン情報.エキスパートパネル開催日,
            症例.EP前レジメン情報.増悪確認日,
            症例.EP後レジメン情報.増悪確認日,
            症例.EP前副作用情報.発症.覚知.日付,
            症例.EP後副作用情報.発症.覚知.日付,
            症例.検体情報.検体採取日.腫瘍組織.,
            症例.背景情報.診断日
          ),
          ~ replace_na(.x, "")
        ),
        across(
          c(
            症例.管理情報.登録日,
            症例.転帰情報.最終生存確認日,
            症例.転帰情報.死亡日,
            EP_treat_date,
            症例.EP前レジメン情報.投与開始日,
            症例.EP前レジメン情報.投与終了日,
            症例.EP後レジメン情報.投与開始日,
            症例.EP後レジメン情報.投与終了日,
            症例.EP前レジメン情報.増悪確認日,
            症例.EP後レジメン情報.増悪確認日,
            症例.EP後レジメン情報.エキスパートパネル開催日,
            症例.EP前副作用情報.発症.覚知.日付,
            症例.EP後副作用情報.発症.覚知.日付,
            症例.検体情報.検体採取日.腫瘍組織.,
            症例.背景情報.診断日
          ),
          ~ str_replace_all(.x, "/", "-")
        ),
        across(
          c(
            症例.転帰情報.死亡日,
            症例.転帰情報.最終生存確認日,
            症例.管理情報.登録日,
            EP_treat_date,
            症例.EP前レジメン情報.投与開始日,
            症例.EP前レジメン情報.投与終了日,
            症例.EP後レジメン情報.投与開始日,
            症例.EP後レジメン情報.投与終了日,
            症例.EP後レジメン情報.エキスパートパネル開催日,
            症例.EP前レジメン情報.増悪確認日,
            症例.EP後レジメン情報.増悪確認日,
            症例.EP前副作用情報.発症.覚知.日付,
            症例.EP後副作用情報.発症.覚知.日付,
            症例.検体情報.検体採取日.腫瘍組織.,
            症例.背景情報.診断日
          ),
          ~ case_when(
            . == as.Date("9999-12-31") ~ "",  # 日付型で比較
            as.character(.) == "9999-12-31" ~ "",  # 文字列でも念のため対応
            TRUE ~ as.character(.)         # それ以外は文字列として残す
          )
        ),
        症例.基本情報.がん種.OncoTree..名称. = case_when(
          is.na(症例.基本情報.がん種.OncoTree..名称.) ~ 症例.基本情報.がん種.OncoTree.LEVEL1.,
          症例.基本情報.がん種.OncoTree..名称. == "" ~ 症例.基本情報.がん種.OncoTree.LEVEL1.,
          TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
        ),
        症例.基本情報.性別.名称. = case_when(
          症例.基本情報.性別.名称. == "女" ~ "Female",
          症例.基本情報.性別.名称. == "未入力・不明" ~ "Unknown",
          症例.基本情報.性別.名称. == "男" ~ "Male"
        ),
        症例.背景情報.喫煙歴有無.名称. = case_when(
          症例.背景情報.喫煙歴有無.名称. == "あり" ~ "Yes",
          症例.背景情報.喫煙歴有無.名称. == "不明" ~ "Unknown",
          症例.背景情報.喫煙歴有無.名称. == "なし" ~ "No"
        ),
        症例.背景情報.ECOG.PS.名称. = case_when(
          症例.背景情報.ECOG.PS.名称. == "不明" ~ "Unknown",
          TRUE ~ 症例.背景情報.ECOG.PS.名称.
        ),
        症例.背景情報.喫煙年数 = case_when(
          症例.背景情報.喫煙年数 == 99 ~ NA_integer_,
          症例.背景情報.喫煙歴有無.名称. == "No" ~ NA_integer_,
          TRUE ~ 症例.背景情報.喫煙年数
        ),
        症例.背景情報.初回治療前のステージ分類.名称. = case_when(
          is.na(症例.背景情報.初回治療前のステージ分類.名称.) ~ "Unknown",
          症例.背景情報.初回治療前のステージ分類.名称. == "0期" ~ "0",
          症例.背景情報.初回治療前のステージ分類.名称. == "Ⅰ期" ~ "1",
          症例.背景情報.初回治療前のステージ分類.名称. == "Ⅱ期" ~ "2",
          症例.背景情報.初回治療前のステージ分類.名称. == "Ⅲ期" ~ "3",
          症例.背景情報.初回治療前のステージ分類.名称. == "Ⅳ期" ~ "4",
          TRUE ~ "Unknown"
        )
      )
    if(!is.null(input$ID_histology)){
      ID_histology_list = data.frame(NULL)
      for(i in 1:length(input$ID_histology[,1])){
        ID_histology_list <- rbind(ID_histology_list,
                                   read.csv(header = TRUE,
                                            file(input$ID_histology[[i, 'datapath']],
                                                 encoding='UTF-8-BOM')))
      }
      clin_tmp$症例.基本情報.がん種.OncoTree..名称._tmp = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
        as.vector(ID_histology_list$Histology[match(x, ID_histology_list$ID)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.基本情報.がん種.OncoTree..名称. = case_when(
          !is.na(症例.基本情報.がん種.OncoTree..名称._tmp) ~ 症例.基本情報.がん種.OncoTree..名称._tmp,
          TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
        )
      )
      clin_tmp$症例.基本情報.がん種.OncoTree._tmp = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
        as.vector(ID_histology_list$Histology[match(x, ID_histology_list$ID)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.基本情報.がん種.OncoTree. = case_when(
          !is.na(症例.基本情報.がん種.OncoTree._tmp) ~ 症例.基本情報.がん種.OncoTree._tmp,
          TRUE ~ 症例.基本情報.がん種.OncoTree.
        )
      )
    }
    if(!is.null(input$regimen_rename)){
      RenList = data.frame(NULL)
      for(i in 1:length(input$regimen_rename[,1])){
        RenList <- rbind(RenList,
                         read.csv(header = TRUE,
                                  file(input$regimen_rename[[i, 'datapath']],
                                       encoding='UTF-8-BOM')))
      }
      for(i in 1:length(RenList$P1)){
        clin_tmp = clin_tmp %>% dplyr::mutate(
          症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = case_when(
            症例.EP前レジメン情報.薬剤名.YJ一般名.EN. == "" &
              症例.EP前レジメン情報.化学療法レジメン名称 == RenList$P1[i] ~ RenList$P2[i],
            TRUE ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN.
          ))
      }
    }
    if(!is.null(input$drug_rename)){
      drug_rename_list = data.frame(NULL)
      for(i in 1:length(input$drug_rename[,1])){
        drug_rename_list <- rbind(drug_rename_list,
                                  read.csv(header = TRUE,
                                           file(input$drug_rename[[i, 'datapath']],
                                                encoding='UTF-8-BOM')))
      }
      clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp = unlist(lapply(list(clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN.), function(x) {
        as.vector(drug_rename_list$Rename[match(x, drug_rename_list$Drug)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = case_when(
          !is.na(症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp) ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp,
          TRUE ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN.
        )
      )
    }
    if(!is.null(input$histology_rename)){
      histology_rename_list = data.frame(NULL)
      for(i in 1:length(input$histology_rename[,1])){
        histology_rename_list <- rbind(histology_rename_list,
                                       read.csv(header = TRUE,
                                                file(input$histology_rename[[i, 'datapath']],
                                                     encoding='UTF-8-BOM')))
      }
      clin_tmp$症例.基本情報.がん種.OncoTree..名称._tmp = unlist(lapply(list(clin_tmp$症例.基本情報.がん種.OncoTree..名称.), function(x) {
        as.vector(histology_rename_list$Rename[match(x, histology_rename_list$Histology)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.基本情報.がん種.OncoTree..名称. = case_when(
          !is.na(症例.基本情報.がん種.OncoTree..名称._tmp) ~ 症例.基本情報.がん種.OncoTree..名称._tmp,
          TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
        )
      )
      clin_tmp$症例.基本情報.がん種.OncoTree._tmp = unlist(lapply(list(clin_tmp$症例.基本情報.がん種.OncoTree.), function(x) {
        as.vector(histology_rename_list$Rename[match(x, histology_rename_list$Histology)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.基本情報.がん種.OncoTree. = case_when(
          !is.na(症例.基本情報.がん種.OncoTree._tmp) ~ 症例.基本情報.がん種.OncoTree._tmp,
          TRUE ~ 症例.基本情報.がん種.OncoTree.
        )
      )
    }
    RenList = read.csv(file(app_path("source/drug_rename.csv"),
                            encoding='UTF-8-BOM'), header = T)
    # for(i in 1:length(RenList$P1)){
    #   clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = str_replace_all(clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN., RenList$P1[i], RenList$P2[i])
    #   clin_tmp$症例.EP後レジメン情報.薬剤名.YJ一般名.EN. = str_replace_all(clin_tmp$症例.EP後レジメン情報.薬剤名.YJ一般名.EN., RenList$P1[i], RenList$P2[i])
    # }
    replace_map <- setNames(RenList$P2, RenList$P1)
    clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN. <- str_replace_all(clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN., replace_map)
    clin_tmp$症例.EP後レジメン情報.薬剤名.YJ一般名.EN. <- str_replace_all(clin_tmp$症例.EP後レジメン情報.薬剤名.YJ一般名.EN., replace_map)

    if(!is.null(input$regimen_rename)){
      RenList = data.frame(NULL)
      for(i in 1:length(input$regimen_rename[,1])){
        RenList <- rbind(RenList,
                         read.csv(header = TRUE,
                                  file(input$regimen_rename[[i, 'datapath']],
                                       encoding='UTF-8-BOM')))
      }
      for(i in 1:length(RenList$P1)){
        clin_tmp = clin_tmp %>% dplyr::mutate(
          症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = case_when(
            症例.EP前レジメン情報.薬剤名.YJ一般名.EN. == "" &
              症例.EP前レジメン情報.化学療法レジメン名称 == RenList$P1[i] ~ RenList$P2[i],
            TRUE ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN.
          ))
      }
    }
    if(!is.null(input$drug_rename)){
      drug_rename_list = data.frame(NULL)
      for(i in 1:length(input$drug_rename[,1])){
        drug_rename_list <- rbind(drug_rename_list,
                                  read.csv(header = TRUE,
                                           file(input$drug_rename[[i, 'datapath']],
                                                encoding='UTF-8-BOM')))
      }
      clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp = unlist(lapply(list(clin_tmp$症例.EP前レジメン情報.薬剤名.YJ一般名.EN.), function(x) {
        as.vector(drug_rename_list$Rename[match(x, drug_rename_list$Drug)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        症例.EP前レジメン情報.薬剤名.YJ一般名.EN. = case_when(
          !is.na(症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp) ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN._tmp,
          TRUE ~ 症例.EP前レジメン情報.薬剤名.YJ一般名.EN.
        )
      )
    }
    incProgress(1 / 4)
    Data_survival = clin_tmp %>%
      dplyr::select(C.CAT調査結果.基本項目.ハッシュID, 症例.管理情報.登録日) %>%
      dplyr::arrange(C.CAT調査結果.基本項目.ハッシュID, desc(症例.管理情報.登録日)) %>%
      dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = TRUE)
    clin_tmp = clin_tmp %>%
      dplyr::select(-症例.管理情報.登録日)
    clin_tmp = left_join(clin_tmp, Data_survival, by="C.CAT調査結果.基本項目.ハッシュID")

    Data_survival = clin_tmp %>%
      dplyr::mutate(final_observe = case_when(
        症例.転帰情報.最終生存確認日 != "" & 症例.転帰情報.最終生存確認日 != "           " & !is.na(症例.転帰情報.最終生存確認日) ~ 症例.転帰情報.最終生存確認日,
        症例.転帰情報.死亡日 != "" & !is.na(症例.転帰情報.死亡日) & 症例.転帰情報.死亡日 != "           " ~ 症例.転帰情報.死亡日,
        症例.EP後レジメン情報.投与終了日 != "" & !is.na(症例.EP後レジメン情報.投与終了日) & 症例.EP後レジメン情報.投与終了日 != "           " ~ 症例.EP後レジメン情報.投与終了日,
        症例.EP後レジメン情報.投与開始日 != "" & !is.na(症例.EP後レジメン情報.投与開始日) & 症例.EP後レジメン情報.投与開始日 != "           " ~ 症例.EP後レジメン情報.投与開始日,
        症例.管理情報.登録日 != "" & !is.na(症例.管理情報.登録日) & 症例.管理情報.登録日 != "           " ~ 症例.管理情報.登録日,
        TRUE ~ NA_character_))
    Data_survival = Data_survival %>%
      dplyr::mutate(censor = case_when(
        症例.転帰情報.死亡日 != ""  & !is.na(症例.転帰情報.死亡日) & 症例.転帰情報.死亡日 != "           " ~ 1,
        TRUE ~ 0)) %>%
      dplyr::arrange(C.CAT調査結果.基本項目.ハッシュID)
    Data_survival_tmp = Data_survival %>%
      dplyr::select(C.CAT調査結果.基本項目.ハッシュID, censor) %>%
      dplyr::arrange(desc(censor)) %>%
      dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = TRUE)
    clin_tmp = left_join(clin_tmp, Data_survival_tmp, by="C.CAT調査結果.基本項目.ハッシュID")
    Data_survival = Data_survival %>%
      dplyr::select(C.CAT調査結果.基本項目.ハッシュID, final_observe) %>%
      dplyr::arrange(desc(final_observe)) %>%
      dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = TRUE)
    date_columns <- c("final_observe", "EP_treat_date", "症例.管理情報.登録日",
                      "症例.EP前レジメン情報.投与開始日", "症例.EP前レジメン情報.投与終了日",
                      "症例.EP後レジメン情報.投与開始日", "症例.EP後レジメン情報.投与終了日",
                      "症例.EP前レジメン情報.増悪確認日", "症例.EP後レジメン情報.増悪確認日",
                      "症例.EP前副作用情報.発症.覚知.日付", "症例.EP後副作用情報.発症.覚知.日付",
                      "症例.背景情報.診断日")
    clin_tmp = left_join(clin_tmp, Data_survival, by="C.CAT調査結果.基本項目.ハッシュID") %>%
      dplyr::mutate(across(all_of(date_columns), as.Date))
    clin_tmp$症例.EP前レジメン情報.化学療法レジメン名称[is.na(clin_tmp$症例.EP前レジメン情報.化学療法レジメン名称)] = ""
    clin_tmp$症例.EP後レジメン情報.化学療法レジメン名称[is.na(clin_tmp$症例.EP後レジメン情報.化学療法レジメン名称)] = ""
    clin_tmp$症例.EP前レジメン情報.最良総合効果.名称.[is.na(clin_tmp$症例.EP前レジメン情報.最良総合効果.名称.)] = "Unknown"
    clin_tmp$症例.EP後レジメン情報.最良総合効果.名称.[is.na(clin_tmp$症例.EP後レジメン情報.最良総合効果.名称.)] = "Unknown"
    clin_tmp = clin_tmp %>%
      dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                    症例.基本情報.年齢,
                    Lymph_met,
                    Brain_met,
                    Lung_met,
                    Bone_met,
                    Liver_met,
                    Other_met,
                    Not_brain_bone_liver_met,
                    EP_option,
                    EP_treat,
                    final_observe,
                    censor,
                    症例.背景情報.病理診断名,
                    症例.背景情報.臨床診断名,
                    症例.検体情報.病理診断名,
                    症例.基本情報.性別.名称.,
                    症例.背景情報.初回治療前のステージ分類.名称.,
                    症例.基本情報.がん種.OncoTree.,
                    症例.基本情報.がん種.OncoTree..名称.,
                    症例.基本情報.がん種.OncoTree.LEVEL1.,
                    症例.検体情報.パネル.名称.,
                    症例.背景情報.ECOG.PS.名称.,
                    症例.背景情報.家族歴有無.名称.,
                    症例.背景情報.喫煙歴有無.名称.,
                    症例.背景情報.喫煙年数,
                    #症例.背景情報.喫煙１日本数,
                    症例.背景情報.アルコール多飲有無.名称.,
                    症例.背景情報.重複がん有無.異なる臓器..名称.,
                    症例.背景情報.多発がん有無.同一臓器..名称.,
                    症例.がん種情報.登録時転移有無.名称.,
                    症例.がん種情報.登録時転移部位.名称.,
                    症例.検体情報.腫瘍細胞含有割合,
                    症例.検体情報.検体採取部位.名称.,
                    症例.検体情報.検体採取日.腫瘍組織.,
                    症例.EP前レジメン情報.治療ライン.名称.,
                    症例.EP前レジメン情報.実施目的.名称.,
                    症例.EP前レジメン情報.化学療法レジメン名称,
                    症例.EP前レジメン情報.薬剤名.YJ一般名.EN.,
                    症例.EP前レジメン情報.投与開始日,
                    症例.EP前レジメン情報.投与終了日,
                    症例.EP前レジメン情報.終了理由.名称.,
                    症例.EP前レジメン情報.レジメン継続区分.名称.,
                    症例.EP前レジメン情報.最良総合効果.名称.,
                    症例.EP後レジメン情報.EPの結果新たな治療薬の選択肢が提示された.名称.,
                    症例.EP後レジメン情報.提示された治療薬を投与した.名称.,
                    症例.EP後レジメン情報.治療ライン,
                    症例.EP後レジメン情報.治療ライン.名称.,
                    症例.EP後レジメン情報.化学療法レジメン名称,
                    症例.EP後レジメン情報.薬剤名.YJ一般名.EN.,
                    症例.EP後レジメン情報.投与開始日,
                    症例.EP後レジメン情報.投与終了日,
                    症例.EP後レジメン情報.終了理由.名称.,
                    症例.EP後レジメン情報.レジメン継続区分.名称.,
                    症例.EP後レジメン情報.最良総合効果.名称.,
                    症例.転帰情報.転帰.名称.,
                    症例.背景情報.診断日,
                    症例.管理情報.登録日,
                    症例.転帰情報.最終生存確認日,
                    症例.転帰情報.死亡日,
                    症例.がん種情報.固形がん.マイクロサテライト不安定性.名称.,
                    症例.がん種情報.固形がん.ミスマッチ修復機能欠損.名称.,
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
                    症例.がん種情報.食道.胃.腸.HER2.名称.,
                    症例.がん種情報.食道.胃.腸.HER2遺伝子増幅..ISH法..名称.,
                    HER2_IHC,
                    MSI_PCR,
                    MMR_IHC,
                    症例.EP前レジメン情報.増悪確認日,
                    症例.EP前レジメン情報.Grade3以上有害事象の有無.名称.,
                    症例.EP前レジメン情報.中止に至った有害事象名.英語.,
                    症例.EP前レジメン情報.最悪Grade.名称.,
                    症例.EP前副作用情報.発症.覚知.日付,
                    症例.EP前副作用情報.CTCAEv5.0名称英語,
                    症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称.,
                    症例.EP後レジメン情報.提示された治療薬を投与しなかった理由.名称.,
                    症例.EP後レジメン情報.提示された治療薬とは異なる薬剤の投与有無.名称.,
                    症例.EP後レジメン情報.増悪確認日,
                    症例.EP後レジメン情報.Grade3以上有害事象の有無.名称.,
                    症例.EP後レジメン情報.中止に至った有害事象名.英語.,
                    症例.EP後レジメン情報.最悪Grade.名称.,
                    症例.EP後副作用情報.発症.覚知.日付,
                    症例.EP後副作用情報.CTCAEv5.0名称英語,
                    症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称.,
                    EP_treat_date
      )
    incProgress(1 / 4)
  })
  return(clin_tmp)
})

if (CCAT_FLAG & file.exists(file.path(app_dir, "source", "clinical_data_whole.qs"))) {
  initial_data_case <- QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(app_dir, "source", "clinical_data_whole.qs")) %>%
    dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
  initial_data_case = initial_data_case %>% dplyr::mutate(
    症例.背景情報.喫煙年数 = case_when(
      症例.背景情報.喫煙年数 == 99 ~ NA_integer_,
      症例.背景情報.喫煙歴有無.名称. %in% c("Unknown", "No") ~ NA_integer_,
      TRUE ~ 症例.背景情報.喫煙年数
    ))
  Data_case_raw <- reactive({ initial_data_case })
} else {
  Data_case_raw =  reactive({
    withProgress(message = "Clinical data loading.", {
      if(!is.null(input$new_analysis) && input$new_analysis =="No, use the previous dataset" &&
         file.exists(file.path(tempdir(), "clinical_data.qs"))){
        clin_tmp = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(tempdir(), "clinical_data.qs")) %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
      } else {
        clin_tmp = Data_case_raw_pre()
        # figure_drug_load_trigger(figure_drug_load_trigger() + 1)
        ID_pre_CGP_lines = Data_drug_raw() %>%
          dplyr::filter(TxCGP == "Pre") %>%
          dplyr::select(ID, CTx_line) %>%
          dplyr::arrange(desc(CTx_line)) %>%
          dplyr::distinct()
        colnames(ID_pre_CGP_lines) = c("ID", "CTx_lines_before_CGP")
        ID_pre_CGP_RECIST = Data_drug_raw() %>%
          dplyr::filter(TxCGP == "Pre") %>%
          dplyr::select(ID, RECIST) %>%
          dplyr::arrange(desc(RECIST)) %>%
          dplyr::distinct()
        colnames(ID_pre_CGP_RECIST) = c("ID", "RECIST")
        ID_pre_CGP_RECIST = ID_pre_CGP_RECIST %>% dplyr::mutate(RECIST = case_when(
          RECIST == "5" ~ "CR",
          RECIST == "4" ~ "PR",
          RECIST == "3" ~ "SD",
          RECIST == "2" ~ "PD",
          RECIST == "1" ~ "NE",
          RECIST == "0" ~ "NE",
          TRUE ~ "NE"
        ))
        ID_pre_CGP_1L_date = Data_drug_raw() %>%
          dplyr::filter(CTx_line == 1) %>%
          dplyr::select(ID, 投与開始日)
        colnames(ID_pre_CGP_1L_date) = c("ID", "start_1L")
        ID_pre_CGP_2L_date = Data_drug_raw() %>%
          dplyr::filter(CTx_line == 2) %>%
          dplyr::select(ID, 投与開始日)
        colnames(ID_pre_CGP_2L_date) = c("ID", "start_2L")
        clin_tmp$CTx_lines_before_CGP = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
          as.vector(ID_pre_CGP_lines$CTx_lines_before_CGP[match(x, ID_pre_CGP_lines$ID)])}))
        clin_tmp$pre_CGP_best_RECIST = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
          as.vector(ID_pre_CGP_RECIST$RECIST[match(x, ID_pre_CGP_RECIST$ID)])}))
        clin_tmp$Date_start_1L = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
          as.vector(ID_pre_CGP_1L_date$start_1L[match(x, ID_pre_CGP_1L_date$ID)])}))
        clin_tmp$Date_start_2L = unlist(lapply(list(clin_tmp$C.CAT調査結果.基本項目.ハッシュID), function(x) {
          as.vector(ID_pre_CGP_2L_date$start_2L[match(x, ID_pre_CGP_2L_date$ID)])}))
        clin_tmp$Date_start_1L = as.Date(clin_tmp$Date_start_1L)
        clin_tmp$Date_start_2L = as.Date(clin_tmp$Date_start_2L)
        clin_tmp$CTx_lines_before_CGP[is.na(clin_tmp$CTx_lines_before_CGP)] = "0"
        clin_tmp$pre_CGP_best_RECIST[is.na(clin_tmp$pre_CGP_best_RECIST)] = "CTx_naive"
        clin_tmp = clin_tmp %>%
          dplyr::mutate(time_enroll_treat =
                          as.integer(difftime(EP_treat_date,
                                              症例.管理情報.登録日,
                                              units = "days")),
                        time_enroll_final =
                          as.integer(difftime(final_observe,
                                              症例.管理情報.登録日,
                                              units = "days")),
                        time_palliative_final =
                          as.integer(difftime(final_observe,
                                              Date_start_1L,
                                              units = "days")),
                        time_2L_final =
                          as.integer(difftime(final_observe,
                                              Date_start_2L,
                                              units = "days")),
                        time_diagnosis_final =
                          as.integer(difftime(final_observe,
                                              症例.背景情報.診断日,
                                              units = "days")),
                        time_palliative_enroll =
                          time_palliative_final -
                          time_enroll_final,
                        time_diagnosis_enroll =
                          time_diagnosis_final -
                          time_enroll_final,
                        time_2L_enroll =
                          time_2L_final -
                          time_enroll_final
          )
        Data_tmp = Data_drug_raw() %>% dplyr::filter(
          TxCGP == "Post") %>%
          dplyr::distinct()
        ID_tmp = unique(Data_tmp$ID)
        clin_tmp$afterCGPtreat = clin_tmp$EP_treat
        clin_tmp$afterCGPtreat[clin_tmp$C.CAT調査結果.基本項目.ハッシュID %in% ID_tmp] <- 1
        clin_tmp = clin_tmp %>%
          dplyr::mutate(treat_group = case_when(
            afterCGPtreat == 1 & EP_treat == 1 ~ "Targetable mutation (+) and recommended treatment done",
            afterCGPtreat == 1 & EP_option == 1 ~ "Targetable mutation (+) and other treatment done",
            afterCGPtreat == 1 & EP_option == 0 ~ "Targetable mutation (-) and treatment done",
            afterCGPtreat == 0 & EP_option == 1 ~ "Targetable mutation (+) and no treatment done",
            afterCGPtreat == 0 & EP_option == 0 ~ "Targetable mutation (-) and no treatment done"
          ),
          treat_group_2 = case_when(
            afterCGPtreat == 1 & EP_treat == 1 ~ "Recommended treatment done",
            afterCGPtreat == 1 & EP_option == 1 ~ "Not-recommended treatment done",
            afterCGPtreat == 1 & EP_option == 0 ~ "Not-recommended treatment done",
            afterCGPtreat == 0 & EP_option == 1 ~ "No treatment done",
            afterCGPtreat == 0 & EP_option == 0 ~ "No treatment done"
          ),
          treat_group_3 = case_when(
            afterCGPtreat == 1 & EP_treat == 1 ~ "Treatment(+)",
            afterCGPtreat == 1 & EP_option == 1 ~ "Treatment(+)",
            afterCGPtreat == 1 & EP_option == 0 ~ "Treatment(+)",
            afterCGPtreat == 0 & EP_option == 1 ~ "Treatment(-)",
            afterCGPtreat == 0 & EP_option == 0 ~ "Treatment(-)"
          ))
        incProgress(1 / 4)
        clin_tmp = clin_tmp %>%
          dplyr::select(-症例.EP前レジメン情報.治療ライン.名称.,
                        -症例.EP前レジメン情報.実施目的.名称.,
                        -症例.EP前レジメン情報.化学療法レジメン名称,
                        -症例.EP前レジメン情報.薬剤名.YJ一般名.EN.,
                        -症例.EP前レジメン情報.投与開始日,
                        -症例.EP前レジメン情報.投与終了日,
                        -症例.EP前レジメン情報.終了理由.名称.,
                        -症例.EP前レジメン情報.レジメン継続区分.名称.,
                        -症例.EP前レジメン情報.最良総合効果.名称.,
                        -症例.EP後レジメン情報.治療ライン,
                        -症例.EP後レジメン情報.治療ライン.名称.,
                        -症例.EP後レジメン情報.化学療法レジメン名称,
                        -症例.EP後レジメン情報.薬剤名.YJ一般名.EN.,
                        -症例.EP後レジメン情報.投与開始日,
                        -症例.EP後レジメン情報.投与終了日,
                        -症例.EP後レジメン情報.終了理由.名称.,
                        -症例.EP後レジメン情報.レジメン継続区分.名称.,
                        -症例.EP後レジメン情報.最良総合効果.名称.,
                        -症例.EP前レジメン情報.増悪確認日,
                        -症例.EP前レジメン情報.Grade3以上有害事象の有無.名称.,
                        -症例.EP前レジメン情報.中止に至った有害事象名.英語.,
                        -症例.EP前レジメン情報.最悪Grade.名称.,
                        -症例.EP前副作用情報.発症.覚知.日付,
                        -症例.EP前副作用情報.CTCAEv5.0名称英語,
                        -症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称.,
                        -症例.EP後レジメン情報.提示された治療薬とは異なる薬剤の投与有無.名称.,
                        -症例.EP後レジメン情報.増悪確認日,
                        -症例.EP後レジメン情報.Grade3以上有害事象の有無.名称.,
                        -症例.EP後レジメン情報.中止に至った有害事象名.英語.,
                        -症例.EP後レジメン情報.最悪Grade.名称.,
                        -症例.EP後副作用情報.発症.覚知.日付,
                        -症例.EP後副作用情報.CTCAEv5.0名称英語,
                        -症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称.
          ) %>%
          dplyr::distinct()
        if(ENV_ != "server")
          QS_SAVE(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), clin_tmp, file=file.path(tempdir(), "clinical_data.qs"))
      }
    })
    return(clin_tmp)
  })
}

Data_case =  reactive({
  clin_tmp = Data_case_raw()
  if(input$histology_detail == "Yes, use oncotree 1st level"){
    clin_tmp$症例.基本情報.がん種.OncoTree. = clin_tmp$症例.基本情報.がん種.OncoTree.LEVEL1.
    clin_tmp$症例.基本情報.がん種.OncoTree..名称. = clin_tmp$症例.基本情報.がん種.OncoTree.LEVEL1.
  }
  if(!is.null(input$histology_group_1)){
    clin_tmp = clin_tmp %>% dplyr::mutate(
      症例.基本情報.がん種.OncoTree. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_1 ~ str_replace(str_replace(input$histology_group_1_name, "\\)", ""), ".*\\(",""),
        TRUE ~ 症例.基本情報.がん種.OncoTree.
      ),
      症例.基本情報.がん種.OncoTree..名称. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_1 ~ input$histology_group_1_name,
        TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
      )

    )
  }
  if(!is.null(input$histology_group_2)){
    clin_tmp = clin_tmp %>% dplyr::mutate(
      症例.基本情報.がん種.OncoTree. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_2 ~ str_replace(str_replace(input$histology_group_2_name, "\\)", ""), ".*\\(",""),
        TRUE ~ 症例.基本情報.がん種.OncoTree.
      ),
      症例.基本情報.がん種.OncoTree..名称. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_2 ~ input$histology_group_2_name,
        TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
      )
    )
  }
  if(!is.null(input$histology_group_3)){
    clin_tmp = clin_tmp %>% dplyr::mutate(
      症例.基本情報.がん種.OncoTree. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_3 ~ str_replace(str_replace(input$histology_group_3_name, "\\)", ""), ".*\\(",""),
        TRUE ~ 症例.基本情報.がん種.OncoTree.
      ),
      症例.基本情報.がん種.OncoTree..名称. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_3 ~ input$histology_group_3_name,
        TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
      )
    )
  }
  if(!is.null(input$histology_group_4)){
    clin_tmp = clin_tmp %>% dplyr::mutate(
      症例.基本情報.がん種.OncoTree. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_4 ~ str_replace(str_replace(input$histology_group_4_name, "\\)", ""), ".*\\(",""),
        TRUE ~ 症例.基本情報.がん種.OncoTree.
      ),
      症例.基本情報.がん種.OncoTree..名称. = case_when(
        症例.基本情報.がん種.OncoTree..名称. %in% input$histology_group_4 ~ input$histology_group_4_name,
        TRUE ~ 症例.基本情報.がん種.OncoTree..名称.
      )
    )
  }
  clin_tmp <- clin_tmp %>%
    # ① cancer name によるがん種名称の再設定
    {
      if (!is.null(input$minimum_pts) & input$minimum_pts != 1) {
        Cancername <- clin_tmp %>%
          distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.がん種.OncoTree.) %>%
          pull(症例.基本情報.がん種.OncoTree.)
        counts <- table(Cancername)
        Cancername <- names(counts[counts >= input$minimum_pts])
        mutate(.,
               症例.基本情報.がん種.OncoTree..名称. = if_else(症例.基本情報.がん種.OncoTree. %in% Cancername,
                                                   症例.基本情報.がん種.OncoTree..名称.,
                                                   症例.基本情報.がん種.OncoTree.LEVEL1.),
               症例.基本情報.がん種.OncoTree. = if_else(症例.基本情報.がん種.OncoTree. %in% Cancername,
                                               症例.基本情報.がん種.OncoTree.,
                                               症例.基本情報.がん種.OncoTree.LEVEL1.)
        )
      } else .
    } %>%
    # ② 各フィルター処理（NULL チェック付き）
    { if (!is.null(input$histology))
      filter(., 症例.基本情報.がん種.OncoTree..名称. %in% input$histology)
      else . } %>%
    { if (!is.null(input$panel))
      filter(., 症例.検体情報.パネル.名称. %in% input$panel)
      else . } %>%
    { if (!is.null(input$PS))
      filter(., 症例.背景情報.ECOG.PS.名称. %in% input$PS)
      else . } %>%
    { if (!is.null(input$sex))
      filter(., 症例.基本情報.性別.名称. %in% input$sex)
      else . } %>%
    { if (!is.null(input$stage))
      filter(., 症例.背景情報.初回治療前のステージ分類.名称. %in% input$stage)
      else . } %>%
    { if (!is.null(input$smoking))
      filter(., 症例.背景情報.喫煙歴有無.名称. %in% input$smoking)
      else . } %>%
    # ③ 登録日のフィルタリング：まず日付変換し、年フィルタ適用
    {
      if (!is.null(input$year)) {
        . <- mutate(., 症例.管理情報.登録日 = as.Date(症例.管理情報.登録日))
        filter(., (as.POSIXlt(症例.管理情報.登録日)$year + 1900) %in% as.integer(input$year) ) %>%
          mutate(., 症例.管理情報.登録日 = as.character(症例.管理情報.登録日))
      } else .
    } %>%
    # ④ 年齢フィルター（Data_case_raw() の最小・最大と比較して、入力範囲が変化している場合のみ）
    {
      if (!is.null(input$age)) {
        min_age <- min(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE)
        max_age <- max(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE)
        if (any(input$age != c(min_age, max_age))) {
          mutate(., 症例.基本情報.年齢 = replace_na(症例.基本情報.年齢, -1)) %>%
            filter(症例.基本情報.年齢 >= input$age[1] & 症例.基本情報.年齢 <= input$age[2])
        } else .
      } else .
    } %>%
    # ⑤ 年齢を基に若年/高齢フラグの追加
    { if (!is.null(input$mid_age))
      mutate(., YoungOld = if_else(症例.基本情報.年齢 <= input$mid_age, "Younger", "Older"))
      else .
    } %>% dplyr::distinct() %>% dplyr::arrange(C.CAT調査結果.基本項目.ハッシュID)
  clin_tmp$year = as.POSIXlt(as.Date(clin_tmp$症例.管理情報.登録日))$year + 1900
  return(clin_tmp)
})
