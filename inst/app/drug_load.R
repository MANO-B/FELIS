if (CCAT_FLAG & file.exists(file.path(app_dir, "source", "drug_data_whole.qs"))) {
  initial_data_drug <- QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(app_dir, "source", "drug_data_whole.qs")) %>%
    dplyr::filter(!ID %in% ID_exclude)
  Data_drug_raw <- reactive({ initial_data_drug })
} else {
  Data_drug_raw =  reactive({
    # req(figure_drug_load_trigger)
    withProgress(message = "Drug data loading.", {
      if(input$new_analysis == "No, use the previous dataset" &
         file.exists(file.path(tempdir(), "drug_data.qs"))){
        Data_drug = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(tempdir(), "drug_data.qs")) %>%
          dplyr::filter(!ID %in% ID_exclude)
      } else {
        clin_tmp = Data_case_raw_pre()
        tmp1 = clin_tmp %>%
          dplyr::filter(症例.EP前レジメン情報.実施目的.名称. %in%
                          c("その他", "緩和")) %>%
          dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                        症例.EP前レジメン情報.薬剤名.YJ一般名.EN.,
                        症例.EP前レジメン情報.化学療法レジメン名称,
                        症例.EP前レジメン情報.治療ライン.名称.,
                        症例.EP前レジメン情報.投与開始日,
                        症例.EP前レジメン情報.投与終了日,
                        症例.EP前レジメン情報.レジメン継続区分.名称.,
                        症例.管理情報.登録日,
                        症例.EP前レジメン情報.実施目的.名称.,
                        final_observe,
                        censor,
                        症例.EP前レジメン情報.終了理由.名称.,
                        症例.EP前レジメン情報.最良総合効果.名称.,
                        症例.EP前レジメン情報.増悪確認日,
                        症例.EP前副作用情報.発症.覚知.日付,
                        症例.EP前副作用情報.CTCAEv5.0名称英語,
                        症例.EP前副作用情報.CTCAEv5.0最悪Grade.名称.
          ) %>%
          dplyr::distinct() %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
        tmp1 = tmp1 %>%
          dplyr::mutate(TTD =
                          as.integer(difftime(final_observe,
                                              症例.EP前レジメン情報.投与開始日,
                                              units = "days")),
                        ToT =
                          as.integer(difftime(症例.EP前レジメン情報.投与終了日,
                                              症例.EP前レジメン情報.投与開始日,
                                              units = "days")),
                        TTF =
                          as.integer(difftime(症例.EP前レジメン情報.増悪確認日,
                                              症例.EP前レジメン情報.投与開始日,
                                              units = "days")),
                        TTAE =
                          as.integer(difftime(症例.EP前副作用情報.発症.覚知.日付,
                                              症例.EP前レジメン情報.投与開始日,
                                              units = "days")),
                        TTE =
                          as.integer(difftime(症例.管理情報.登録日,
                                              症例.EP前レジメン情報.投与開始日,
                                              units = "days")),
                        ToT_censor = 1,
                        TTF_censor = 1,
                        TTAE_censor = 1,
                        TTE_censor = 1
          )
        colnames(tmp1) = c("ID",
                           "Drug",
                           "レジメン",
                           "治療ライン",
                           "投与開始日",
                           "投与終了日",
                           "Ongoing",
                           "登録日",
                           "実施目的",
                           "final_observe",
                           "censor",
                           "終了理由",
                           "最良総合効果",
                           "増悪確認日",
                           "副作用発症日",
                           "副作用名称",
                           "副作用Grade",
                           "TTD",
                           "ToT",
                           "TTF",
                           "TTAE",
                           "TTE",
                           "ToT_censor",
                           "TTF_censor",
                           "TTAE_censor",
                           "TTE_censor"
        )
        tmp1$TxCGP = "Pre"
        tmp1 = tmp1 %>% dplyr::mutate(
          ToT_censor = case_when(
            is.na(投与終了日) & Ongoing == "継続中" & TxCGP =="Pre" ~ 0,
            TRUE ~ 1
          ),
          TTF_censor = case_when(
            is.na(増悪確認日) & Ongoing == "継続中" & TxCGP =="Pre" ~ 0,
            TRUE ~ 1
          ),
          ToT = case_when(
            is.na(投与終了日) & Ongoing == "継続中" & TxCGP =="Pre" ~ as.integer(difftime(登録日,
                                                                                  投与開始日,
                                                                                  units = "days")),
            TRUE ~ ToT
          ),
          TTF = case_when(
            is.na(増悪確認日) & Ongoing == "継続中" & TxCGP =="Pre" ~ as.integer(difftime(登録日,
                                                                                  投与開始日,
                                                                                  units = "days")),
            TRUE ~ TTF
          )
        )
        tmp2 = clin_tmp %>%
          dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                        症例.EP後レジメン情報.薬剤名.YJ一般名.EN.,
                        症例.EP後レジメン情報.化学療法レジメン名称,
                        症例.EP後レジメン情報.治療ライン.名称.,
                        症例.EP後レジメン情報.投与開始日,
                        症例.EP後レジメン情報.投与終了日,
                        症例.EP後レジメン情報.レジメン継続区分.名称.,
                        症例.管理情報.登録日,
                        症例.EP前レジメン情報.実施目的.名称.,
                        final_observe,
                        censor,
                        症例.EP後レジメン情報.終了理由.名称.,
                        症例.EP後レジメン情報.最良総合効果.名称.,
                        症例.EP後レジメン情報.増悪確認日,
                        症例.EP後副作用情報.発症.覚知.日付,
                        症例.EP後副作用情報.CTCAEv5.0名称英語,
                        症例.EP後副作用情報.CTCAEv5.0最悪Grade.名称.) %>%
          dplyr::distinct()
        tmp2$症例.EP前レジメン情報.実施目的.名称. = "緩和"
        tmp2 = tmp2 %>%
          dplyr::mutate(TTD =
                          as.integer(difftime(final_observe,
                                              症例.EP後レジメン情報.投与開始日,
                                              units = "days")),
                        ToT =
                          as.integer(difftime(症例.EP後レジメン情報.投与終了日,
                                              症例.EP後レジメン情報.投与開始日,
                                              units = "days")),
                        TTF =
                          as.integer(difftime(症例.EP後レジメン情報.増悪確認日,
                                              症例.EP後レジメン情報.投与開始日,
                                              units = "days")),
                        TTAE =
                          as.integer(difftime(症例.EP後副作用情報.発症.覚知.日付,
                                              症例.EP後レジメン情報.投与開始日,
                                              units = "days")),
                        TTE =
                          as.integer(difftime(症例.管理情報.登録日,
                                              症例.EP後レジメン情報.投与開始日,
                                              units = "days")),
                        ToT_censor = 1,
                        TTF_censor = 1,
                        TTAE_censor = 1,
                        TTE_censor = 1
          )
        colnames(tmp2) = c("ID",
                           "Drug",
                           "レジメン",
                           "治療ライン",
                           "投与開始日",
                           "投与終了日",
                           "Ongoing",
                           "登録日",
                           "実施目的",
                           "final_observe",
                           "censor",
                           "終了理由",
                           "最良総合効果",
                           "増悪確認日",
                           "副作用発症日",
                           "副作用名称",
                           "副作用Grade",
                           "TTD",
                           "ToT",
                           "TTF",
                           "TTAE",
                           "TTE",
                           "ToT_censor",
                           "TTF_censor",
                           "TTAE_censor",
                           "TTE_censor"
        )
        tmp2$TxCGP = "Post"
        Data_drug = rbind(tmp1, tmp2)

        if(!is.null(input$ID_drug)){
          Data_drug = data.frame(NULL)
          for(i in 1:length(input$ID_drug[,1])){
            tmp = read.csv(header = TRUE, file(input$ID_drug[[i, 'datapath']],
                                               encoding='UTF-8-BOM'))
            if("症例.EP前レジメン情報.投与開始日" %in% colnames(tmp)){
              colnames(tmp) = c("ID",
                                "登録日",
                                "治療ライン",
                                "実施目的",
                                "レジメン",
                                "Drug",
                                "投与開始日",
                                "投与終了日",
                                "Ongoing",
                                "終了理由",
                                "最良総合効果",
                                "最終生存確認日",
                                "死亡日",
                                "増悪確認日",
                                "副作用発症日",
                                "副作用名称",
                                "副作用Grade"
              )
              tmp$投与開始日 = as.Date(tmp$投与開始日)
              tmp$投与終了日 = as.Date(tmp$投与終了日)
              tmp = tmp %>%
                dplyr::mutate(
                  ToT = as.integer(difftime(tmp$投与終了日,
                                            tmp$投与開始日,
                                            units = "days")),
                  censor = case_when(
                    死亡日 != ""  & !is.na(死亡日) & 死亡日 != "           " ~ 1,
                    TRUE ~ 0
                  ),
                  final_observe = case_when(
                    最終生存確認日 != "" & 最終生存確認日 != "           " & !is.na(最終生存確認日) ~ 最終生存確認日,
                    死亡日 != "" & !is.na(死亡日) & 死亡日 != "           " ~ 死亡日,
                    投与終了日 != "" & !is.na(投与終了日) & 投与終了日 != "           " ~ 投与終了日,
                    投与開始日 != "" & !is.na(投与開始日) & 投与開始日 != "           " ~ 投与開始日,
                    登録日 != "" & !is.na(登録日) & 登録日 != "           " ~ 登録日,
                    TRUE ~ NA_character_))
              tmp$登録日 = as.Date(tmp$登録日)
              tmp$最終生存確認日 = as.Date(tmp$最終生存確認日)
              tmp$死亡日 = as.Date(tmp$死亡日)
              tmp$final_observe = as.Date(tmp$final_observe)
              tmp = tmp %>%
                dplyr::mutate(censor = ifelse(is.na(censor), 0, censor))
              tmp = tmp %>% dplyr::select(-最終生存確認日, -死亡日)
              tmp$TxCGP = "Pre"
              tmp = tmp %>%
                dplyr::mutate(TTD =
                                as.integer(difftime(final_observe,
                                                    投与開始日,
                                                    units = "days")),
                              ToT =
                                as.integer(difftime(投与終了日,
                                                    投与開始日,
                                                    units = "days")),
                              TTF =
                                as.integer(difftime(増悪確認日,
                                                    投与開始日,
                                                    units = "days")),
                              TTAE =
                                as.integer(difftime(副作用発症日,
                                                    投与開始日,
                                                    units = "days")),
                              TTE =
                                as.integer(difftime(登録日,
                                                    投与開始日,
                                                    units = "days")),
                              ToT_censor = 1,
                              TTF_censor = 1,
                              TTAE_censor = 1,
                              TTE_censor = 1
                )
              tmp = tmp %>% dplyr::mutate(
                ToT_censor = case_when(
                  is.na(投与終了日) & Ongoing == "継続中" & TxCGP =="Pre" ~ 0,
                  TRUE ~ 1
                ),
                TTF_censor = case_when(
                  is.na(増悪確認日) & Ongoing == "継続中" & TxCGP =="Pre" ~ 0,
                  TRUE ~ 1
                ),
                ToT = case_when(
                  is.na(投与終了日) & Ongoing == "継続中" & TxCGP =="Pre" ~ as.integer(difftime(登録日,
                                                                                        投与開始日,
                                                                                        units = "days")),
                  TRUE ~ ToT
                ),
                TTF = case_when(
                  is.na(増悪確認日) & Ongoing == "継続中" & TxCGP =="Pre" ~ as.integer(difftime(登録日,
                                                                                        投与開始日,
                                                                                        units = "days")),
                  TRUE ~ TTF
                )
              )
            } else{
              colnames(tmp) = c("ID",
                                "登録日",
                                "治療ライン",
                                "レジメン",
                                "Drug",
                                "投与開始日",
                                "投与終了日",
                                "Ongoing",
                                "終了理由",
                                "最良総合効果",
                                "最終生存確認日",
                                "死亡日",
                                "増悪確認日",
                                "副作用発症日",
                                "副作用名称",
                                "副作用Grade")
              tmp$投与開始日 = as.Date(tmp$投与開始日)
              tmp$投与終了日 = as.Date(tmp$投与終了日)
              tmp = tmp %>%
                dplyr::mutate(
                  ToT = as.integer(difftime(tmp$投与終了日,
                                            tmp$投与開始日,
                                            units = "days")),
                  censor = case_when(
                    死亡日 != ""  & !is.na(死亡日) & 死亡日 != "           " ~ 1,
                    TRUE ~ 0
                  ),
                  final_observe = case_when(
                    最終生存確認日 != "" & 最終生存確認日 != "           " & !is.na(最終生存確認日) ~ 最終生存確認日,
                    死亡日 != "" & !is.na(死亡日) & 死亡日 != "           " ~ 死亡日,
                    投与終了日 != "" & !is.na(投与終了日) & 投与終了日 != "           " ~ 投与終了日,
                    投与開始日 != "" & !is.na(投与開始日) & 投与開始日 != "           " ~ 投与開始日,
                    登録日 != "" & !is.na(登録日) & 登録日 != "           " ~ 登録日,
                    TRUE ~ NA_character_))
              tmp$登録日 = as.Date(tmp$登録日)
              tmp$最終生存確認日 = as.Date(tmp$最終生存確認日)
              tmp$死亡日 = as.Date(tmp$死亡日)
              tmp$final_observe = as.Date(tmp$final_observe)
              tmp = tmp %>%
                dplyr::mutate(censor = ifelse(is.na(censor), 0, censor))
              tmp = tmp %>% dplyr::select(-最終生存確認日, -死亡日)
              tmp$実施目的 = "緩和"
              tmp$TxCGP = "Post"
              tmp = tmp %>%
                dplyr::mutate(TTD =
                                as.integer(difftime(final_observe,
                                                    投与開始日,
                                                    units = "days")),
                              ToT =
                                as.integer(difftime(投与終了日,
                                                    投与開始日,
                                                    units = "days")),
                              TTF =
                                as.integer(difftime(増悪確認日,
                                                    投与開始日,
                                                    units = "days")),
                              TTAE =
                                as.integer(difftime(副作用発症日,
                                                    投与開始日,
                                                    units = "days")),
                              TTE =
                                as.integer(difftime(登録日,
                                                    投与開始日,
                                                    units = "days")),
                              ToT_censor = 1,
                              TTF_censor = 1,
                              TTAE_censor = 1,
                              TTE_censor = 1
                )
            }
            Data_drug = rbind(Data_drug, tmp)
          }
        }
        Data_drug = Data_drug %>% dplyr::mutate(RECIST = case_when(
          最良総合効果 == "CR" ~ "5",
          最良総合効果 == "PR" ~ "4",
          最良総合効果 == "SD" ~ "3",
          最良総合効果 == "PD" ~ "2",
          最良総合効果 == "NE" ~ "1",
          最良総合効果 == "Unknown" ~ "0",
          TRUE ~ "0"
        ))
        # 治療ライン → 数値に変換
        convert_line_to_num <- function(x) {
          dplyr::case_when(
            str_detect(x, "１次") ~ 1,
            str_detect(x, "２次") ~ 2,
            str_detect(x, "３次") ~ 3,
            str_detect(x, "４次") ~ 4,
            str_detect(x, "５次治療以降") ~ 5,
            TRUE ~ 99
          )
        }
        Data_drug = Data_drug %>%
          dplyr::filter(!is.na(治療ライン) | !is.na(投与開始日)) %>%
          dplyr::mutate(
            治療ライン数値 = convert_line_to_num(治療ライン),
            ToT_censor = case_when(
              ToT < 0 | is.na(ToT) ~ NA,
              TRUE ~ ToT_censor
            ),
            TTF_censor = case_when(
              TTF < 0 | is.na(TTF) ~ NA,
              TRUE ~ TTF_censor
            ),
            TTAE_censor = case_when(
              TTAE < 0 | is.na(TTAE) ~ NA,
              TRUE ~ TTAE_censor
            ),
            TTE_censor = case_when(
              TTE < 0 | is.na(TTE) ~ NA,
              TRUE ~ TTE_censor
            ),
            TTD = case_when(
              TTD < 0 | is.na(TTD) ~ NA,
              TRUE ~ TTD
            ),
            ToT = case_when(
              ToT < 0 | is.na(ToT) ~ NA,
              TRUE ~ ToT
            ),
            TTF = case_when(
              TTF < 0 | is.na(TTF) ~ NA,
              TRUE ~ TTF
            ),
            TTAE = case_when(
              TTAE < 0 | is.na(TTAE) ~ NA,
              TRUE ~ TTAE
            ),
            TTE = case_when(
              TTE < 0 | is.na(TTE) ~ NA,
              TRUE ~ TTE
            )
          )
        reorder_treatment_fast <- function(df) {
          dt <- as.data.table(df)

          # 日付の有無でフラグを設定
          dt[, has_date := !is.na(投与開始日)]

          # 日付あり・なしで分離
          dt_with_date <- dt[has_date == TRUE]
          dt_no_date <- dt[has_date == FALSE]

          # 日付ありデータを投与開始日でソート
          if (nrow(dt_with_date) > 0) {
            setorder(dt_with_date, 投与開始日)
          }

          # 日付なしデータを治療ライン数値でソート
          if (nrow(dt_no_date) > 0) {
            setorder(dt_no_date, 治療ライン数値)
          }

          if (nrow(dt_no_date) == 0) {
            # 日付なしデータがない場合
            result <- dt_with_date
          } else if (nrow(dt_with_date) == 0) {
            # 日付ありデータがない場合
            result <- dt_no_date
          } else {
            # 両方ある場合の高速マージ
            # 各日付なしレコードの挿入位置を計算
            insert_positions <- sapply(dt_no_date$治療ライン数値, function(line) {
              same_line_pos <- which(dt_with_date$治療ライン数値 == line)
              if (length(same_line_pos) > 0) {
                return(max(same_line_pos))
              } else {
                lower_pos <- which(dt_with_date$治療ライン数値 < line)
                return(if (length(lower_pos) > 0) max(lower_pos) else 0)
              }
            })

            # 挿入位置に基づいてソート順を決定
            dt_no_date[, insert_pos := insert_positions]
            dt_no_date[, sort_key := insert_pos + seq_len(nrow(dt_no_date)) / 1000]
            dt_with_date[, sort_key := seq_len(nrow(dt_with_date))]

            # 結合してソート
            result <- rbindlist(list(dt_with_date, dt_no_date), fill = TRUE)
            setorder(result, sort_key)
          }
          # CTx_lineを設定
          result[, CTx_line := seq_len(.N)]
          # 存在する列のみを削除
          cols_to_remove <- intersect(names(result), c("has_date", "insert_pos", "sort_key"))
          if (length(cols_to_remove) > 0) {
            result[, (cols_to_remove) := NULL]
          }
          return(as.data.frame(result))
        }

        # Classification function
        classify_ae <- function(ae_term, rules) {
          ae_lower <- tolower(ae_term)
          # If already contains "- Other, specify", extract the category part
          if (grepl("- Other, specify", ae_term)) {
            category <- str_extract(ae_term, "^[^-]+")
            category <- trimws(category)
            return(category)
          }
          # Match against each category's rules
          for (category in names(rules)) {
            keywords <- rules[[category]]
            for (keyword in keywords) {
              if (grepl(keyword, ae_lower, fixed = TRUE)) {
                return(category)
              }
            }
          }
          return(NA)
        }
        Data_drug <- Data_drug %>%
          mutate(
            副作用名称_clean = tolower(trimws(副作用名称)),
            # スペル修正マッピング
            副作用名称 = case_when(
              副作用名称_clean == "neutorophil count decreased" ~ "Neutrophil count decreased",
              副作用名称_clean == "neutoropenia" ~ "Neutropenia",
              副作用名称_clean == "neutophil count decreased" ~ "Neutrophil count decreased",
              副作用名称_clean == "proteiuria" ~ "Proteinuria",
              副作用名称_clean == "meningtits" ~ "Meningitis",
              副作用名称_clean == "alanine aminotranseferase" ~ "Alanine aminotransferase increased",
              副作用名称_clean == "myathenia gravis" ~ "Myasthenia gravis",
              副作用名称_clean == "thorombocytepenia" ~ "Thrombocytopenia",
              副作用名称_clean == "diarehhea|diarrea" ~ "Diarrhea",
              副作用名称_clean == "cardiomyocitis" ~ "Cardiomyopathy",
              副作用名称_clean == "bilarty tract infection" ~ "Biliary tract infection",
              副作用名称_clean == "ild" ~ "interstitial lung disease",
              副作用名称_clean == "rash acnelform" ~ "Rash acneiform",
              副作用名称_clean == "transeaminase" ~ "Transaminase elevation",
              副作用名称_clean == "ast、alt上昇|肝障害" ~ "Transaminase elevation",
              副作用名称_clean == "lymphocye count decreased" ~ "Lymphocyte count decreased",
              副作用名称_clean == "fn\\n" ~ "Febrile neutropenia",
              副作用名称_clean == "neutrophil count decrease" ~ "Neutrophil count decreased",
              TRUE ~ 副作用名称
            )
          ) %>%
          dplyr::select(-副作用名称_clean)
        Data_drug <- Data_drug %>%
          mutate(
            CTCAE = sapply(副作用名称, function(x) classify_ae(x, classification_rules)),
            Adverse_effect = ifelse(!is.na(副作用名称), "AE (+)", "AE (-)"),
            副作用Grade = case_when(
              Adverse_effect == "AE (-)" ~ NA,
              TRUE ~ 副作用Grade
            ),
            Overall_AE_grage = 副作用Grade
          )
        # 全体に適用（IDごと）
        Data_drug_unique <- Data_drug %>%
          dplyr::distinct(ID, レジメン, 治療ライン数値, 投与開始日, .keep_all = T) %>%
          group_split(ID) %>%
          map_df(reorder_treatment_fast) %>%
          ungroup()
        for (soc in names(classification_rules)) {
          # Drugを整形したData_drug
          Data_drug_clean <- Data_drug %>%
            dplyr::select(ID, レジメン, 治療ライン数値, 投与開始日, CTCAE, 副作用Grade) %>%
            dplyr::filter(!is.na(CTCAE)) %>%
            dplyr::filter(CTCAE == soc) %>%
            dplyr::arrange(desc(副作用Grade)) %>%
            dplyr::distinct(ID, レジメン, 治療ライン数値, 投与開始日, CTCAE,.keep_all = T)
          Data_drug_clean[[soc]] = Data_drug_clean$副作用Grade
          Data_drug_clean = Data_drug_clean %>%
            dplyr::select(-副作用Grade)
          # 左結合してDrug列を追加
          Data_drug_unique <- Data_drug_unique %>%
            left_join(Data_drug_clean, by = c("ID", "レジメン", "治療ライン数値", "投与開始日", "CTCAE"))
          Data_drug_unique[[soc]] = ifelse(Data_drug_unique[[soc]] == "", "AE (-)", Data_drug_clean[[soc]])
          Data_drug_unique[[soc]] = ifelse(is.na(Data_drug_unique[[soc]]), "AE (-)", Data_drug_clean[[soc]])
        }
        # Drugを整理
        Data_drug_clean <- Data_drug %>%
          dplyr::select(ID, レジメン, 治療ライン数値, 投与開始日, Drug) %>%
          filter(!is.na(Drug) & Drug !="") %>%
          group_by(ID, レジメン, 治療ライン数値, 投与開始日) %>%
          summarise(Drug = paste(sort(unique(Drug)), collapse = ","), .groups = "drop")
        Data_drug_unique <- Data_drug_unique %>%
          dplyr::select(-Drug) %>%
          left_join(Data_drug_clean, by = c("ID", "レジメン", "治療ライン数値", "投与開始日")) %>%
          dplyr::mutate(Drug = ifelse(Drug == "", NA, Drug))
        # AEの有無を整理
        Data_drug_clean <- Data_drug %>%
          dplyr::select(ID, レジメン, 治療ライン数値, 投与開始日, Adverse_effect) %>%
          filter(!is.na(Adverse_effect) & Adverse_effect =="AE (+)") %>%
          distinct(ID, レジメン, 治療ライン数値, 投与開始日, .keep_all = T)
        Data_drug_unique <- Data_drug_unique %>%
          dplyr::select(-Adverse_effect) %>%
          left_join(Data_drug_clean, by = c("ID", "レジメン", "治療ライン数値", "投与開始日")) %>%
          dplyr::mutate(Adverse_effect = ifelse(is.na(Adverse_effect), "AE (-)", Adverse_effect))
        # 最高Gradeを整理
        Data_drug_clean <- Data_drug %>%
          dplyr::select(ID, レジメン, 治療ライン数値, 投与開始日, Overall_AE_grage) %>%
          filter(!is.na(Overall_AE_grage) & Overall_AE_grage !="") %>%
          dplyr::arrange(desc(Overall_AE_grage)) %>%
          distinct(ID, レジメン, 治療ライン数値, 投与開始日, .keep_all = T)
        Data_drug_unique <- Data_drug_unique %>%
          dplyr::select(-Overall_AE_grage) %>%
          left_join(Data_drug_clean, by = c("ID", "レジメン", "治療ライン数値", "投与開始日")) %>%
          dplyr::mutate(Overall_AE_grage = ifelse(Overall_AE_grage == "", NA, Overall_AE_grage))
        # 副作用発症日を整理
        Data_drug_clean <- Data_drug %>%
          dplyr::select(ID, レジメン, 治療ライン数値, 投与開始日, 副作用発症日) %>%
          dplyr::filter(!is.na(副作用発症日)) %>%
          dplyr::arrange(副作用発症日) %>%
          dplyr::distinct(ID, レジメン, 治療ライン数値, 投与開始日, .keep_all = T)
        Data_drug_unique = Data_drug_unique %>% dplyr::select(-副作用発症日)
        Data_drug_unique <- Data_drug_unique %>%
          left_join(Data_drug_clean, by = c("ID", "レジメン", "治療ライン数値", "投与開始日"))
        setDT(Data_drug_unique)
        # 結合とTTNT計算を一度に実行
        Data_drug <- as.data.table(Data_drug_unique[, TTNT := as.numeric(shift(投与開始日, type = "lead") - 投与開始日), by = ID])
        Data_drug = Data_drug %>% dplyr::mutate(
          最良総合効果 = case_when(
            RECIST == "5" ~ "CR",
            RECIST == "4" ~ "PR",
            RECIST == "3" ~ "SD",
            RECIST == "2" ~ "PD",
            RECIST == "1" ~ "NE",
            RECIST == "0" ~ "NE",
            TRUE ~ "NE"
          )
        ) %>%
          dplyr::select(-副作用名称, -CTCAE)
        incProgress(1 / 4)
        QS_SAVE(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), Data_drug, file=file.path(tempdir(), "drug_data.qs"))
      }
    })
    return(Data_drug)
  })
}

Data_drug_raw_rename =  reactive({
  Data_drug = Data_drug_raw()
  if(!is.null(input$drug_combination_rename)){
    drug_combination_rename_list = data.frame(NULL)
    for(i in 1:length(input$drug_combination_rename[,1])){
      drug_combination_rename_list <- rbind(drug_combination_rename_list,
                                            read.csv(header = TRUE,
                                                     file(input$drug_combination_rename[[i, 'datapath']],
                                                          encoding='UTF-8-BOM')))
    }
    Data_drug$Drug_tmp = unlist(lapply(list(Data_drug$Drug), function(x) {
      as.vector(drug_combination_rename_list$Rename[match(x, drug_combination_rename_list$Drug)])}))
    Data_drug = Data_drug %>% dplyr::mutate(
      Drug = case_when(
        !is.na(Drug_tmp) ~ Drug_tmp,
        TRUE ~ Drug
      )
    )
    Data_drug = Data_drug %>% dplyr::select(-Drug_tmp)
  }
  Data_drug$CTx_line = as.character(Data_drug$CTx_line)
  Data_drug = Data_drug %>% dplyr::mutate(
    CTx_line = case_when(
      CTx_line %in% c("1", "2", "3", "4") ~ CTx_line,
      TRUE ~ "5~"
    )
  )
  return(Data_drug)
})
