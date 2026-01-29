survival_CGP_analysis_logic <- function() {
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
                      time_enroll_final,
                      censor,
                      CTx_lines_before_CGP,
                      pre_CGP_best_RECIST,
                      treat_group,
                      treat_group_2,
                      treat_group_3,
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
                      time_enroll_treat,
                      year
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
      Data_cluster_ID_list = Data_cluster_ID() %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster)
      Data_case_target = left_join(Data_case_target,
                                   Data_cluster_ID_list,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_case_target$cluster[is.na(Data_case_target$cluster)] = max(Data_case_target$cluster, na.rm = T) + 1
      Data_drug = Data_drug_raw()
      OUTPUT_DATA$figure_surv_CGP_Data_drug = Data_drug
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

      Data_report_TMB = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
        dplyr::select(Tumor_Sample_Barcode, TMB) %>%
        dplyr::distinct(Tumor_Sample_Barcode, .keep_all=T)
      colnames(Data_report_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB")
      Data_case_target = left_join(Data_case_target, Data_report_TMB,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")

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
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
        dplyr::select(Tumor_Sample_Barcode, Evidence_level)
      colnames(Data_Best_Evidence_Level) = c("C.CAT調査結果.基本項目.ハッシュID", "Best_Evidence_Level")
      Data_case_target = left_join(Data_case_target, Data_Best_Evidence_Level,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_case_target$Best_Evidence_Level[is.na(Data_case_target$Best_Evidence_Level)] = "None"

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

      Data_survival = Data_case_target %>% dplyr::filter(
        !is.na(time_enroll_final) &
          is.finite(time_enroll_final) &
          time_enroll_final > 0 &
          !is.na(censor) &
          is.finite(censor)
      )
      Data_MAF_target = Data_MAF_target %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_survival$C.CAT調査結果.基本項目.ハッシュID)
      OUTPUT_DATA$figure_surv_CGP_Data_MAF_target = Data_MAF_target
      incProgress(1 / 13)

      Data_survival_interactive = Data_survival
      OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive = Data_survival_interactive
      candidate_genes = sort(unique(c(Data_MAF_target$Hugo_Symbol,
                                      paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                                      paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                                      paste0(input$special_gene, "_NOS"))))
      candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
      candidate_genes = candidate_genes[!is.na(candidate_genes)]
      Top_gene = unique(names(sort(table(Data_MAF_target$Hugo_Symbol), decreasing = T)))
      OUTPUT_DATA$figure_surv_CGP_Top_gene = Top_gene[Top_gene %in% candidate_genes]
      OUTPUT_DATA$figure_surv_CGP_candidate_genes = candidate_genes

      candidate_drugs = sort(unique(c(Data_drug$Drug)))
      OUTPUT_DATA$figure_surv_CGP_candidate_drugs = candidate_drugs[!is.na(candidate_drugs)]
      OUTPUT_DATA$figure_surv_CGP_candidate_lines = sort(unique(Data_survival_interactive$CTx_lines_before_CGP))
      OUTPUT_DATA$figure_surv_CGP_candidate_RECIST = sort(unique(Data_survival_interactive$pre_CGP_best_RECIST))
      OUTPUT_DATA$figure_surv_CGP_candidate_Age = sort(unique(Data_survival_interactive$YoungOld))
      OUTPUT_DATA$figure_surv_CGP_candidate_Sex = sort(unique(Data_survival_interactive$症例.基本情報.性別.名称.))
      OUTPUT_DATA$figure_surv_CGP_candidate_Histology = sort(unique(Data_survival_interactive$Cancers))
      OUTPUT_DATA$figure_surv_CGP_candidate_cluster = sort(unique(Data_survival_interactive$cluster))
      OUTPUT_DATA$figure_surv_CGP_candidate_PS = sort(unique(Data_survival_interactive$症例.背景情報.ECOG.PS.名称.))
      OUTPUT_DATA$figure_surv_CGP_candidate_Panel = sort(unique(Data_survival_interactive$症例.検体情報.パネル.名称.))
      OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level = sort(unique(Data_survival_interactive$Best_Evidence_Level))
      OUTPUT_DATA$figure_surv_CGP_candidate_Year = sort(unique(Data_survival_interactive$year))
      OUTPUT_DATA$figure_surv_CGP_candidate_meta = c('Lymph_met','Brain_met','Lung_met','Bone_met','Liver_met')

      incProgress(1 / 13)

      oncogenic_genes = Data_MAF_target %>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        dplyr::distinct() %>%
        dplyr::count(Hugo_Symbol) %>%
        dplyr::arrange(-n)
      colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
      oncogenic_genes = rbind(oncogenic_genes %>% dplyr::filter(gene_mutation %in% input$gene),
                              oncogenic_genes %>% dplyr::filter(!gene_mutation %in% input$gene))
      oncogenic_genes = oncogenic_genes[1:min(length(oncogenic_genes$all_patients), input$gene_no),]
      OUTPUT_DATA$figure_surv_CGP_oncogenic_genes = oncogenic_genes
      incProgress(1 / 13)

      Data_forest = Data_case_target %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                      YoungOld,
                      Lymph_met,
                      Brain_met,
                      Lung_met,
                      Bone_met,
                      Liver_met,
                      #Other_met,
                      EP_option,
                      EP_treat,
                      症例.基本情報.性別.名称.,
                      症例.基本情報.がん種.OncoTree.,
                      症例.基本情報.がん種.OncoTree.LEVEL1.,
                      症例.背景情報.ECOG.PS.名称.,
                      #症例.背景情報.喫煙歴有無.名称.,
                      #症例.背景情報.アルコール多飲有無.名称.,
                      time_enroll_final,
                      censor,
                      CTx_lines_before_CGP,
                      pre_CGP_best_RECIST,
                      症例.検体情報.パネル.名称.,
                      症例.基本情報.年齢,
                      cluster)
      colnames(Data_forest) =
        c('ID', 'Age',
          'Lymph_met',
          'Brain_met',
          'Lung_met',
          'Bone_met',
          'Liver_met',
          #'Other_met',
          'EP_option',
          'EP_treat',
          'Sex', 'Histology', 'Histology_dummy', 'PS',
          #'Smoking_history', 'Alcoholic_history',
          'time_enroll_final', 'censor', "Lines", "Best_effect", "Panel","Age_raw","Cluster")
      #Data_forest$Histology_dummy = paste0(Data_forest$Histology_dummy,"_NOS")
      col_added = 0
      if(input$HER2 != "No"){
        Data_forest$HER2_IHC = Data_case_target$HER2_IHC
        Data_forest = Data_forest %>%
          dplyr::filter(
            HER2_IHC != "Unknown"
          )
        col_added = col_added + 1
      }
      if(input$MSI != "No"){
        Data_forest$MSI_PCR = Data_case_target$MSI_PCR
        Data_forest = Data_forest %>%
          dplyr::filter(
            MSI_PCR != "Unknown"
          )
        col_added = col_added + 1
      }
      if(input$MMR != "No"){
        Data_forest$MMR_IHC = Data_case_target$MMR_IHC
        Data_forest = Data_forest %>%
          dplyr::filter(
            MMR_IHC != "Unknown"
          )
        col_added = col_added + 1
      }

      Data_forest = Data_forest %>%
        dplyr::filter(PS != "Unknown" & Sex != "Unknown") %>%
        dplyr::mutate(
          PS = case_when(
            PS == 0 ~ "0",
            PS == 1 ~ "1",
            TRUE ~ "2_4"
          ),
          Lines = case_when(
            Lines == "0" ~ "0",
            Lines == "1" ~ "1",
            TRUE ~ "2~"
          ),
          Panel = case_when(
            Panel %in% c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP") ~ "Solid",
            TRUE ~ "Liquid"
          )
        )
      Data_forest$Panel = factor(Data_forest$Panel)
      Data_forest$Cluster = factor(Data_forest$Cluster)
      Data_forest$Best_effect[Data_forest$Best_effect == "Unknown"] = "NE"
      Data_forest$Best_effect[Data_forest$Best_effect == "NE"] = "NE"
      Data_forest$Best_effect[is.na(Data_forest$Best_effect)] = "NE"
      Data_forest$Best_effect = factor(Data_forest$Best_effect, levels = c("SD","CR","PR","PD","CTx_naive","NE"))
      if(length(unique(Data_forest$Panel))>1){
        Data_forest$Panel =  relevel(Data_forest$Panel, ref="Solid")
      }
      frequent_organs = sort(table(Data_forest$Histology),decreasing = T)
      frequent_organs_names = names(frequent_organs[frequent_organs>=50])
      Data_forest = Data_forest %>% dplyr::mutate(
        Histology = case_when(
          Histology %in% frequent_organs_names ~ Histology,
          TRUE ~ Histology_dummy
        )
      )
      Data_ML =  Data_forest %>%
        dplyr::select(-Age, -time_enroll_final, -censor, -EP_option, -ID, -Histology_dummy)
      Data_forest =  Data_forest %>%
        dplyr::select(-Age_raw)
      gene_names = unique(c(input$gene,
                            as.character(oncogenic_genes$gene_mutation)))
      gene_names = gene_names[gene_names %in% unique(Data_MAF_target$Hugo_Symbol)]

      gene_names = gene_names[1:min(12, length(gene_names))]
      OUTPUT_DATA$figure_surv_CGP_gene_names = gene_names
      # 1. 全遺伝子の変異サンプルIDを一括取得
      mutation_lists <- Data_MAF_target %>%
        filter(Hugo_Symbol %in% gene_names) %>%
        group_by(Hugo_Symbol) %>%
        summarise(ID_mutation = list(unique(Tumor_Sample_Barcode)), .groups = 'drop')
      # 2. Data_forestの一意IDを事前取得
      unique_ids <- Data_forest %>%
        select(ID) %>%
        distinct()
      # 4. 全遺伝子の変異ステータスを一括作成（NA対応版）
      mutation_matrix <- sapply(gene_names, function(gene) {
        # 該当遺伝子の変異リストを取得
        gene_row <- mutation_lists[mutation_lists$Hugo_Symbol == gene, ]
        if(nrow(gene_row) == 0) {
          # 遺伝子がMAFデータに存在しない場合
          rep("mut(-)", nrow(unique_ids))
        } else {
          current_mutations <- gene_row$ID_mutation[[1]]
          if(length(current_mutations) == 0) {
            # 変異IDが空の場合
            rep("mut(-)", nrow(unique_ids))
          } else {
            # 正常な場合の変異ステータス判定
            ifelse(unique_ids$ID %in% current_mutations, "mut(+)", "mut(-)")
          }
        }
      })
      # 5. 列名設定
      colnames(mutation_matrix) <- gene_names
      # 6. データフレーム化（文字列として保持）
      mutation_df <- data.frame(unique_ids, mutation_matrix, stringsAsFactors = FALSE)
      # 7. 一括結合
      Data_forest <- left_join(Data_forest, mutation_df, by = "ID")

      Factor_names = colnames(Data_forest)
      Factor_names = Factor_names[!Factor_names %in% c('ID', 'time_enroll_final', 'censor', 'EP_treat')]

      Data_forest_tmp = Data_forest
      gene_cols <- (19 + col_added + 1):(19 + col_added + length(gene_names))
      mut_counts <- sapply(gene_cols, function(col) {
        sum(Data_forest_tmp[[col]] == "mut(+)" & Data_forest_tmp$censor == 1)
      })
      non_mut_counts <- sapply(gene_cols, function(col) {
        sum(Data_forest_tmp[[col]] != "mut(+)" & Data_forest_tmp$censor == 1)
      })
      genes_to_remove <- which(mut_counts < 3 | non_mut_counts < 3)
      if(length(genes_to_remove) > 0) {
        remove_indices <- 15 + col_added + genes_to_remove
        Factor_names <- Factor_names[-remove_indices]
      }
      disease_counts <- Data_forest_tmp %>%
        group_by(Histology) %>%
        summarise(censor_count = sum(censor == 1), .groups = 'drop') %>%
        filter(censor_count < 3)
      if(nrow(disease_counts) > 0) {
        Data_forest_tmp <- Data_forest_tmp %>%
          mutate(Histology = ifelse(Histology %in% disease_counts$Histology,
                                    Histology_dummy,
                                    Histology))
      }
      Data_forest_tmp$Histology = factor(Data_forest_tmp$Histology)
      setDT(Data_forest_tmp)

      # 再レベル設定（Histology, Cluster）
      for (col in c("Histology", "Cluster")) {
        if (length(unique(Data_forest_tmp[[col]])) < 2) {
          Factor_names <- setdiff(Factor_names, col)
        } else {
          Data_forest_tmp[[col]] <- relevel(Data_forest_tmp[[col]], ref = names(sort(table(Data_forest_tmp[[col]]), decreasing = TRUE))[1])
        }
      }

      # 2値カテゴリ用: レベル数チェック or イベント数チェック
      check_binary <- function(var, level1) {
        s1 <- sum(Data_forest_tmp[[var]] == level1 & Data_forest_tmp$censor == 1, na.rm = TRUE)
        s2 <- sum(Data_forest_tmp[[var]] != level1 & Data_forest_tmp$censor == 1, na.rm = TRUE)
        if (s1 < 3 || s2 < 3) {
          Factor_names <<- setdiff(Factor_names, var)
        }
      }

      # EP_option, Age, Sex, Panel, PS, Lines, Best_effect, 各転移
      check_binary("EP_option", 1)
      check_binary("Age", "Older")
      check_binary("Sex", "Male")
      check_binary("Panel", "Solid")
      check_binary("PS", "0")
      check_binary("Lines", "1")
      check_binary("Best_effect", "SD")
      for (met in c("Lymph_met", "Lung_met", "Brain_met", "Bone_met", "Liver_met")) {
        check_binary(met, "Yes")
      }

      # 変換（条件付き）
      if (sum(Data_forest_tmp$PS == "2_4" & Data_forest_tmp$censor == 1) < 3) {
        Data_forest_tmp[PS %in% c("1", "2_4"), Lines := "1_4"]
      }
      if (sum(Data_forest_tmp$Lines == "0" & Data_forest_tmp$censor == 1) < 3) {
        Data_forest_tmp[Lines %in% c("0", "1"), Lines := "0~1"]
      }
      if (sum(Data_forest_tmp$Lines == "2~" & Data_forest_tmp$censor == 1) < 3) {
        Data_forest_tmp[Lines %in% c("1", "2~"), Lines := "1~"]
      }

      # HER2 / MSI / MMR チェック（Shiny input付き）
      if (input$HER2 != "No") check_binary("HER2_IHC", "Positive")
      if (input$MSI != "No") check_binary("MSI_PCR", "Positive")
      if (input$MMR != "No") check_binary("MMR_IHC", "Positive")

      setDT(Data_forest_tmp)

      # ----- 同じ値しか持たない列を除去 -----
      Factor_names <- Factor_names[
        vapply(Factor_names, function(x) uniqueN(Data_forest_tmp[[x]]) > 1, logical(1))
      ]

      # Histology_dummy を除外
      Factor_names <- setdiff(Factor_names, "Histology_dummy")

      # 保存
      Factor_names_univariant <- Factor_names

      # ----- 高相関列の除去 -----
      if (length(Factor_names) > 1) {
        # 因子型にして数値化（corに使うため）
        dt_factor <- Data_forest_tmp[, lapply(.SD, function(x) as.integer(as.factor(x))),
                                     .SDcols = Factor_names]

        # 相関行列
        cor_mat <- cor(dt_factor, use = "pairwise.complete.obs")

        # 上三角行列だけを見る
        cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
        too_high <- which(abs(cor_mat) > 0.95, arr.ind = TRUE)

        # 除外対象列（列番号で後ろの列）
        to_remove <- unique(colnames(cor_mat)[too_high[, "col"]])

        Factor_names <- setdiff(Factor_names, to_remove)
      }
      incProgress(1 / 13)

      OUTPUT_DATA$figure_surv_CGP_Data_forest_tmp = Data_forest_tmp
      OUTPUT_DATA$figure_surv_CGP_candidate_EP_option = sort(unique(Data_survival_interactive$EP_option))
      OUTPUT_DATA$figure_surv_CGP_candidate_EP_treat = sort(unique(Data_survival_interactive$treat_group_2))
      OUTPUT_DATA$figure_surv_CGP_candidate_Panel = sort(unique(Data_survival_interactive$症例.検体情報.パネル.名称.))
      OUTPUT_DATA$figure_surv_CGP_Factor_names = Factor_names
      OUTPUT_DATA$figure_surv_CGP_Factor_names_univariant = Factor_names_univariant
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}

output$figure_survival_treatment_reach_1 = renderGirafe({
  Data_clean <- ID_select_surv_CGP()
  req(Data_clean)
  if(nrow(Data_clean) > 0 & sum(Data_clean$EP_treat) > 0){
    # 競合リスク解析（cumulative incidence function）
    cif_result <- cuminc(
      ftime = Data_clean$time_to_event/365.25*12,
      fstatus = Data_clean$event,
      cencode = 0  # 打ち切りコード
    )
    # 治療到達の累積発生率を取得
    cif_treat <- data.frame(
      time = cif_result$`1 1`$time,
      prob = cif_result$`1 1`$est
    )
    # 同一時点のデータを集約して処理
    cif_treat_clean <- cif_treat %>%
      arrange(time) %>%
      group_by(time) %>%
      summarise(prob = last(prob), .groups = 'drop') %>%  # 同一時点では最後の値を使用
      arrange(time)
    # 瞬間治療到達率の計算
    cif_treat <- cif_treat_clean %>%
      mutate(
        time_diff = c(1, diff(time)),  # 時間差（最初は1日とする）
        prob_diff = c(prob[1], diff(prob)),  # 確率差
        # 瞬間治療到達率（ハザード様の概念）
        instant_rate = ifelse(time_diff > 0,
                              prob_diff / time_diff * 100,  # %/日単位
                              0),  # 時間差が0なら0
        instant_rate = pmax(instant_rate, 0)  # 負の値を0に
      )
    # 移動平均の計算
    window_size <- 3
    # 十分なデータポイントがある場合のみ移動平均を計算
    if(nrow(cif_treat) >= window_size) {
      cif_treat <- cif_treat %>%
        mutate(
          ma_instant_rate = rollmean(instant_rate, k = min(window_size, nrow(cif_treat)),
                                     fill = NA, align = "center"),
        )
    } else {
      cif_treat <- cif_treat %>%
        mutate(
          ma_instant_rate = instant_rate,
        )
    }
    # より滑らかな移動平均のための局所回帰（LOESS）
    if(nrow(cif_treat) > 10) {
      loess_fit <- loess(instant_rate ~ time, data = cif_treat, span = 0.3)
      cif_treat$smooth_rate <- predict(loess_fit)
    } else {
      cif_treat$smooth_rate <- cif_treat$instant_rate
    }
    # LOESS平滑化曲線からの最大値算出
    max_instant_idx <- which.max(cif_treat$smooth_rate)
    max_instant_rate_loess <- cif_treat$smooth_rate[max_instant_idx]
    max_instant_time <- cif_treat$time[max_instant_idx]
    treat_max <- max(cif_result$`1 1`$est)
    treat_time_max <- cif_result$`1 1`$time[which.max(cif_result$`1 1`$est)]
    # 可視化（改良版）
    p1 <- ggplot(cif_treat, aes(x = time)) +
      geom_line_interactive(aes(y = instant_rate, color = "instantaneous rate"),
                            alpha = 0.6, size = 0.5) +
      geom_point_interactive(
        aes(
          y = instant_rate,
          tooltip = paste0("time: ", round(time, 1), " months\ninstant rate: ", round(instant_rate, 1), "%")
        ),
        size = 0, alpha = 0
      ) +
      geom_line_interactive(aes(y = ma_instant_rate, color = "moving average (3-month window)"),
                            alpha = 0.6, size = 1) +
      geom_point_interactive(
        aes(
          y = ma_instant_rate,
          tooltip = paste0("time: ", round(time, 1), " months\ninstant rate: ", round(ma_instant_rate, 1), "%")
        ),
        size = 0, alpha = 0
      ) +
      geom_line_interactive(aes(y = smooth_rate, color = "LOESS smoothing"),
                            alpha = 1, size = 1.2) +
      geom_point_interactive(
        aes(
          y = smooth_rate,
          tooltip = paste0("time: ", round(time, 1), " months\ninstant rate: ", round(smooth_rate, 1), "%")
        ),
        size = 0, alpha = 0
      ) +
      geom_point_interactive(aes(x = max_instant_time, y = max_instant_rate_loess),
                             color = "red", size = 3, shape = 21, fill = "white", stroke = 2) +
      geom_text_interactive(aes(x = max_instant_time * 1.4, y = max_instant_rate_loess),
                            label = paste0("Max: ", round(max_instant_rate_loess, 1), "%/month, at ",
                                           sprintf("%.1f", max_instant_time), " month"),
                            vjust = 0.5, hjust = -0.5, size = 3.5, color = "red", fontface = "bold") +
      scale_color_manual_interactive(values = c("instantaneous rate" = "lightgray",
                                                "moving average (3-month window)" = "blue",
                                                "LOESS smoothing" = "red")) +
      labs(
        title = "Recommended treatment receiving rate (considering competitive risk)",
        subtitle = "Derivative of cumulative incidence rate",
        x = "Months after CGP test",
        y = "Treatment reach rate (%/month)",
        color = "Smoothing method",
        caption = "Competition risk: End of observation due to death"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "bottom"
      ) +
      scale_x_continuous(breaks = seq(0, max(6, cif_treat$time, na.rm = TRUE),
                                      by = max(6, ceiling(max(cif_treat$time)/60)*6)))
    ggiraph::girafe(ggobj = p1, width_svg = 8, height_svg = 8)
  }
})

output$figure_survival_treatment_reach_2 = renderGirafe({
  Data_clean <- ID_select_surv_CGP()
  req(Data_clean)
  if(nrow(Data_clean) > 0 & sum(Data_clean$EP_treat) > 0){
    # 競合リスク解析（cumulative incidence function）
    cif_result <- cuminc(
      ftime = Data_clean$time_to_event/365.25*12,
      fstatus = Data_clean$event,
      cencode = 0  # 打ち切りコード
    )
    # 治療到達の累積発生率を取得
    cif_treat <- data.frame(
      time = cif_result$`1 1`$time,
      prob = cif_result$`1 1`$est
    )
    treat_max <- max(cif_result$`1 1`$est)
    treat_time_max <- cif_result$`1 1`$time[which.max(cif_result$`1 1`$est)]
    # 可視化（改良版）
    # 累積発生率も併せて表示
    p2 <- ggplot(cif_treat) +
      geom_line_interactive(
        aes(
          x = time,
          y = prob * 100
        ),
        color = "darkgreen", size = 1.2) +
      geom_point_interactive(
        aes(
          x = time,
          y = prob * 100,
          tooltip = paste0("time: ", round(time, 1), " months\nCumulative rate: ", round(prob * 100, 1), "%")
        ),
        size = 0, alpha = 0
      ) +
      geom_point_interactive(aes(x = treat_time_max, y = treat_max * 100),
                             color = "red", size = 3, shape = 21, fill = "white", stroke = 2) +
      geom_text_interactive(aes(x = treat_time_max * 0.95, y = treat_max * 70),
                            label = paste0("Max: ", round(treat_max * 100, 1), "%\nat ",
                                           sprintf("%.1f", treat_time_max), " month"),
                            vjust = -1.5, hjust = 0.5, size = 3.5, color = "red", fontface = "bold") +
      labs(
        title = "Cumulative treatment reach rate (considering competitive risk)",
        subtitle = paste0("Total cases: ", nrow(Data_clean),", Treated cases: ", sum(Data_clean$EP_treat == 1), " ",
                          sprintf("(%.1f%%)", sum(Data_clean$EP_treat == 1)/nrow(Data_clean)*100),
                          ", Death: ", sum(Data_clean$censor == 1), " ",
                          sprintf("(%.1f%%)", sum(Data_clean$censor == 1)/nrow(Data_clean)*100),
                          "\nMedian observation: ", sprintf("%.1f", median(Data_clean$time_enroll_final)/365.25*12), " months"),
        x = "Months after CGP test",
        y = "Cumulative treatment reach rate (%)",
        caption = "Cumulative Incidence Function"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(0, max(6, cif_treat$time, na.rm = TRUE), 6)) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))
    ggiraph::girafe(ggobj = p2, width_svg = 8, height_svg = 8)
  }
})

output$figure_survival_CGP_5_1 = render_gt({
  req(input$figure_survival_CGP_5_1_var,
      figure_survival_CGP_5_1_trigger,
      figure_survival_CGP_5_1_trigger_2,
      OUTPUT_DATA$figure_surv_CGP_gene_names,
      OUTPUT_DATA$figure_surv_CGP_Factor_names,
      OUTPUT_DATA$figure_surv_CGP_Factor_names_univariant,
      OUTPUT_DATA$figure_surv_CGP_Data_forest_tmp)
  Factor_names = OUTPUT_DATA$figure_surv_CGP_Factor_names
  Factor_names_univariant = OUTPUT_DATA$figure_surv_CGP_Factor_names_univariant
  Data_forest_tmp = OUTPUT_DATA$figure_surv_CGP_Data_forest_tmp
  gene_names = OUTPUT_DATA$figure_surv_CGP_gene_names
  # 入力に応じてFactor_namesを調整
  if (input$figure_survival_CGP_5_1_var == "gene") {
    Factor_names <- Factor_names[Factor_names != "Cluster"]
    Factor_names_univariant <- Factor_names_univariant[Factor_names_univariant != "Cluster"]
  } else if (input$figure_survival_CGP_5_1_var == "cluster") {
    Factor_names <- Factor_names[!Factor_names %in% gene_names]
    Factor_names_univariant <- Factor_names_univariant[!Factor_names_univariant %in% gene_names]
  } else if (input$figure_survival_CGP_5_1_var == "None") {
    Factor_names <- Factor_names[!Factor_names %in% c(gene_names, "Cluster")]
    Factor_names_univariant <- Factor_names_univariant[!Factor_names_univariant %in% c(gene_names, "Cluster")]
  }
  max_samples = isolate(input$figure_survival_CGP_5_1_max_samples)
  max_samples = ifelse(is.null(max_samples), 5000, max_samples)
  if (nrow(Data_forest_tmp) > max_samples){
    sample_indices = sample(nrow(Data_forest_tmp),max_samples)
    Data_forest_tmp = Data_forest_tmp[sample_indices, ]
  } else {
    Data_forest_tmp = Data_forest_tmp
  }
  # 調整した変数リストでテーブル生成関数を呼び出す
  if (length(Factor_names) > 0) {
    create_gt_table(Data_forest_tmp, Factor_names, Factor_names_univariant)
  }
})

output$figure_survival_CGP_3_rawdata = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive, "FELIS downloaded raw data in Survival and clinical information tab")
})

ID_select_surv_CGP = reactive({
  req(OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive,
      OUTPUT_DATA$figure_surv_CGP_Data_MAF_target,
      OUTPUT_DATA$figure_surv_CGP_Data_drug,
      input$RMST_CGP)

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive
  Data_MAF_target = OUTPUT_DATA$figure_surv_CGP_Data_MAF_target
  Data_drug = OUTPUT_DATA$figure_surv_CGP_Data_drug
  Data_survival_interactive$Panel = Data_survival_interactive$症例.検体情報.パネル.名称.
  # 転移部位のマッピング
  metastasis_mapping <- c(
    "Lymph_met" = "Lymph_met",
    "Brain_met" = "Brain_met",
    "Lung_met" = "Lung_met",
    "Bone_met" = "Bone_met",
    "Liver_met" = "Liver_met"
  )

  # 初期IDセット
  ID_1 <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)

  # 入力プレフィックス
  input_prefix <- "gene_survival_reach_1_"

  # 遺伝子フィルタ（P_1: 必須変異1）
  p1_input <- input[[paste0(input_prefix, "P_1")]]
  if(!all(is.null(p1_input))) {
    ID_1 <- intersect(ID_1, (Data_MAF_target %>%
                               dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
  }

  # 遺伝子フィルタ（P_2: 必須変異2）
  p2_input <- input[[paste0(input_prefix, "P_2")]]
  if(!all(is.null(p2_input))) {
    ID_1 <- intersect(ID_1, (Data_MAF_target %>%
                               dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
  }

  # 遺伝子除外（W: 除外変異）
  w_input <- input[[paste0(input_prefix, "W")]]
  if(!all(is.null(w_input))) {
    ID_1 <- setdiff(ID_1, (Data_MAF_target %>%
                             dplyr::filter(Hugo_Symbol %in% w_input))$Tumor_Sample_Barcode)
  }

  # 臨床データフィルタ
  clinical_filters <- list(
    L = "CTx_lines_before_CGP",
    R = "pre_CGP_best_RECIST",
    A = "YoungOld",
    S = "症例.基本情報.性別.名称.",
    H = "Cancers",
    C = "cluster",
    EP_option = "EP_option",
    Panel = "症例.検体情報.パネル.名称.",
    P = "症例.背景情報.ECOG.PS.名称.",
    Best_Evidence_Level = "Best_Evidence_Level",
    Year = "year"
  )

  for(filter_key in names(clinical_filters)) {
    filter_input <- input[[paste0(input_prefix, filter_key)]]
    if(!all(is.null(filter_input))) {
      column_name <- clinical_filters[[filter_key]]
      filter_expr <- paste0(column_name, " %in% filter_input")
      ID_1 <- intersect(ID_1, (Data_survival_interactive %>%
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
    ID_1 <- intersect(ID_1, met_ids)
  }

  # 薬剤フィルタ（D: CGP前の薬剤）
  d_input <- input[[paste0(input_prefix, "D")]]
  if(!all(is.null(d_input))) {
    ID_1 <- intersect(ID_1, (Data_drug %>%
                               dplyr::filter(TxCGP == "Pre" & Drug %in% d_input))$ID)
  }
  nd_input <- input[[paste0(input_prefix, "ND")]]
  if(!all(is.null(nd_input))) {
    ID_1 <- setdiff(ID_1, (Data_drug %>%
                               dplyr::filter(TxCGP == "Pre" & Drug %in% nd_input))$ID)
  }

  # フィルタされたデータを返す
  return(Data_survival_interactive %>%
           dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>%
           dplyr::mutate(
             time_enroll_treat = case_when(
               time_enroll_treat <= 0 & EP_treat == 1 ~ time_enroll_final/2,
               !is.na(time_enroll_treat) & EP_treat == 1 ~ time_enroll_treat,
               is.na(time_enroll_treat) & EP_treat == 1 ~ time_enroll_final/2,
               TRUE ~ time_enroll_final
             ),
             # イベントタイプの定義：0=打ち切り, 1=治療, 2=死亡
             event_type = case_when(
               EP_treat == 1 ~ 1,           # 治療あり
               censor == 1 ~ 2,             # 死亡
               TRUE ~ 0                     # 打ち切り
             ),
             # 時間変数の設定（治療ありの場合は治療までの時間、そうでなければ観察終了時間）
             time_to_event = ifelse(EP_treat == 1, time_enroll_treat, time_enroll_final)
           )
  )
})

output$figure_survival_CGP_1 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive,
      OUTPUT_DATA$figure_surv_CGP_Data_MAF_target,
      OUTPUT_DATA$figure_surv_CGP_Data_drug,
      input$RMST_CGP)

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive
  Data_MAF_target = OUTPUT_DATA$figure_surv_CGP_Data_MAF_target
  Data_drug = OUTPUT_DATA$figure_surv_CGP_Data_drug
  Data_survival_interactive$Panel = Data_survival_interactive$症例.検体情報.パネル.名称.
  Data_survival_interactive$PS = Data_survival_interactive$症例.背景情報.ECOG.PS.名称.
  Data_survival_interactive$Sex = Data_survival_interactive$症例.基本情報.性別.名称.
  # 転移部位のマッピングを定義
  metastasis_mapping <- c(
    "Lymph_met" = "Lymph_met",
    "Brain_met" = "Brain_met",
    "Lung_met" = "Lung_met",
    "Bone_met" = "Bone_met",
    "Liver_met" = "Liver_met"
  )

  # ID抽出関数を定義
  extract_cgp_group_ids <- function(group_num) {
    # 初期IDセット
    IDs <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)

    # 動的に入力名を構築
    input_prefix <- paste0("gene_survival_CGP_", group_num, "_")

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
      S = "Sex",
      EP_option = "EP_option",
      EP_treat = "treat_group_2",
      H = "Cancers",
      C = "cluster",
      P = "PS",
      Panel = "Panel",
      Best_Evidence_Level = "Best_Evidence_Level",
      Year = "year"
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

    # 薬剤フィルタ（D: CGP後の薬剤として使用した場合）
    d_input <- input[[paste0(input_prefix, "D")]]
    if(!all(is.null(d_input))) {
      IDs <- intersect(IDs, (Data_drug %>%
                               dplyr::filter(TxCGP == "Post" & Drug %in% d_input))$ID)
    }
    # 薬剤フィルタ（ND: CGP後の薬剤として使用しなかった場合）
    nd_input <- input[[paste0(input_prefix, "ND")]]
    if(!all(is.null(nd_input))) {
      IDs <- setdiff(IDs, (Data_drug %>%
                               dplyr::filter(TxCGP == "Post" & Drug %in% nd_input))$ID)
    }

    return(IDs)
  }

  # Group1とGroup2のIDを取得
  ID_1 <- extract_cgp_group_ids(1)
  ID_2 <- extract_cgp_group_ids(2)

  # データ作成
  Data_survival_1 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>%
    dplyr::mutate(Group = 1)

  Data_survival_2 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>%
    dplyr::mutate(Group = 2)

  Data_survival <- rbind(Data_survival_1, Data_survival_2)

  if(!all(is.null(input$propensity_survival_CGP))) {
    Data_survival$Group_0_1 =　Data_survival$Group - 1
    formula_PS = formula(
      paste0("Group_0_1 ~", paste(input$propensity_survival_CGP, collapse = "+"))
    )
    matching = NA
    if(!all(Data_survival$Group_0_1 == Data_survival$Group_0_1[1])){
      if("Cancers" %in% input$propensity_survival_CGP && input$propensity_survival_cancer_complete_match_CGP == "Yes"){
        matching = MatchIt::matchit(formula_PS,
                                    data = Data_survival,
                                    method = "nearest",
                                    exact = "Cancers",
                                    distance = "logit",
                                    std.caliper = TRUE,
                                    caliper = 0.2, replace = FALSE)
      } else {
        matching = MatchIt::matchit(formula_PS,
                                    data = Data_survival,
                                    method = "nearest",
                                    distance = "logit",
                                    std.caliper = TRUE,
                                    caliper = 0.2, replace = FALSE)
      }


    }
    if(!is.na(matching)[[1]]){
      bal <- cobalt::bal.tab(matching, un = TRUE)
      bal_df <- bal$Balance
      bal_df <- bal_df[rownames(bal_df)!="distance", ]
      max_smd <- max(abs(bal_df$Diff.Adj), na.rm=TRUE)
      subtitle_balance <- paste0("PSM balance: max |SMD| = ",
                                 format(max_smd, digits = 2, nsmall = 2))
      # -------------------------------------------------------
      # ★追加: PS分布のプロット (マッチング前後比較)
      # -------------------------------------------------------
      # cobalt::bal.plot を使うと、Unadjusted(前) と Adjusted(後) を自動で並べてくれます

      p_ps_dist <- cobalt::bal.plot(
        matching,
        var.name = "distance", # 傾向スコア(logit)の分布を指定
        which = "both",        # マッチング前(Unadjusted)と後(Adjusted)の両方を表示
        type = "histogram",    # ヒストグラム
        mirror = TRUE,         # 上下対象(Mirrored)にして見やすく
        colors = c("#F8766D", "#00BFC4") # ggplot風の色
      ) +
        ggplot2::labs(
          title = "Propensity Score Distribution (Before & After Matching)",
          caption = "Unadjusted = Original / Adjusted = Matched"
        ) +
        ggplot2::theme(legend.position = "top")

      # PDFとして保存
      ggsave(file.path(tempdir(), "PS_distribution.pdf"), p_ps_dist, width = 8, height = 6)

      Data_matched <- MatchIt::match.data(matching)
      p <- cobalt::love.plot(
        bal,
        stats = "mean.diffs",
        abs = TRUE,
        thresholds = c(m = .1),
        var.order = "unadjusted",
        drop.distance = TRUE
      )
      out_path <- file.path(tempdir(), "love_plot_PSM.pdf")
      ggsave(filename = out_path, plot = p, width = 7, height = min(nrow(bal$Balance)/9+2, 45))
      Data_matched$pair_id <- Data_matched$subclass
      # CGPサバイバル解析実行
      survival_compare_and_plot_match(
        data = Data_matched,
        time_var = "time_enroll_final",
        status_var = "censor",
        group_var = "Group",
        input_rmst_cgp = input$RMST_CGP,
        plot_title = paste0("PS-matched survival, ", subtitle_balance),
        pair_var = "pair_id",
        n_boot = 2000,
        seed = 123
      )
    } else {
      # CGPサバイバル解析実行
      survival_compare_and_plot(
        data = Data_survival,
        time_var = "time_enroll_final",
        status_var = "censor",
        group_var = "Group",
        input_rmst_cgp = input$RMST_CGP,
        plot_title = "Survival analisys based on variants and drugs"
      )
    }
  } else if(!all(is.null(input$IPW_survival_CGP))) {
    Data_survival$Group_0_1 = Data_survival$Group - 1
    formula_PS = formula(
      paste0("Group_0_1 ~", paste(input$IPW_survival_CGP, collapse = "+"))
    )
    Data_iptw <- Data_survival
    subtitle_balance <- NA

    # 2群が混在しているときだけ実行
    if (!all(Data_iptw$Group_0_1 == Data_iptw$Group_0_1[1])) {

      # -------------------------------------------------------
      # 1) PS推定（ロジスティック回帰）
      # -------------------------------------------------------
      fit_ps <- glm(formula_PS, data = Data_iptw, family = binomial())
      ps <- predict(fit_ps, type = "response")

      # PSクリッピング
      eps <- 1e-6
      ps <- pmin(pmax(ps, eps), 1 - eps)

      # -------------------------------------------------------
      # 2) IPTW（治療重み）の計算
      # -------------------------------------------------------
      trt <- Data_iptw$Group_0_1
      p_trt <- mean(trt == 1, na.rm = TRUE)

      if(input$IPW_survival_CGP_ATE_ATT == "ATE"){
        w_trt <- ifelse(trt == 1, p_trt / ps, (1 - p_trt) / (1 - ps)) # ATE
      } else {
        w_trt <- ifelse(trt == 1, 1, ps / (1 - ps)) # ATT
      }

      # -------------------------------------------------------
      # 3) 外れ値（巨大重み）除外：治療重みに基づいてトリミング
      # -------------------------------------------------------
      thr <- suppressWarnings(as.numeric(input$IPW_threshold))
      if (is.na(thr) || thr <= 0) thr <- Inf

      keep <- is.finite(w_trt) & (w_trt <= thr)

      # データをフィルタリング（行が減る）
      Data_iptw <- Data_iptw[keep, , drop = FALSE]
      w_trt <- w_trt[keep]
      Data_iptw$.ps <- ps[keep]

      # -------------------------------------------------------
      # ★追加・修正: IPCWを実行するかどうかの分岐
      # -------------------------------------------------------

      # デフォルトは「打ち切り重み = 1」（補正なし）
      w_censor <- rep(1, nrow(Data_iptw))
      ipcw_status_text <- "" # タイトルに表示する用

      # UIで "Yes" が選択されている場合のみ計算を実行
      if (!is.null(input$IPCW_survival_CGP) && input$IPCW_survival_CGP == "Yes") {

        # A. 打ち切りステータスの作成 (1=打ち切り, 0=イベント発生or観察終了)
        Data_iptw$status_for_ipc <- 1 - Data_iptw$censor

        # B. 打ち切りモデルの構築 (Cox比例ハザード)
        # ※共変量はIPWと同じものを使用
        formula_censor <- formula(
          paste0("Surv(time_enroll_final, status_for_ipc) ~ Group_0_1 +",
                 paste(input$IPW_survival_CGP, collapse = "+"))
        )

        # 計算エラー回避（万が一モデルが収束しない場合など）
        tryCatch({
          fit_censor <- survival::coxph(formula_censor, data = Data_iptw)

          # C. 「打ち切られない確率」の推定 (累積ハザードから換算)
          chaz_censor <- predict(fit_censor, type = "expected")
          prob_uncensored <- exp(-chaz_censor)

          # D. 逆確率重み
          w_censor_calc <- 1 / prob_uncensored

          # クリッピング（上限10）
          w_censor <- pmin(w_censor_calc, as.numeric(input$IPW_threshold))

          # -------------------------------------------------------
          # ★追加: IPCWの適切性評価 (Diagnostics)
          # -------------------------------------------------------

          # (1) 打ち切りバランステスト (Censoring Balance)
          # 目的: 「打ち切られた群」と「残った群」の背景差が重みで補正されたか確認
          # treat引数に「打ち切りステータス(status_for_ipc)」を渡します

          covs_censor <- Data_iptw %>% dplyr::select(all_of(input$IPW_survival_CGP))

          bal_censor <- cobalt::bal.tab(
            x = covs_censor,
            treat = Data_iptw$status_for_ipc, # 1=Censored, 0=Event/End
            weights = w_censor,               # IPCW重み
            method = "weighting",
            un = TRUE
          )

          # Love Plot (Censoring)
          p_love_censor <- cobalt::love.plot(
            bal_censor,
            stats = "mean.diffs",
            abs = TRUE,
            thresholds = c(m = .1),
            var.order = "unadjusted",
            drop.distance = TRUE,
            title = "Covariate Balance for Censoring (IPCW)"
          ) +
            ggplot2::labs(subtitle = "Comparing 'Censored' vs 'Not Censored'")
          out_path <- file.path(tempdir(), "check_IPCW.pdf")
          ggsave(
            filename = out_path,
            plot = p_hist_censor,
            width = 7,
            height = min(nrow(bal_censor$Balance) / 9 + 2, 40)
          )


          # (2) 重みの分布確認 (Weight Histogram)
          # 目的: 極端な重みがないか、平均が1付近かを確認
          df_weights <- data.frame(w = w_censor)
          p_hist_censor <- ggplot2::ggplot(df_weights, aes(x = w)) +
            ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = "Distribution of Censoring Weights (IPCW)",
              subtitle = paste0("Mean: ", round(mean(w_censor), 3),
                                ", Max: ", round(max(w_censor), 3),
                                " (Capped at 10)"),
              x = "Inverse Probability of Censoring Weight",
              y = "Count"
            )
          out_path <- file.path(tempdir(), "IPCW_weight.pdf")
          ggsave(
            filename = out_path,
            plot = p_hist_censor,
            width = 7,
            height = 6
          )

          ipcw_status_text <- "+ IPCW"

        }, error = function(e) {
          showNotification("IPCW calculation failed. Using standard IPW only.", type = "warning")
        })
      }

      # -------------------------------------------------------
      # 4) 最終的な重みの合成
      # -------------------------------------------------------
      # Noの場合: w_trt * 1 = w_trt
      # Yesの場合: w_trt * w_censor
      Data_iptw$iptw_final <- w_trt * w_censor

      # -------------------------------------------------------
      # ★追加: 重みの正規化 (Normalization)
      # -------------------------------------------------------
      # 目的: No at risk が仮想人数で膨れ上がるのを防ぎ、実人数レベルに戻す

      sum_w <- sum(Data_iptw$iptw_final)
      n_obs <- nrow(Data_iptw)

      # 正規化された重み = 元の重み * (実人数 / 重みの総和)
      Data_iptw$iptw_final_norm <- Data_iptw$iptw_final * (n_obs / sum_w)

      # バランス評価用（常に治療重みのみ）
      Data_iptw$w_trt_only <- w_trt

      # -------------------------------------------------------
      # ★追加 1: 傾向スコアの分布プロット (Mirrored Histogram)
      # -------------------------------------------------------
      # 目的: 2群の重なり (Overlap) を視覚的に確認する
      # ps: 傾向スコア, Group_0_1: 0 or 1
      # -------------------------------------------------------
      # -------------------------------------------------------
      # ★追加 2: E-value の算出
      # -------------------------------------------------------

      # 関数: CoxモデルからE-valueを計算して表示用テキストを返す
      calc_evalue_from_cox <- function(model, label) {
        # 1. モデルから要約統計量を取得
        summ <- summary(model)

        # Group変数の行を特定（行名に "Group" を含むものを探す）
        target_rows <- grep("Group", rownames(summ$conf.int), value = TRUE)
        if (length(target_rows) == 0) {
          target_row_name <- rownames(summ$conf.int)[1] # 見つからなければ1行目
        } else {
          target_row_name <- target_rows[1] # 見つかればその先頭
        }

        # 2. 値の抽出（属性を完全に削除して純粋な数値にする）
        # as.numeric() で属性を剥がし、[1] で確実に1つの値にします
        hr_val      <- as.numeric(summ$conf.int[target_row_name, "exp(coef)"])[1]
        ci_low_val  <- as.numeric(summ$conf.int[target_row_name, "lower .95"])[1]
        ci_high_val <- as.numeric(summ$conf.int[target_row_name, "upper .95"])[1]

        # 3. E-value計算用の内部関数
        # 公式: E = HR + sqrt(HR * (HR - 1))
        # ※ HR < 1 の場合は、逆数 (1/HR) をとってから計算するのがルールです
        calc_formula <- function(val) {
          if (is.na(val) || !is.finite(val)) return(NA)
          if (val < 1) val <- 1 / val
          return(val + sqrt(val * (val - 1)))
        }

        # 4. 計算実行
        # 点推定値のE-value
        e_point <- round(calc_formula(hr_val), 2)

        # 信頼区間限界のE-value (より1に近い側の限界値を使います)
        # HR < 1 (予防的) なら、信頼区間の上限 (ci_high) が 1 に近い側
        # HR > 1 (有害)   なら、信頼区間の下限 (ci_low)  が 1 に近い側
        if (hr_val < 1) {
          e_ci <- round(calc_formula(ci_high_val), 2)
          # 信頼区間が1をまたいでいる（有意でない）場合、E-valueは定義上 1 になります
          if (ci_high_val >= 1) e_ci <- 1
        } else {
          e_ci <- round(calc_formula(ci_low_val), 2)
          if (ci_low_val <= 1) e_ci <- 1
        }

        # 5. 結果をリストで返す
        list(
          text = paste0(label, " HR: ", round(hr_val, 2),
                        " (", round(ci_low_val, 2), "-", round(ci_high_val, 2), ")\n",
                        "E-value: Point = ", e_point, ", CI Limit = ", e_ci),
          # 以前のコードとの互換性のため val も返しますが、中身は単純な行列にします
          val = matrix(c(e_point, e_ci), ncol = 2, dimnames = list(NULL, c("point", "lower")))
        )
      }

      # IPWのみ (治療重みだけ)
      fit_cox_ipw <- coxph(Surv(time_enroll_final, censor) ~ Group_0_1,
                           data = Data_iptw, weights = w_trt, robust = TRUE)
      ev_ipw <- calc_evalue_from_cox(fit_cox_ipw, "[Cox model]")

      # プロット用データ作成
      df_plot <- data.frame(
        PS = Data_iptw$.ps,
        Group = factor(Data_iptw$Group_0_1, labels = c("Group1", "Group2")),
        Weight = Data_iptw$iptw_final_norm # 正規化済み重みを使用（分布確認用）
      )

      p_hist_ps <- ggplot(df_plot, aes(x = PS, fill = Group)) +
        # 上向き (Group2)
        geom_histogram(data = subset(df_plot, Group == "Group2"),
                       aes(y = ..count..),
                       binwidth = 0.05, color = "black", alpha = 0.7) +
        # 下向き (Group1) - yをマイナスにする
        geom_histogram(data = subset(df_plot, Group == "Group1"),
                       aes(y = -..count..),
                       binwidth = 0.05, color = "black", alpha = 0.7) +
        scale_fill_manual(values = c("Group1" = "#F8766D", "Group2" = "#00BFC4")) +
        geom_hline(yintercept = 0, color = "white") +
        labs(
          title = "Distribution of Propensity Scores (Mirrored Histogram)",
          subtitle = paste0("Upper: Group2 / Lower: Group1\n", ev_ipw$text),
          x = "Propensity Score",
          y = "Count"
        ) +
        theme_minimal() +
        theme(legend.position = "top")

      # 保存 (UIで表示する場合は renderPlot 等に渡す)
      ggsave(file.path(tempdir(), "PS_distribution.pdf"), p_hist_ps, width = 7, height = 5)


      # -------------------------------------------------------
      # 5) バランス評価（cobalt）
      # ※ここにはIPCWの影響を入れないのが一般的（Baseline Balanceのため）
      # -------------------------------------------------------
      covs <- Data_iptw %>% dplyr::select(all_of(input$IPW_survival_CGP))

      bal <- cobalt::bal.tab(
        x = covs,
        treat = Data_iptw$Group_0_1,
        weights = Data_iptw$w_trt_only, # 治療重みのみ
        method = "weighting",
        estimand = "ATE",
        un = TRUE
      )

      bal_df <- bal$Balance
      max_smd <- max(abs(bal_df$Diff.Adj), na.rm = TRUE)
      subtitle_balance <- paste0(
        "IPTW balance: max |SMD| = ",
        format(max_smd, digits = 2, nsmall = 2),
        ", threshold(w) <= ", thr
      )

      p_love <- cobalt::love.plot(
        bal,
        stats = "mean.diffs",
        abs = TRUE,
        thresholds = c(m = .1),
        var.order = "unadjusted",
        drop.distance = TRUE,
        title = "Covariate Balance (Love Plot)"
      )
      out_path <- file.path(tempdir(), "love_plot_PSM.pdf")
      ggsave(
        filename = out_path,
        plot = p,
        width = 7,
        height = min(nrow(bal$Balance) / 9 + 2, 40)
      )

      # -------------------------------------------------------
      # ★追加: IPTW重みの分布プロット (Histogram)
      # -------------------------------------------------------

      # プロット用データフレーム作成
      # w_trt はトリミング後のデータに対応している必要があります
      df_w_trt <- data.frame(
        Weight = w_trt,
        Group = factor(Data_iptw$Group_0_1, labels = c("Group1", "Group2"))
      )

      # 統計量（最大値など）を字幕に入れると便利です
      summary_txt <- paste0(
        "Max Weight: ", round(max(w_trt), 2),
        " / Mean Weight: ", round(mean(w_trt), 2),
        "\n(Threshold applied: ", thr, ")"
      )

      p_hist_iptw <- ggplot2::ggplot(df_w_trt, ggplot2::aes(x = Weight, fill = Group)) +
        ggplot2::geom_histogram(bins = 30, position = "identity", alpha = 0.6, color = "white") +
        ggplot2::scale_fill_manual(values = c("Group1" = "#F8766D", "Group2" = "#00BFC4")) +
        ggplot2::labs(
          title = "Distribution of IPTW Weights",
          subtitle = summary_txt,
          x = "Weight",
          y = "Count"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top")

      # -------------------------------------------------------
      # 保存処理
      # -------------------------------------------------------
      out_path <- file.path(tempdir(), "IPW_weight.pdf")
      ggsave(
        filename = out_path,
        plot = p,
        width = 7,
        height = 6
      )


      # -------------------------------------------------------
      # 6) 生存解析（合成重みを使用）
      # -------------------------------------------------------
      # タイトル作成
      final_plot_title <- paste0("IPW ", ipcw_status_text, " weighted survival, ", subtitle_balance)

      survival_compare_and_plot_match(
        data = Data_iptw,
        time_var = "time_enroll_final",
        status_var = "censor",
        group_var = "Group",
        input_rmst_cgp = input$RMST_CGP,
        plot_title = final_plot_title,
        n_boot = 2000,
        seed = 123,
        weight_var = "iptw_final_norm" # ★ここを正規化版に変更
      )
    } else {
      # CGPサバイバル解析実行
      survival_compare_and_plot(
        data = Data_survival,
        time_var = "time_enroll_final",
        status_var = "censor",
        group_var = "Group",
        input_rmst_cgp = input$RMST_CGP,
        plot_title = "Survival analisys based on variants and drugs"
      )
    }
  } else {
    # CGPサバイバル解析実行
    survival_compare_and_plot(
      data = Data_survival,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "Group",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Survival analisys based on variants and drugs"
    )
  }
})

output$dl_love_plot_PSM <- downloadHandler(
  filename = function() {
    "love_plot_PSM.pdf"
  },
  content = function(file) {
    # Check the existence of the source file in temp directory
    src <- file.path(tempdir(), "love_plot_PSM.pdf")

    if (!file.exists(src)) {
      # Show a message to the user via notification instead of validate()
      showNotification("PDF is not available yet. Please generate the plot first.", type = "error")
      # Stop the download process quietly
      return(NULL)
    }

    # Copy file to the download target
    file.copy(src, file, overwrite = TRUE)
  }
)
output$dl_weight_IPW <- downloadHandler(
  filename = function() {
    "IPW_weight.pdf"
  },
  content = function(file) {
    src <- file.path(tempdir(), "IPW_weight.pdf")

    # ファイルが存在するかチェック (IPCWが実行されていない場合は警告)
    if (!file.exists(src)) {
      showNotification("IPW weight distribution not available. Please analyse IPW first.", type = "error")
      return(NULL)
    }

    file.copy(src, file, overwrite = TRUE)
  }
)
output$dl_weight_IPCW <- downloadHandler(
  filename = function() {
    "IPCW_weight.pdf"
  },
  content = function(file) {
    src <- file.path(tempdir(), "IPCW_weight.pdf")

    # ファイルが存在するかチェック (IPCWが実行されていない場合は警告)
    if (!file.exists(src)) {
      showNotification("IPCW weight distribution not available. Please analyse IPCW first.", type = "error")
      return(NULL)
    }

    file.copy(src, file, overwrite = TRUE)
  }
)
output$dl_check_IPCW <- downloadHandler(
  filename = function() {
    "IPCW_diagnostics.pdf"
  },
  content = function(file) {
    src <- file.path(tempdir(), "check_IPCW.pdf")

    # ファイルが存在するかチェック (IPCWが実行されていない場合は警告)
    if (!file.exists(src)) {
      showNotification("IPCW diagnostics not available. Please enable IPCW first.", type = "error")
      return(NULL)
    }

    file.copy(src, file, overwrite = TRUE)
  }
)
output$dl_PS_distribution <- downloadHandler(
  filename = function() {
    "PS_distribution.pdf"
  },
  content = function(file) {
    src <- file.path(tempdir(), "PS_distribution.pdf")

    # ファイルが存在するかチェック (IPCWが実行されていない場合は警告)
    if (!file.exists(src)) {
      showNotification("PS distribution not available. Please enable IPCW first.", type = "error")
      return(NULL)
    }

    file.copy(src, file, overwrite = TRUE)
  }
)

output$figure_survival_CGP_4 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CGP_oncogenic_genes,
      OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive,
      OUTPUT_DATA$figure_surv_CGP_Data_MAF_target)

  oncogenic_genes = OUTPUT_DATA$figure_surv_CGP_oncogenic_genes
  Data_survival = OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive
  Data_MAF_target = OUTPUT_DATA$figure_surv_CGP_Data_MAF_target
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
      traditional_fit = survfit(Surv(event = censor,
                                     time = time_enroll_final) ~ gene_mut,
                                data = Data_survival)
      if(nrow(Data_survival) != traditional_fit[[1]][[1]]){
        tau0 = ((max((Data_survival %>% dplyr::filter(gene_mut == 0))$time_enroll_final) - 1)/365.25)
        tau1 = ((max((Data_survival %>% dplyr::filter(gene_mut == 1))$time_enroll_final) - 1)/365.25)
        tau = floor(min(tau0, tau1, input$RMST_CGP) * 10) / 10
        verify <- survRM2::rmst2(
          time = Data_survival$time_enroll_final,
          status = Data_survival$censor,
          arm = Data_survival$gene_mut,
          tau = tau * 365.25,
          alpha = 0.05
        )

        diff_0 = survdiff(Surv(time_enroll_final, censor)~gene_mut,
                          data=Data_survival, rho=0)
        diff_1 = survdiff(Surv(time_enroll_final, censor)~gene_mut,
                          data=Data_survival, rho=1)
        diff_2 = coxph(Surv(time = time_enroll_final, censor)~gene_mut,
                       data=Data_survival)
        mut_gene = paste(oncogenic_genes$gene_mutation[i],
                         ", top ", i, " gene", sep="")
        tmp = data.frame(summary(traditional_fit)$table)
        gene_table$positive_median[i] = round(tmp$median[2] / 365.25 * 12, digits = 1)
        gene_table$positive_LL[i] = round(tmp$X0.95LCL[2] / 365.25 * 12, digits = 1)
        gene_table$positive_UL[i] = round(tmp$X0.95UCL[2] / 365.25 * 12, digits = 1)
        gene_table$negative_median[i] = round(tmp$median[1] / 365.25 * 12, digits = 1)
        gene_table$negative_LL[i] = round(tmp$X0.95LCL[1] / 365.25 * 12, digits = 1)
        gene_table$negative_UL[i] = round(tmp$X0.95UCL[1] / 365.25 * 12, digits = 1)
        gene_table$diff_median[i] = round(verify$unadjusted.result[1], digits = 2)
        gene_table$diff_LL[i] = round(verify$unadjusted.result[4], digits = 2)
        gene_table$diff_UL[i] = round(verify$unadjusted.result[7], digits = 2)
        hs[[k]] = surv_curv_entry(traditional_fit, Data_survival,
                                  paste0(tau, "-year RMST diff.: ",
                                         gene_table$diff_median[i],
                                         " (", gene_table$diff_LL[i], "-", gene_table$diff_UL[i], ") days"),
                                  c("Not mut", paste0(oncogenic_genes$gene_mutation[i], " mut")), diff_0, diff_1, diff_2)
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
  gene_table = gene_table %>% dplyr::mutate(
    Gene = factor(gene_table$Gene, levels = Gene_arrange))
  p_mid <-
    gene_table |>
    ggplot(aes(y = fct_rev(Gene))) +
    theme_classic() +
    geom_point(aes(x=diff_median), shape=15, size=3) +
    geom_linerange(aes(xmin=diff_LL, xmax=diff_UL)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    labs(x="Median OS (months) and difference in restricted mean survival time in 2 years (days)", y="Genes with oncogenic alteration") +
    coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                    xlim=c( - max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5,
                            max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)) +
    annotate("text",
             x = (-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5)/2,
             y = length(Gene_arrange) + 1,
             label = "mut=short survival") +
    annotate("text",
             x = (max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)/2,
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
        estimate_lab_3 = "RMST difference",
        patients = "Patients"
      ), gene_table)
  Gene_arrange = gene_table$Gene
  gene_table$Gene = factor(gene_table$Gene, levels = Gene_arrange)
  OUTPUT_DATA$figure_surv_CGP_gene_table = gene_table
  withProgress(message = "Drawing survival curves...", {

    p_left <- gene_table %>% ggplot(aes(y = fct_rev(Gene)))
    p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                 fontface = ifelse(gene_table$Gene == "Gene", "bold", "bold.italic"))
    p_left <- p_left + geom_text(aes(x = 1.8, label = estimate_lab_1), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_1 == "Survival, mut (+)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 3.4, label = estimate_lab_2), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_2 == "Survival, mut (-)", "bold", "plain"))
    p_left <- p_left + geom_text(aes(x = 5.0, label = estimate_lab_3), hjust = 0,
                                 fontface = ifelse(gene_table$estimate_lab_3 == "RMST difference", "bold", "plain"))
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
    OUTPUT_DATA$figure_surv_CGP_forest_plot = p_left + p_mid + p_right + plot_layout(design = layout)
    arrange_ggsurvplots(hs,print=TRUE,ncol=4,nrow=8,surv.plot.height = 0.8,risk.table.height = 0.2)
  })
})

output$figure_survival_CGP_3 = renderPlot({
  req(OUTPUT_DATA$figure_surv_CGP_forest_plot)
  withProgress(message = "Drawing forest plot...", {
    OUTPUT_DATA$figure_surv_CGP_forest_plot
  })
})
output$figure_survival_CGP_3_data = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$figure_surv_CGP_gene_table)
  create_datatable_with_confirm(OUTPUT_DATA$figure_surv_CGP_gene_table, "FELIS downloaded forest plot data in Survival and mutations, forest plot tab")
})

output$figure_surv_CGP = renderPlot({
  req(input$color_var_surv_CGP, input$RMST_CGP,
      OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive
      )
  Data_survival_interactive = OUTPUT_DATA$figure_surv_CGP_Data_survival_interactive
  Data_survival_interactive$Panel = Data_survival_interactive$症例.検体情報.パネル.名称.
  if(input$color_var_surv_CGP == "entire"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "1",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "All enrolled patients"
    )
  } else if(input$color_var_surv_CGP == "treat_option"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "EP_option",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Expert panel recommended treatment option existed or not"
    )
  } else if(input$color_var_surv_CGP == "treat_option_pos"){
    survival_compare_and_plot(
      data = Data_survival_interactive %>%
        dplyr::filter(EP_option == 1),
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "treat_group_2",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Only patients with expert panel recommended treatment option"
    )
  } else if(input$color_var_surv_CGP == "treat_group"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "treat_group",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "EP option and treatment"
    )
  } else if(input$color_var_surv_CGP == "treat_group_2"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "treat_group_2",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "EP option and treatment"
    )
  } else if(input$color_var_surv_CGP == "treat_group_3"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "treat_group_3",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "EP option and treatment"
    )
  } else if(input$color_var_surv_CGP == "treat_group_4"){
    survival_compare_and_plot(
      data = Data_survival_interactive %>%
        dplyr::mutate(CTx_lines_before_CGP = case_when(
          CTx_lines_before_CGP %in% c("0", "1", "2", "3") ~ CTx_lines_before_CGP,
          TRUE ~ "4~"
        )),
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "CTx_lines_before_CGP",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "CTx lines before CGP"
    )
  } else if(input$color_var_surv_CGP == "treat_group_5"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "pre_CGP_best_RECIST",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Best CTx effect before CGP"
    )
  } else if(input$color_var_surv_CGP == "treat_group_6"){
    survival_compare_and_plot(
      data = Data_survival_interactive %>%
        dplyr::filter(treat_group_2 != "No treatment done"),
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "treat_group_2",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Recommended vs not-recommended treatment"
    )
  } else if(input$color_var_surv_CGP == "treat_group_7"){
    Data_survival_interactive$PS = Data_survival_interactive$症例.背景情報.ECOG.PS.名称.
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "PS",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "ECOG Performance status"
    )
  } else if(input$color_var_surv_CGP == "treat_group_8"){
    Data_survival_interactive = Data_survival_interactive %>%
      dplyr::filter(症例.背景情報.ECOG.PS.名称. != "Unknown") %>%
      dplyr::mutate(PS = case_when(
        症例.背景情報.ECOG.PS.名称. == "0" ~ "0",
        症例.背景情報.ECOG.PS.名称. == "1" ~ "1",
        TRUE ~ "2_4"
      ))
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "PS",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "ECOG Performance status (unknown patients excluded)"
    )
  } else if(input$color_var_surv_CGP == "treat_group_9"){
    Data_survival_tmp = Data_survival_interactive %>% dplyr::filter(
      Cancers %in% names(sort(table(Data_survival_interactive$Cancers), decreasing = T))[1:min(7,length(unique(Data_survival_interactive$Cancers)))])
    Data_survival_tmp2 = Data_survival_interactive
    Data_survival_tmp2$Cancers = " ALL"
    Data_survival_tmp2 = rbind(Data_survival_tmp, Data_survival_tmp2)
    survival_compare_and_plot(
      data = Data_survival_tmp2,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "Cancers",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Diagnosis"
    )
  } else if(input$color_var_surv_CGP == "treat_group_10"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "Panel",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Panel"
    )
  } else if(input$color_var_surv_CGP == "treat_group_11"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "Best_Evidence_Level",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Best Evidence Level"
    )
  } else if(input$color_var_surv_CGP == "treat_group_12"){
    survival_compare_and_plot(
      data = Data_survival_interactive,
      time_var = "time_enroll_final",
      status_var = "censor",
      group_var = "year",
      input_rmst_cgp = input$RMST_CGP,
      plot_title = "Test year"
    )
  }
})
