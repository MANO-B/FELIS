custering_analysis_logic <- function() {
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
                      症例.検体情報.腫瘍細胞含有割合,
                      症例.検体情報.検体採取部位.名称.,
                      final_observe,
                      censor,
                      time_enroll_final,
                      time_palliative_final,
                      time_palliative_enroll
        )
      incProgress(1 / 13)
      Data_case_target$Cancers = Data_case_target$症例.基本情報.がん種.OncoTree.
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
      Data_MAF_target = Data_MAF_target  %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        Hugo_Symbol,
                        .keep_all = TRUE)
      incProgress(1 / 13)
      Cancername = sort(unique(Data_case_target$Cancers))
      Total_pts = as.vector(table((Data_case_target %>%
                                     dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = T))$Cancers))
      Gene_list = unique(names(sort(table(Data_MAF_target$Hugo_Symbol),
                                    decreasing = T)))
      if(!is.null(Gene_list)){
        Data_mutation_cord = Data_cluster_ID() %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
          dplyr::select(-C.CAT調査結果.基本項目.ハッシュID, -cluster, -driver_mutations)
        req(Data_mutation_cord)
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Tumor_Sample_Barcode %in% unique(Data_cluster_ID()$C.CAT調査結果.基本項目.ハッシュID))
        Data_case_target = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% unique(Data_cluster_ID()$C.CAT調査結果.基本項目.ハッシュID))
        Data_case_target = left_join(Data_case_target,
                                     Data_cluster_ID() %>% dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster, driver_mutations),
                                     by = "C.CAT調査結果.基本項目.ハッシュID")
        Data_case_target$cluster[is.na(Data_case_target$cluster)] = max(Data_case_target$cluster, na.rm = T) + 1
        Data_report_tmp = Data_report() %>%
          dplyr::filter(
            !str_detect(Hugo_Symbol, ",") &
              Hugo_Symbol != "" &
              Variant_Classification != "expression"
          )
        Data_report_tmp = Data_report_tmp %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
        Data_disease = (Data_case_target %>%
                          dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID,
                                          症例.基本情報.がん種.OncoTree.))
        colnames(Data_disease) = c("Tumor_Sample_Barcode", "Cancers")
        Data_MAF_target = left_join(Data_MAF_target,
                                    Data_disease,
                                    by = "Tumor_Sample_Barcode")
        Data_MAF_target$Cancers = as.factor(Data_MAF_target$Cancers)
        Data_report_tmp = left_join(Data_report_tmp,
                                    Data_disease,
                                    by = "Tumor_Sample_Barcode")
        Data_report_tmp$Cancers = as.factor(Data_report_tmp$Cancers)
        Data_report_tmp = Data_report_tmp %>% dplyr::mutate(
          Evidence_level = case_when(
            Evidence_level == "F" & Hugo_Symbol == "TMB" ~ "A",
            amino.acid.change == "high" & Hugo_Symbol == "MSI" ~ "A",
            Evidence_level == "" ~ "G",
            TRUE ~ Evidence_level),
          Drug = case_when(
            Evidence_level == "A" &
              Hugo_Symbol == "TMB" &
              Drug == "" ~ "pembrolizumab",
            Evidence_level == "A" &
              Hugo_Symbol == "TMB" &
              is.na(Drug) ~ "pembrolizumab",
            Evidence_level == "A" &
              Hugo_Symbol == "MSI" &
              Drug == "" ~ "pembrolizumab",
            Evidence_level == "A" &
              Hugo_Symbol == "MSI" &
              is.na(Drug) ~ "pembrolizumab",
            TRUE ~ Drug),
          Resistance = case_when(
            Evidence_level %in% c("R2*","R1*","R3*","R2","R3","R","R1") ~ 1,
            TRUE ~ 0)
        ) %>%
          dplyr::distinct(Tumor_Sample_Barcode, Evidence_level, Drug, .keep_all = T)
        incProgress(1 / 13)



        Data_evidence_table_tmp = Data_report_tmp %>%
          dplyr::filter(Evidence_level %in% c("A","B","C")) %>%
          dplyr::select(Hugo_Symbol,
                        # Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
                        amino.acid.change, Evidence_level, Drug, Cancers)
        Data_evidence_table_tmp$tmp = paste(Data_evidence_table_tmp$Hugo_Symbol,"_",
                                            Data_evidence_table_tmp$amino.acid.change,"_",
                                            Data_evidence_table_tmp$Evidence_level,"_",
                                            Data_evidence_table_tmp$Drug,"_",
                                            Data_evidence_table_tmp$Cancers)
        Data_evidence_table = Data_evidence_table_tmp %>%
          dplyr::distinct()
        Data_evidence_table_tmp2 = data.frame(table(Data_evidence_table_tmp$tmp))
        colnames(Data_evidence_table_tmp2) = c("tmp", "Frequency")
        Data_evidence_table = Data_evidence_table %>%
          dplyr::left_join(Data_evidence_table_tmp2, by="tmp") %>%
          dplyr::arrange(desc(Frequency))
        OUTPUT_DATA$clustering_Data_evidence_table = Data_evidence_table

        Data_report_tmp = Data_report_tmp %>%
          dplyr::arrange(Evidence_level) %>%
          dplyr::filter(Resistance != 1) %>%
          dplyr::filter(!Evidence_level %in% c("F", "G")) %>%
          dplyr::distinct(Tumor_Sample_Barcode, .keep_all = T)
        Data_report_tmp$Level = as.factor(Data_report_tmp$Evidence_level)
        data_figure = Data_report_tmp %>% dplyr::select(Cancers,Level,Tumor_Sample_Barcode) %>% dplyr::filter(!is.na(Cancers)&!is.na(Level))
        # Levelごとの件数を Cancername に基づいて補正する関数
        get_level_counts <- function(data, level, cancer_names, total_pts) {
          counts <- table(data %>% filter(Level == level) %>% pull(Cancers))
          result <- setNames(rep(0, length(cancer_names)), cancer_names)
          result[names(counts)] <- counts
          fun_zero(result, total_pts)
        }

        # 各レベルの集計
        levels <- c("A", "B", "C", "D", "E")
        level_counts <- lapply(levels, function(lvl) {
          get_level_counts(data_figure, lvl, Cancername, Total_pts)
        })
        names(level_counts) <- paste0("Num_", levels)
        Num_A = level_counts$Num_A
        Num_B = level_counts$Num_B
        Num_C = level_counts$Num_C
        Num_D = level_counts$Num_D
        Num_E = level_counts$Num_E


        # OptionとTxの処理も関数化
        get_flag_counts <- function(data, flag_col, cancer_names, total_pts) {
          counts <- table(data %>% filter(.data[[flag_col]] == 1) %>% pull(Cancers))
          result <- setNames(rep(0, length(cancer_names)), cancer_names)
          result[names(counts)] <- counts
          fun_zero(result, total_pts)
        }

        Num_Option <- get_flag_counts(Data_case_target, "EP_option", Cancername, Total_pts)
        Num_Receive_Tx <- get_flag_counts(Data_case_target, "EP_treat", Cancername, Total_pts)
        Num_Receive_Tx_per_Option <- fun_zero(Num_Receive_Tx, Num_Option)

        Total_pts_short = rep(0, length(Cancername))
        Total_pts_long = rep(0, length(Cancername))
        Num_pre_CGP = rep(0, length(Cancername))
        Num_post_CGP = rep(0, length(Cancername))

        Total_pts_Bone_no = rep(0, length(Cancername))
        Total_pts_Bone = rep(0, length(Cancername))
        Total_pts_Brain_no = rep(0, length(Cancername))
        Total_pts_Brain = rep(0, length(Cancername))
        Total_pts_Lung_no = rep(0, length(Cancername))
        Total_pts_Lung = rep(0, length(Cancername))
        Total_pts_Liver_no = rep(0, length(Cancername))
        Total_pts_Liver = rep(0, length(Cancername))

        Disease_cluster = data.frame(matrix(rep(0, length(Cancername) * max(Data_case_target$cluster,na.rm = T)),
                                            ncol = max(Data_case_target$cluster,na.rm = T)))
        colnames(Disease_cluster) = seq(1,max(Data_case_target$cluster,na.rm = T))
        rownames(Disease_cluster) = Cancername

        Data_case_target$diagnosis = Data_case_target$症例.基本情報.がん種.OncoTree.
        Data_survival = Data_case_target

        for (i in seq_along(Cancername)) {
          Data_survival_tmp <- Data_survival %>%
            filter(Cancers == Cancername[i])

          if (nrow(Data_survival_tmp) > 0) {

            ## --- pre_CGP ---
            if (sum(Data_survival_tmp$time_palliative_enroll, na.rm = TRUE) > 0) {
              survival_simple <- survfit(
                Surv(time_palliative_enroll, rep(1, nrow(Data_survival_tmp))) ~ 1,
                data = Data_survival_tmp
              )
              Num_pre_CGP[i] <- median(survival_simple)[[1]]
            } else {
              Num_pre_CGP[i] <- 0
            }

            ## --- post_CGP ---
            if (sum(Data_survival_tmp$time_enroll_final, na.rm = TRUE) > 0) {
              survival_simple <- survfit(
                Surv(time_enroll_final, censor) ~ 1,
                data = Data_survival_tmp
              )
              Num_post_CGP[i] <- median(survival_simple)[[1]]
            } else {
              Num_post_CGP[i] <- 0
            }

            ## --- Bone, Brain, Lung, Liver counts ---
            # ここで一度だけカウント
            counts <- Data_survival_tmp %>%
              summarise(
                Bone_yes  = sum(Bone_met  == "Yes", na.rm = TRUE),
                Bone_no   = sum(Bone_met  == "No",  na.rm = TRUE),
                Brain_yes = sum(Brain_met == "Yes", na.rm = TRUE),
                Brain_no  = sum(Brain_met == "No",  na.rm = TRUE),
                Lung_yes  = sum(Lung_met  == "Yes", na.rm = TRUE),
                Lung_no   = sum(Lung_met  == "No",  na.rm = TRUE),
                Liver_yes = sum(Liver_met == "Yes", na.rm = TRUE),
                Liver_no  = sum(Liver_met == "No",  na.rm = TRUE)
              )

            Total_pts_Bone[i]     <- counts$Bone_yes
            Total_pts_Bone_no[i]  <- counts$Bone_no
            Total_pts_Brain[i]    <- counts$Brain_yes
            Total_pts_Brain_no[i] <- counts$Brain_no
            Total_pts_Lung[i]     <- counts$Lung_yes
            Total_pts_Lung_no[i]  <- counts$Lung_no
            Total_pts_Liver[i]    <- counts$Liver_yes
            Total_pts_Liver_no[i] <- counts$Liver_no
          }
        }
        incProgress(1 / 13)
        Data_survival$timing = "Pre"
        for(i in 1:length(Cancername)){
          Data_survival[Data_survival$Cancers == unique(Data_survival$Cancers)[i],]$timing = ifelse(
            Data_survival[Data_survival$Cancers == unique(Data_survival$Cancers)[i],]$time_palliative_enroll < Num_pre_CGP[i],
            "Pre", "Post"
          )
        }
        Diseases = Cancername
        OUTPUT_DATA$clustering_Diseases = Diseases
        Summary_Table <- Diseases %>%
          purrr::map_dfr(function(d) {
            df <- Data_case_target %>% filter(症例.基本情報.がん種.OncoTree. == d)

            # 腫瘍部位のカウント
            site_counts <- table(df$症例.検体情報.検体採取部位.名称.)
            sample_is_primary_tumor   <- site_counts["原発巣"]   %||% 0
            sample_is_metastatic_tumor <- site_counts["転移巣"]   %||% 0
            sample_is_unknown_tumor   <- sum(site_counts[c("不明", "")], na.rm = TRUE)

            # 腫瘍細胞含有割合
            purity <- df$症例.検体情報.腫瘍細胞含有割合
            tumor_purity_0_25   <- sum(!is.na(purity) & purity < 25)
            tumor_purity_25_50  <- sum(!is.na(purity) & purity >= 25 & purity < 50)
            tumor_purity_50_75  <- sum(!is.na(purity) & purity >= 50 & purity < 75)
            tumor_purity_75_100 <- sum(!is.na(purity) & purity >= 75 & purity <= 100)
            tumor_purity_NA     <- sum(is.na(purity))

            # 年齢統計
            age <- df$症例.基本情報.年齢
            age_quant <- quantile(age, na.rm = TRUE)
            age_median <- age_quant[3]
            age_lowest <- age_quant[1]
            age_highest <- age_quant[5]
            age_range_25 <- age_quant[2]
            age_range_75 <- age_quant[4]

            # 性別
            male         <- sum(df$症例.基本情報.性別.名称. == "Male", na.rm = TRUE)
            female       <- sum(df$症例.基本情報.性別.名称. == "Female", na.rm = TRUE)
            sex_unknown  <- sum(df$症例.基本情報.性別.名称. == "Unknown", na.rm = TRUE)

            # driver mutations
            oncogenic_driver_plus  <- sum(df$driver_mutations > 0, na.rm = TRUE)
            oncogenic_driver_minus <- sum(df$driver_mutations == 0, na.rm = TRUE)

            # TMB統計
            tmb <- df$TMB
            tmb_quant <- quantile(tmb, na.rm = TRUE)
            TMB_median <- tmb_quant[3]
            TMB_lowest <- tmb_quant[1]
            TMB_highest <- tmb_quant[5]
            TMB_range_25 <- tmb_quant[2]
            TMB_range_75 <- tmb_quant[4]

            # cluster集計
            tmp_cluster <- table(df$cluster)
            entropy_val <- shannon.entropy(as.numeric(tmp_cluster), max(Data_case_target$cluster,na.rm = T))

            # Disease_cluster更新（別途処理が必要）
            for (j in seq_along(tmp_cluster)) {
              Disease_cluster[d, as.numeric(names(tmp_cluster)[j])] <<- tmp_cluster[j]
            }

            tibble(
              Diseases = d,
              sample_is_primary_tumor   = sample_is_primary_tumor,
              sample_is_metastatic_tumor = sample_is_metastatic_tumor,
              sample_is_unknown_tumor   = sample_is_unknown_tumor,
              tumor_purity_0_25   = tumor_purity_0_25,
              tumor_purity_25_50  = tumor_purity_25_50,
              tumor_purity_50_75  = tumor_purity_50_75,
              tumor_purity_75_100 = tumor_purity_75_100,
              tumor_purity_NA     = tumor_purity_NA,
              age_median = age_median,
              age_lowest = age_lowest,
              age_highest = age_highest,
              age_range_25 = age_range_25,
              age_range_75 = age_range_75,
              male = male,
              female = female,
              sex_unknown = sex_unknown,
              oncogenic_driver_plus = oncogenic_driver_plus,
              oncogenic_driver_minus = oncogenic_driver_minus,
              TMB_median = TMB_median,
              TMB_lowest = TMB_lowest,
              TMB_highest = TMB_highest,
              TMB_range_25 = TMB_range_25,
              TMB_range_75 = TMB_range_75,
              entropy = entropy_val
            )
          })
        Summary_Table$brain = Total_pts_Brain
        Summary_Table$brain_no = Total_pts_Brain_no
        Summary_Table$bone = Total_pts_Bone
        Summary_Table$bone_no = Total_pts_Bone_no
        Summary_Table$lung = Total_pts_Lung
        Summary_Table$lung_no = Total_pts_Lung_no
        Summary_Table$liver = Total_pts_Liver
        Summary_Table$liver_no = Total_pts_Liver_no

        Summary_Table$option = Num_Option * 100
        Summary_Table$treat = Num_Receive_Tx * 100
        Summary_Table$time_before_CGP = Num_pre_CGP
        Summary_Table$time_after_CGP = Num_post_CGP
        incProgress(1 / 13)

        testSummary = data.frame(as.vector(t(matrix(rep(1:max(Data_case_target$cluster,na.rm = T),
                                                        length(Diseases)),
                                                    ncol =length(Diseases)))))
        colnames(testSummary) = "Cluster"
        testSummary$Histology = rep(Diseases, max(Data_case_target$cluster,na.rm = T))
        testSummary$All_patients = rep(rep(0, length(Diseases)), max(Data_case_target$cluster,na.rm = T))
        testSummary$Positive_patients = rep(rep(0, length(Diseases)), max(Data_case_target$cluster,na.rm = T))
        testSummary$OddsRatio = rep(rep(1, length(Diseases)), max(Data_case_target$cluster,na.rm = T))
        testSummary$Pvalue = rep(rep(1, length(Diseases)), max(Data_case_target$cluster,na.rm = T))

        selected_genes = rep("", max(Data_survival$cluster,na.rm = T))
        cluster_set = sort(unique(Data_survival$cluster))
        for(i in 1:length(cluster_set)){
          pos_num = data.frame(Diseases)
          pos_num$n = 0
          colnames(pos_num) = c("diagnosis", "n")
          Data_cluster_plus = Data_survival %>%
            dplyr::filter(cluster == cluster_set[i]) %>%
            dplyr::count(diagnosis) %>%
            dplyr::arrange(-n)
          Data_cluster_plus = dplyr::full_join(Data_cluster_plus, pos_num, by = join_by(diagnosis))
          Data_cluster_plus = Data_cluster_plus[,1:2]
          colnames(Data_cluster_plus) = c("diagnosis", "n")
          Data_cluster_plus$n = replace_na(Data_cluster_plus$n, 0)
          Data_cluster_plus = Data_cluster_plus %>% dplyr::arrange(diagnosis)
          Data_cluster_plus = Data_cluster_plus$n
          pos_num = data.frame(Diseases)
          pos_num$n = 0
          colnames(pos_num) = c("diagnosis", "n")
          Data_cluster_minus = Data_survival %>%
            dplyr::filter(cluster != cluster_set[i]) %>%
            dplyr::count(diagnosis) %>%
            dplyr::arrange(-n)
          Data_cluster_minus = dplyr::full_join(Data_cluster_minus, pos_num, by = join_by(diagnosis))
          Data_cluster_minus = Data_cluster_minus[,1:2]
          colnames(Data_cluster_minus) = c("diagnosis", "n")
          Data_cluster_minus$n = replace_na(Data_cluster_minus$n, 0)
          Data_cluster_minus = Data_cluster_minus %>% dplyr::arrange(diagnosis)
          Data_cluster_minus = Data_cluster_minus$n
          No_plus_neg = length((Data_survival %>% dplyr::filter(cluster == cluster_set[i]))[,1][[1]]) -
            Data_cluster_plus
          No_minus_neg = length((Data_survival %>% dplyr::filter(cluster != cluster_set[i]))[,1][[1]]) -
            Data_cluster_minus
          testSummary[((1+(i-1)*length(Diseases)):(i*length(Diseases))),3] = length((Data_survival %>% dplyr::filter(cluster == cluster_set[i]))[,1][[1]])
          testSummary[((1+(i-1)*length(Diseases)):(i*length(Diseases))),4] = Data_cluster_plus
          for(j in 1:length(Data_cluster_plus)){
            mutation_matrix = matrix(c(Data_cluster_plus[j],
                                       Data_cluster_minus[j],
                                       No_plus_neg[j],
                                       No_minus_neg[j]),
                                     nrow = 2,
                                     dimnames = list(Cancer = c("CancerType", "OtherCancerType"),
                                                     Cluster = c(paste("cluster", cluster_set[i]), "Other clusters")))
            testSummary[(j+(i-1)*length(Diseases)),5] =
              sprintf("%.3f", fisher.test(mutation_matrix)$estimate[[1]])
            testSummary[(j+(i-1)*length(Diseases)),6] =
              sprintf("%.3f", fisher.test(mutation_matrix)$p.value)
          }
          top_genes = testSummary[((1+(i-1)*length(Diseases)):(i*length(Diseases))),] %>%
            dplyr::arrange(desc(OddsRatio)) %>%
            dplyr::filter(OddsRatio > 1 & Pvalue < 0.05)
          if(length(top_genes$OddsRatio) >= 1){
            top_genes = top_genes[1:min(3, length(top_genes$OddsRatio)),]$Histology
          } else{
            top_genes = ""
          }
          selected_genes[i] = paste0(cluster_set[i], ":", paste(top_genes,collapse = ";"))
        }
        testSummary[,"OddsRatio"] = ifelse(is.infinite(testSummary[,"OddsRatio"]), 99999, testSummary[,"OddsRatio"])
        OUTPUT_DATA$clustering_testSummary_disease = testSummary

        incProgress(1 / 13)

        Data_cluster = Data_case_target %>% dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster)
        colnames(Data_cluster) = c("Tumor_Sample_Barcode", "cluster")
        Data_MAF_target = left_join(Data_MAF_target, Data_cluster, by = "Tumor_Sample_Barcode")
        Mutations = sort(unique(Data_MAF_target$Hugo_Symbol))
        testSummary = data.frame(as.vector(t(matrix(rep(1:max(Data_MAF_target$cluster,na.rm = T),
                                                        length(Mutations)),
                                                    ncol =length(Mutations)))))
        colnames(testSummary) = "Cluster"
        testSummary$mutation = rep(Mutations, max(Data_MAF_target$cluster,na.rm = T))
        testSummary$All_patients = rep(rep(0, length(Mutations)), max(Data_MAF_target$cluster,na.rm = T))
        testSummary$Positive_patients = rep(rep(0, length(Mutations)), max(Data_MAF_target$cluster,na.rm = T))
        testSummary$OddsRatio = rep(rep(1, length(Mutations)), max(Data_MAF_target$cluster,na.rm = T))
        testSummary$Pvalue = rep(rep(1, length(Mutations)), max(Data_MAF_target$cluster,na.rm = T))

        cluster_set = sort(unique(Data_MAF_target$cluster))
        for(i in 1:length(cluster_set)){
          pos_num = data.frame(Mutations)
          pos_num$n = 0
          colnames(pos_num) = c("Hugo_Symbol", "n")
          Data_cluster_plus = Data_MAF_target %>%
            dplyr::filter(cluster == cluster_set[i]) %>%
            dplyr::count(Hugo_Symbol) %>%
            dplyr::arrange(-n)
          Data_cluster_plus = dplyr::full_join(Data_cluster_plus, pos_num, by = join_by(Hugo_Symbol))
          Data_cluster_plus = Data_cluster_plus[,1:2]
          colnames(Data_cluster_plus) = c("Hugo_Symbol", "n")
          Data_cluster_plus$n = replace_na(Data_cluster_plus$n, 0)
          Data_cluster_plus = Data_cluster_plus %>% dplyr::arrange(Hugo_Symbol)
          Data_cluster_plus = Data_cluster_plus$n
          pos_num = data.frame(Mutations)
          pos_num$n = 0
          colnames(pos_num) = c("Hugo_Symbol", "n")
          Data_cluster_minus = Data_MAF_target %>%
            dplyr::filter(cluster != cluster_set[i]) %>%
            dplyr::count(Hugo_Symbol) %>%
            dplyr::arrange(-n)
          Data_cluster_minus = dplyr::full_join(Data_cluster_minus, pos_num, by = join_by(Hugo_Symbol))
          Data_cluster_minus = Data_cluster_minus[,1:2]
          colnames(Data_cluster_minus) = c("Hugo_Symbol", "n")
          Data_cluster_minus$n = replace_na(Data_cluster_minus$n, 0)
          Data_cluster_minus = Data_cluster_minus %>% dplyr::arrange(Hugo_Symbol)
          Data_cluster_minus = Data_cluster_minus$n
          No_plus_neg = length((Data_MAF_target %>%
                                  dplyr::filter(cluster == cluster_set[i]) %>%
                                  dplyr::distinct(Tumor_Sample_Barcode))[,1][[1]]) -
            Data_cluster_plus
          No_minus_neg = length((Data_MAF_target %>%
                                   dplyr::filter(cluster != cluster_set[i]) %>%
                                   dplyr::distinct(Tumor_Sample_Barcode))[,1][[1]]) -
            Data_cluster_minus
          testSummary[((1+(i-1)*length(Mutations)):(i*length(Mutations))),3] =
            length((Data_MAF_target %>%
                      dplyr::filter(cluster == cluster_set[i]) %>%
                      dplyr::distinct(Tumor_Sample_Barcode))[,1][[1]])
          testSummary[((1+(i-1)*length(Mutations)):(i*length(Mutations))),4] = Data_cluster_plus
          for(j in 1:length(Data_cluster_plus)){
            mutation_matrix = matrix(c(Data_cluster_plus[j],
                                       Data_cluster_minus[j],
                                       No_plus_neg[j],
                                       No_minus_neg[j]),
                                     nrow = 2,
                                     dimnames = list(Cancer = c("CancerType", "OtherCancerType"),
                                                     Cluster = c(paste("cluster", cluster_set[i]), "Other clusters")))
            testSummary[(j+(i-1)*length(Mutations)),5] =
              sprintf("%.3f", fisher.test(mutation_matrix)$estimate[[1]])
            testSummary[(j+(i-1)*length(Mutations)),6] =
              sprintf("%.3f", fisher.test(mutation_matrix)$p.value)
          }
          top_genes = testSummary[((1+(i-1)*length(Mutations)):(i*length(Mutations))),] %>%
            dplyr::arrange(desc(OddsRatio)) %>%
            dplyr::filter(OddsRatio > 1 & Pvalue < 0.05)
          if(length(top_genes$OddsRatio) >= 1){
            top_genes = top_genes[1:min(3, length(top_genes$OddsRatio)),]$mutation
          } else{
            top_genes = ""
          }
          selected_genes[i] = paste0(selected_genes[i], "/", paste(top_genes,collapse = ";"))
        }
        incProgress(1 / 13)
        OUTPUT_DATA$clustering_selected_genes = selected_genes
        testSummary[,"OddsRatio"] = ifelse(is.infinite(testSummary[,"OddsRatio"]), 99999, testSummary[,"OddsRatio"])
        OUTPUT_DATA$clustering_testSummary_mutation = testSummary

        Disease_cluster$Disease = rownames(Disease_cluster)
        Disease_cluster <- transform(Disease_cluster, Disease= factor(Disease, levels = sort(unique(Disease_cluster$Disease), decreasing = TRUE)))
        colnames(Disease_cluster) = c(seq(1,max(Data_survival$cluster,na.rm = T)), "Disease")
        Data_entropy = gather(Disease_cluster, key = cluster, value = samples, -Disease)
        Data_entropy$cluster = as.numeric(Data_entropy$cluster)
        Summary_Table = transform(Summary_Table, Diseases= factor(Diseases, levels = sort(unique(Disease_cluster$Disease), decreasing = FALSE)))
        Data_entropy = transform(Data_entropy, Disease= factor(Disease, levels = sort(unique(Disease_cluster$Disease), decreasing = FALSE)))
        Data_entropy$Patients = Total_pts
        Data_entropy$percent <- Data_entropy$samples / Data_entropy$Patients
        Data_entropy$tooltip_text <- paste0(
          "Cancer type: ", Data_entropy$Disease, "\n",
          "Cluster: ", Data_entropy$cluster, "\n",
          "Total Patients: ", Data_entropy$Patients, "\n",
          "Cases: ", Data_entropy$samples, "\n",
          "Percent: ", percent(Data_entropy$percent, accuracy = 0.1)
        )
        OUTPUT_DATA$clustering_Summary_Table = Summary_Table
        OUTPUT_DATA$clustering_Data_entropy = Data_entropy

        x <- data.frame(
          Cancers   = Diseases,
          Patients = Total_pts,
          Age = Summary_Table$age_median,
          Level_A = (Num_A * 100),
          Level_AB = (Num_A + Num_B) * 100,
          Level_ABC = (Num_A + Num_B + Num_C) * 100,
          Treatment_option = Num_Option * 100,
          Receive_Recom_Tx = Num_Receive_Tx * 100,
          Receive_Tx_rate_with_opt = Num_Receive_Tx_per_Option * 100,
          Survival_pre_CGP = Num_pre_CGP,
          Survival_post_CGP = Num_post_CGP
        )
        x$Cancers = factor(x$Cancers)
        x_max_patients = max(x$Patients)

        y <- data.frame(
          Cancers   = rep(Diseases,5),
          Patients   = rep(Total_pts,5),
          Patients_with_treatment_option   = rep(round(Total_pts*Num_Option, digits = 0),5),
          Patients_received_recommended_treatment   = rep(round(Total_pts*Num_Receive_Tx, digits = 0),5),
          Level = c(rep("A", length(Diseases)),
                    rep("B", length(Diseases)),
                    rep("C", length(Diseases)),
                    rep("D", length(Diseases)),
                    rep("E", length(Diseases))),
          Patients_with_evidence_level = c(round(Total_pts*Num_A, digits = 0),
                                           round(Total_pts*Num_B, digits = 0),
                                           round(Total_pts*Num_C, digits = 0),
                                           round(Total_pts*Num_D, digits = 0),
                                           round(Total_pts*Num_E, digits = 0)),
          weight = c(Num_A, Num_B, Num_C, Num_D, Num_E)
        )
        y$Cancers = factor(y$Cancers)
        y$Level = factor(y$Level)
        y$percent <- y$Patients_with_evidence_level / y$Patients
        y$tooltip <- paste0(
          "Cancer type: ", y$Cancers, "\n",
          "Total Patients: ", y$Patients, "\n",
          "Evidence Level: ", y$Level, "\n",
          "Cases: ", y$Patients_with_evidence_level, "\n",
          "Percent: ", percent(y$percent, accuracy = 0.1)
        )
        incProgress(1 / 13)

        get_level_counts_cluster <- function(data, level, clusters, total_pts) {
          counts <- table(data %>% filter(Level == level) %>% pull(cluster))
          result <- setNames(rep(0, length(clusters)), clusters)
          result[names(counts)] <- counts
          fun_zero(result, total_pts)
        }
        data_figure_cluster = left_join(data_figure, Data_cluster, by = "Tumor_Sample_Barcode")
        All_cluster = sort(unique(Data_case_target$cluster))
        Total_pts_cluster = as.vector(table((Data_case_target %>%
                                               dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = T))$cluster))

        levels <- c("A", "B", "C", "D", "E")
        level_counts_cluster <- lapply(levels, function(lvl) {
          get_level_counts_cluster(data_figure_cluster, lvl, All_cluster, Total_pts_cluster)
        })
        names(level_counts_cluster) <- paste0("Num_", levels)
        z <- data.frame(
          Clusters   = rep(All_cluster,5),
          Patients   = rep(Total_pts_cluster,5),
          Level = c(rep("A", length(All_cluster)),
                    rep("B", length(All_cluster)),
                    rep("C", length(All_cluster)),
                    rep("D", length(All_cluster)),
                    rep("E", length(All_cluster))),
          Patients_with_evidence_level = c(round(Total_pts_cluster*level_counts_cluster$Num_A, digits = 0),
                                           round(Total_pts_cluster*level_counts_cluster$Num_B, digits = 0),
                                           round(Total_pts_cluster*level_counts_cluster$Num_C, digits = 0),
                                           round(Total_pts_cluster*level_counts_cluster$Num_D, digits = 0),
                                           round(Total_pts_cluster*level_counts_cluster$Num_E, digits = 0)),
          weight = c(level_counts_cluster$Num_A, level_counts_cluster$Num_B, level_counts_cluster$Num_C, level_counts_cluster$Num_D, level_counts_cluster$Num_E)
        )
        z$Clusters = factor(z$Clusters)
        z$Level = factor(z$Level)
        z$percent <- z$Patients_with_evidence_level / z$Patients
        z$tooltip <- paste0(
          "Cluster: ", z$Clusters, "\n",
          "Total Patients: ", z$Patients, "\n",
          "Evidence Level: ", z$Level, "\n",
          "Cases: ", z$Patients_with_evidence_level, "\n",
          "Percent: ", percent(z$percent, accuracy = 0.1)
        )
        OUTPUT_DATA$clustering_y = y
        OUTPUT_DATA$clustering_z = z


        Data_case_target$cluster = as.factor(Data_case_target$cluster)
        OUTPUT_DATA$clustering_Data_case_target = Data_case_target
        Data_cord_tmp = Data_case_target %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% Data_cluster_ID()$C.CAT調査結果.基本項目.ハッシュID)
        Data_mutation_cord$cluster = Data_cord_tmp$cluster
        Data_mutation_cord$Cancers = Data_cord_tmp$Cancers
        Data_mutation_cord$EP_option = as.factor(Data_cord_tmp$EP_option)
        Data_mutation_cord$EP_treat = as.factor(Data_cord_tmp$EP_treat)
        Data_mutation_cord$tooltip_text <- paste0(
          "Cluster: ", Data_mutation_cord$cluster, "\n",
          "Histology: ", Data_mutation_cord$Cancers, "\n",
          "Treatment option recommendation: ", Data_mutation_cord$EP_option, "\n",
          "Recommended treatment done: ", Data_mutation_cord$EP_treat
        )
        OUTPUT_DATA$clustering_Data_mutation_cord = Data_mutation_cord
      }
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}


output$table_basic_data = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$clustering_Summary_Table)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_Summary_Table, "FELIS downloaded forest plot data in Basic data tab")
})


output$table_disease = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$clustering_testSummary_disease)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_testSummary_disease, "FELIS downloaded disease data in UMAP clustering based on mutations tab")
})
output$table_mutation = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$clustering_testSummary_mutation)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_testSummary_mutation, "FELIS downloaded mutation cluster data in UMAP clustering based on mutations tab")
})
output$table_drug_evidence_y = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$clustering_y)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_y, "FELIS downloaded disease data in Frequency of patients with targeted therapy tab")
})
output$table_drug_evidence_z = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$clustering_z)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_z, "FELIS downloaded mutation data in Frequency of patients with targeted therapy tab")
})
output$figure_entropy = renderGirafe({
  req(OUTPUT_DATA$clustering_Summary_Table)
  req(OUTPUT_DATA$clustering_Data_entropy)
   g1 = ggplot(OUTPUT_DATA$clustering_Summary_Table, aes(x=Diseases, y=entropy, fill="black")) +
    geom_point(stat = "identity") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none")

  g3 = ggplot(OUTPUT_DATA$clustering_Data_entropy, aes(x=Disease, y=samples,tooltip = tooltip_text, fill=cluster)) +
    geom_bar_interactive(stat = "identity", position = "fill") +
    scale_fill_gradientn_interactive(colours = c("green", "black", "magenta", "darkred", "orange", "blue", "yellow")) +
    coord_flip() +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
  combined <- cowplot::plot_grid(g1, g3, rel_widths = c(1, 2), nrow = 1)
  height_svg = 1.2 + length(unique(OUTPUT_DATA$clustering_Data_entropy$Disease))*0.15

  girafe(ggobj = combined, width_svg = 12, height_svg = height_svg)
})

output$figure_base = renderPlot({
  req(OUTPUT_DATA$clustering_Data_case_target,
      OUTPUT_DATA$clustering_Summary_Table)
  withProgress(message = paste0("Drawing...", sample(nietzsche)[1]), {
    Data_case_target = OUTPUT_DATA$clustering_Data_case_target
    Data_tmp = Data_case_target
    Summary_Table = OUTPUT_DATA$clustering_Summary_Table
    Data_tmp$Age = Data_tmp$症例.基本情報.年齢
    Data_tmp <- transform(Data_tmp, diagnosis= factor(diagnosis, levels = sort(unique(Data_tmp$diagnosis), decreasing = TRUE)))

    chart_3 = ggplot(data = Data_tmp, aes(x = diagnosis, y = Age)) +
      geom_boxplot(outlier.colour = NA) +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank())
    Diseases = Summary_Table$Diseases
    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$sex = Summary_Table$male
    Data_tmp_1$Sex = "Male"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$sex = Summary_Table$female
    Data_tmp_2$Sex = "Female"
    Data_tmp_3 =data.frame(Summary_Table$Diseases)
    Data_tmp_3$sex = Summary_Table$sex_unknown
    Data_tmp_3$Sex = "Unknown"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2, Data_tmp_3)
    Total_pts_num = sum(Data_tmp$sex)
    Data_tmp$sex = fun_zero(Data_tmp$sex, Total_pts_num)
    cancer_freq_order = c("Unknown", "Female", "Male")
    Data_tmp <- transform(Data_tmp, Sex= factor(Sex, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Diseases), decreasing = TRUE)))

    chart_4 = ggplot(data = Data_tmp, aes(x = Diseases, y = sex, fill=Sex)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(#legend.position = "none",
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Sex))])) +
      ylab("Sex") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$driver = Summary_Table$oncogenic_driver_plus
    Data_tmp_1$Driver = "Yes"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$driver = Summary_Table$oncogenic_driver_minus
    Data_tmp_2$Driver = "No"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
    Total_pts_num = sum(Data_tmp$driver)
    Data_tmp$driver = fun_zero(Data_tmp$driver, Total_pts_num)
    cancer_freq_order = c("Unknown", "No", "Yes")
    Data_tmp <- transform(Data_tmp, Driver= factor(Driver, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Diseases), decreasing = TRUE)))

    chart_5 = ggplot(data = Data_tmp, aes(x = Diseases, y = driver, fill=Driver)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(#legend.position = "none",
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Driver))])) +
      ylab("Oncogenic muts detected") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp = Data_case_target
    Data_tmp <- transform(Data_tmp, diagnosis= factor(diagnosis, levels = sort(unique(Diseases), decreasing = TRUE)))
    chart_6 = ggplot(data = Data_tmp, aes(x = diagnosis, y = TMB)) +
      geom_boxplot(outlier.colour = NA) +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(limits = c(0, NA))

    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$brain_meta = Summary_Table$brain
    Data_tmp_1$Brain_meta = "(+)"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$brain_meta = Summary_Table$brain_no
    Data_tmp_2$Brain_meta = "(-)"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
    Total_pts_num = sum(Data_tmp$brain_meta)
    Data_tmp$brain_meta = fun_zero(Data_tmp$brain_meta, Total_pts_num)
    cancer_freq_order = c("(-)", "(+)")
    Data_tmp <- transform(Data_tmp, Brain_meta= factor(Brain_meta, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_7 = ggplot(data = Data_tmp, aes(x = Diseases, y = brain_meta, fill=Brain_meta)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Brain_meta))])) +
      ylab("Brain metastasis") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$bone_meta = Summary_Table$bone
    Data_tmp_1$Bone_meta = "(+)"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$bone_meta = Summary_Table$bone_no
    Data_tmp_2$Bone_meta = "(-)"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
    Total_pts_num = sum(Data_tmp$bone_meta)
    Data_tmp$bone_meta = fun_zero(Data_tmp$bone_meta, Total_pts_num)
    cancer_freq_order = c("(-)", "(+)")
    Data_tmp <- transform(Data_tmp, Bone_meta= factor(Bone_meta, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_8 = ggplot(data = Data_tmp, aes(x = Diseases, y = bone_meta, fill=Bone_meta)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Bone_meta))])) +
      ylab("Bone metastasis") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$lung_meta = Summary_Table$lung
    Data_tmp_1$Lung_meta = "(+)"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$lung_meta = Summary_Table$lung_no
    Data_tmp_2$Lung_meta = "(-)"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
    Total_pts_num = sum(Data_tmp$lung_meta)
    Data_tmp$lung_meta = fun_zero(Data_tmp$lung_meta, Total_pts_num)
    cancer_freq_order = c("(-)", "(+)")
    Data_tmp <- transform(Data_tmp, Lung_meta= factor(Lung_meta, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_9 = ggplot(data = Data_tmp, aes(x = Diseases, y = lung_meta, fill=Lung_meta)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Lung_meta))])) +
      ylab("Lung metastasis") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp_1 =data.frame(Summary_Table$Diseases)
    Data_tmp_1$liver_meta = Summary_Table$liver
    Data_tmp_1$Liver_meta = "(+)"
    Data_tmp_2 =data.frame(Summary_Table$Diseases)
    Data_tmp_2$liver_meta = Summary_Table$liver_no
    Data_tmp_2$Liver_meta = "(-)"
    Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
    Total_pts_num = sum(Data_tmp$liver_meta)
    Data_tmp$liver_meta = fun_zero(Data_tmp$liver_meta, Total_pts_num)
    cancer_freq_order = c("(-)", "(+)")
    Data_tmp <- transform(Data_tmp, Liver_meta= factor(Liver_meta, levels = cancer_freq_order))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_10 = ggplot(data = Data_tmp, aes(x = Diseases, y = liver_meta, fill=Liver_meta)) +
      geom_col(position="fill") +
      coord_flip() +
      theme_classic() +
      theme(# legend.position = "none",
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,1), labels = percent) +
      ylab("Liver metastasis") +
      scale_fill_manual(values = rev(brewer.pal(5, "Paired")[1:length(unique(Data_tmp$Liver_meta))]), name = "Metastasis") +
      guides(fill = guide_legend(reverse = TRUE))

    Data_tmp =data.frame(Summary_Table$Diseases)
    Data_tmp$option = Summary_Table$option / 100
    Data_tmp$Option = "(+)"
    Data_tmp <- transform(Data_tmp, Option= factor(Option))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_11 = ggplot(data = Data_tmp, aes(x = Diseases, y = option, fill=Option)) +
      geom_point() +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,NA), labels = percent, name = "Pts with CTx recommendation")

    Data_tmp =data.frame(Summary_Table$Diseases)
    Data_tmp$treat = Summary_Table$treat / 100
    Data_tmp$Treat = "(+)"
    Data_tmp <- transform(Data_tmp, Treat= factor(Treat))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_12 = ggplot(data = Data_tmp, aes(x = Diseases, y = treat, fill=Treat)) +
      geom_point() +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(lim = c(0,NA), labels = percent, name = "Pts received recommended CTx")

    Data_tmp =data.frame(Summary_Table$Diseases)
    Data_tmp$time_before_CGP = Summary_Table$time_before_CGP
    Data_tmp$Time_before_CGP = "(+)"
    Data_tmp <- transform(Data_tmp, Time_before_CGP= factor(Time_before_CGP))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_13 = ggplot(data = Data_tmp, aes(x = Diseases, y = time_before_CGP, fill=Time_before_CGP)) +
      geom_col() +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(limits = c(0, NA)) +
      ylab("Median time from CTx to CGP (days)") +
      scale_fill_manual(values = brewer.pal(3, "Paired")[1])

    Data_tmp =data.frame(Summary_Table$Diseases)
    Data_tmp$time_after_CGP = Summary_Table$time_after_CGP
    Data_tmp$Time_after_CGP = "(+)"
    Data_tmp <- transform(Data_tmp, Time_after_CGP= factor(Time_after_CGP))
    Data_tmp <- transform(Data_tmp, Diseases= factor(Diseases, levels = sort(unique(Data_case_target$diagnosis), decreasing = TRUE)))

    chart_14 = ggplot(data = Data_tmp, aes(x = Diseases, y = time_after_CGP, fill=Time_after_CGP)) +
      geom_col() +
      coord_flip() +
      theme_classic() +
      theme(legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_y_continuous(limits = c(0, NA)) +
      ylab("Median time from CGP to death (days)") +
      scale_fill_manual(values = brewer.pal(3, "Paired")[1])
  })
  grid.arrange(chart_3, chart_4, chart_5, chart_6, chart_7, chart_8, chart_9, chart_10, chart_11, chart_12, chart_13, chart_14, ncol = 4)
})

output$figure_drug_evidence = renderGirafe({
  req(input$figure_drug_evidence_var)
  req(OUTPUT_DATA$clustering_y)
  if(input$figure_drug_evidence_var == "Cancers"){
    req(OUTPUT_DATA$clustering_y)
    max_weight <- 1
    girafe(ggobj =
             ggplot(OUTPUT_DATA$clustering_y, aes(
               x = reorder(Cancers, desc(Cancers)),
               y = weight,
               fill = reorder(Level, desc(Level)),
               tooltip = tooltip,     # ツールチップを渡す
               data_id = paste0(Cancers, "_", Level)  # ハイライト用ID
             )) +
             geom_bar_interactive(stat = "identity") +
             geom_text(
               data = OUTPUT_DATA$clustering_y %>%
                 dplyr::group_by(Cancers) %>%
                 dplyr::summarise(total_patients = unique(Patients)),
               aes(
                 x = Cancers,
                 y = max_weight * 1.02,  # 横方向に102%の位置
                 label = total_patients
               ),
               inherit.aes = FALSE,
               hjust = 0                # 左寄せ
             ) +
             scale_y_continuous(
               limits = c(0, 1.15),                       # 割合なら0〜1
               breaks = seq(0, 1, 0.25),
               labels = scales::percent_format(accuracy = 1)
             ) +
             scale_fill_nejm() +
             theme_bw() +
             labs(
               y = "Frequency of druggable mutation",
               x = "Cancer type",
               fill = "Evidence level"
             ) +
             coord_flip() +
             guides(fill = guide_legend(reverse = TRUE)),
           width_svg = 8,
           height_svg = 0.6 + length(unique(OUTPUT_DATA$clustering_y$Cancers))*0.15,
           options = list(
             opts_tooltip(opacity = 0.9, css = "background-color:white; padding:5px; border:1px solid black;"),
             opts_hover(css = "fill-opacity:0.7;cursor:pointer;")
           ))
  } else if(input$figure_drug_evidence_var == "cluster"){
    req(OUTPUT_DATA$clustering_z)
    girafe(ggobj =
             ggplot(OUTPUT_DATA$clustering_z, aes(
               x = reorder(Clusters, desc(Clusters)),
               y = weight,
               fill = reorder(Level, desc(Level)),
               tooltip = tooltip,     # ツールチップを渡す
               data_id = paste0(Clusters, "_", Level)  # ハイライト用ID
             )) +
             geom_bar_interactive(stat = "identity") +
             geom_text(
               data = OUTPUT_DATA$clustering_z %>% dplyr::group_by(Clusters) %>%
                 dplyr::summarise(total_patients = unique(Patients)),
               aes(x = Clusters, y = 1, label = total_patients),
               inherit.aes = FALSE,
               hjust = -0.3) +
             expand_limits(y = 1.05) +
             scale_y_continuous(labels = percent) +
             scale_fill_nejm() +
             theme_bw() +
             labs(
               y = "Frequency of druggable mutation",
               x = "Cluster",
               fill = "Evidence level"
             ) +
             coord_flip() +
             guides(fill = guide_legend(reverse = TRUE)),
           width_svg = 8,
           height_svg = 0.6 + length(unique(OUTPUT_DATA$clustering_z$Clusters))*0.15,
           options = list(
             opts_tooltip(opacity = 0.9, css = "background-color:white; padding:5px; border:1px solid black;"),
             opts_hover(css = "fill-opacity:0.7;cursor:pointer;")
           ))
  }
})

output$table_basic_data_raw = DT::renderDataTable(server = FALSE,{
  req(OUTPUT_DATA$clustering_Data_evidence_table)
  create_datatable_with_confirm(OUTPUT_DATA$clustering_Data_evidence_table, "FELIS downloaded raw data in UMAP clustering based on mutations tab")
})

output$figure_cluster_subtype = renderGirafe({
  req(OUTPUT_DATA$clustering_Data_case_target,
      OUTPUT_DATA$clustering_selected_genes,
      OUTPUT_DATA$clustering_Data_mutation_cord,
      input$color_var_cluster)
  Data_case_target = OUTPUT_DATA$clustering_Data_case_target
  selected_genes = OUTPUT_DATA$clustering_selected_genes
  Data_case_target$EP_option = as.factor(Data_case_target$EP_option)
  Data_case_target$EP_treat = as.factor(Data_case_target$EP_treat)
  g = ggplot(OUTPUT_DATA$clustering_Data_mutation_cord,
             aes(V1,V2,tooltip = tooltip_text,color = Data_case_target[[input$color_var_cluster]])) +
    geom_point_interactive(aes(tooltip = tooltip_text), size = 1) +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = "Unsupervised clustering of oncogenic alterations") +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
  if(input$color_var_cluster == "cluster"){
    g = g +
      guides(color = guide_legend(title = "TOP3 histology/mutated gene", nrow=30)) +
      scale_color_discrete_interactive(labels = selected_genes)
  } else {
    g = g +
      guides(color = guide_legend(title = ifelse(input$color_var_cluster == "YoungOld", "Age",
                                                 ifelse(input$color_var_cluster == "Cancers", "Histology",
                                                        ifelse(input$color_var_cluster == "EP_option", "Treatment option recommended","Recommended treatment done")))))
  }
  g_main <- g + theme(legend.position = "none")
  legend <- cowplot::get_legend(g)
  combined <- cowplot::plot_grid(g_main, legend, rel_widths = c(1, 1), nrow = 1)
  girafe(ggobj = combined, width_svg = 16, height_svg = 8)
})
