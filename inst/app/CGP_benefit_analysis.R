CGP_benefit_analysis_logic <- function() {
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
                      症例.管理情報.登録日,
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
      Data_drug = Data_drug_raw()
      if(nrow(Data_case_target)>0){
        Data_case_target = Data_case_target %>%
          dplyr::distinct(.keep_all = TRUE, C.CAT調査結果.基本項目.ハッシュID)
        max_samples = isolate(input$nomogram_max_samples)
        max_samples = ifelse(is.null(max_samples), 10000, max_samples)
        if (nrow(Data_case_target) > max_samples){
          sample_indices = sample(nrow(Data_case_target),max_samples)
          Data_case_target = Data_case_target[sample_indices, ]
        } else {
          Data_case_target = Data_case_target
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
        ID_evidence_ABC = unique((Data_report() %>%
                                    dplyr::filter(
                                      Evidence_level %in% input$evidence_level) %>%
                                    dplyr::arrange(Tumor_Sample_Barcode))$Tumor_Sample_Barcode)
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

        Data_survival = Data_case_target %>%
          dplyr::filter(time_enroll_final > 0 & !is.na(time_enroll_final) & is.finite(time_enroll_final))
        Cancername = sort(unique(Data_survival$Cancers))
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
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
                        症例.背景情報.喫煙歴有無.名称.,
                        症例.背景情報.アルコール多飲有無.名称.,
                        症例.管理情報.登録日,
                        time_enroll_final,
                        time_diagnosis_enroll,
                        censor,
                        CTx_lines_before_CGP,
                        pre_CGP_best_RECIST,
                        症例.検体情報.パネル.名称.,
                        症例.基本情報.年齢)
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
            'Smoking_history', 'Alcoholic_history',
            'Enroll_date',
            'time_enroll_final',
            'time_diagnosis_enroll',
            'censor', "Lines", "Best_effect", "Panel","Age_raw")
        #Data_forest$Histology_dummy = paste0(Data_forest$Histology_dummy,"_NOS")
        Data_forest$time_diagnosis_enroll = ifelse(Data_forest$time_diagnosis_enroll > 365.25*1, ">1-year", "<=1-year")
        Data_forest$Enroll_date <- format(as.Date(Data_forest$Enroll_date), "%Y")
        col_added = 4
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
          dplyr::filter(PS != "Unknown" & Sex != "Unknown" & Smoking_history != "Unknown" & Alcoholic_history != "Unknown" & !is.na(time_diagnosis_enroll) & !is.na(Enroll_date)) %>%
          dplyr::mutate(
            PS = case_when(
              PS == 0 ~ "0",
              PS == 1 ~ "1",
              TRUE ~ "2_4"
            ),
            Panel = case_when(
              Panel %in% c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP") ~ Panel,
              TRUE ~ "Liquid"
            ),
            Lines = case_when(
              Lines == "0" ~ "0",
              Lines == "1" ~ "1",
              TRUE ~ "2~"
            )
          )
        Data_forest[Data_forest$Best_effect == "NE"]$Best_effect = "Unknown"
        Data_forest$Best_effect = factor(Data_forest$Best_effect, levels = c("SD","CR","PR","PD","Unknown","CTx_naive"))
        frequent_organs = sort(table(Data_forest$Histology),decreasing = T)
        frequent_organs_names = names(frequent_organs[frequent_organs>=max(50, input$minimum_pts)])
        Data_forest = Data_forest %>% dplyr::mutate(
          Histology = case_when(
            Histology %in% frequent_organs_names ~ Histology,
            TRUE ~ Histology_dummy
          )
        )
        oncogenic_genes = Data_MAF_target %>%
          dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
          dplyr::distinct() %>%
          dplyr::count(Hugo_Symbol) %>%
          dplyr::arrange(-n)
        colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
        oncogenic_genes = rbind(oncogenic_genes %>% dplyr::filter(gene_mutation %in% input$gene),
                                oncogenic_genes %>% dplyr::filter(!gene_mutation %in% input$gene))
        oncogenic_genes = oncogenic_genes[1:min(length(oncogenic_genes$all_patients), input$gene_no),]
        Gene_data = data.frame(oncogenic_genes$gene_mutation)
        colnames(Gene_data) = c("Gene")


        Data_forest =  Data_forest %>%
          dplyr::select(-Age_raw)
        gene_names = unique(c(input$gene,
                              as.character(Gene_data$Gene)))
        gene_names = gene_names[gene_names %in% unique(Data_MAF_target$Hugo_Symbol)]
        gene_names = gene_names[1:min(input$gene_no, length(gene_names))]
        for(i in 1:length(gene_names)){
          ID_mutation = unique((Data_MAF_target %>% dplyr::filter(Hugo_Symbol == gene_names[i]))$Tumor_Sample_Barcode)
          Data_forest_tmp = Data_forest %>%
            dplyr::select(ID) %>%
            dplyr::distinct() %>%
            dplyr::mutate(
              XXX = case_when(
                ID %in% ID_mutation ~ "mut(+)",
                TRUE ~ "mut(-)"
              )
            )
          colnames(Data_forest_tmp) = c("ID", gene_names[i])
          Data_forest = left_join(Data_forest, Data_forest_tmp, by=c("ID"))
        }
        incProgress(1 / 13)

        Data_forest_tmp__ = Data_forest %>%
          dplyr::select(-time_enroll_final, -censor, -EP_option, -ID)
        Disease_tmp = unique(sort(Data_forest_tmp__$Histology))
        for(i in 1:length(Disease_tmp)){
          Data_forest_tmp_2_ = Data_forest_tmp__ %>% dplyr::filter(Histology == Disease_tmp[i])
          if(sum(Data_forest_tmp_2_$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::mutate(Histology = case_when(
                Histology == Disease_tmp[i] ~ Histology_dummy,
                TRUE ~ Histology
              )
              )
          }
        }
        for(i in length(gene_names):1){
          kk = 14+col_added+i
          if(sum(Data_forest_tmp__[[kk]] == "mut(+)" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          } else if(sum(Data_forest_tmp__[[kk]] != "mut(+)" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          }
        }
        if(sum(Data_forest_tmp__$time_diagnosis_enroll == ">1-year" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        } else if(sum(Data_forest_tmp__$time_diagnosis_enroll != ">1-year" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        }
        if(sum(Data_forest_tmp__$PS == "0" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS != "0" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS == "2_4" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "1"] = "1_4"
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "2_4"] = "1_4"
        }
        if((sum(Data_forest_tmp__$Panel == "NCC OncoPanel" & Data_forest_tmp__$EP_treat == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "FoundationOne CDx" & Data_forest_tmp__$EP_treat == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "GenMineTOP" & Data_forest_tmp__$EP_treat == 1) < 3)){
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "NCC OncoPanel"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "FoundationOne CDx"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "GenMineTOP"] = "Solid"
          if(length(unique(Data_forest_tmp__$Panel))>1){
            Data_forest_tmp__$Panel = factor(Data_forest_tmp__$Panel)
            Data_forest_tmp__$Panel = relevel(Data_forest_tmp__$Panel, ref=names(sort(table(Data_forest_tmp__$Panel),decreasing = T))[[1]])
          }
          if(sum(Data_forest_tmp__$Panel == "Solid" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          } else if(sum(Data_forest_tmp__$Panel != "Solid" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          }
        }
        if("Panel" %in% colnames(Data_forest_tmp__)){
          Panel_data = data.frame(sort(unique(Data_forest_tmp__$Panel)))
          colnames(Panel_data) = c("Panel_type")
          save(file = file.path(tempdir(), "Panel_data.rda"), Panel_data)
        }
        if(sum(Data_forest_tmp__$Lines == "1" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines != "1" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines == "0" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "0"] = "0~1"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "0~1"
        } else if(sum(Data_forest_tmp__$Lines == "2~" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "2~"] = "1~"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "1~"
        }
        if(sum(Data_forest_tmp__$Enroll_date == "2019" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019"] = "2019/20"
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2020"] = "2019/20"
          if(sum(Data_forest_tmp__$Enroll_date == "2019/20" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019/20"] = "2019-21"
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2021"] = "2019-21"
            if(sum(Data_forest_tmp__$Enroll_date == "2019-21" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            } else if(sum(Data_forest_tmp__$Enroll_date != "2019-21" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            }
          } else if(sum(Data_forest_tmp__$Enroll_date != "2019/20" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Enroll_date)
          }
        } else if(sum(Data_forest_tmp__$Enroll_date != "2019" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Enroll_date)
        }
        if("Enroll_date" %in% colnames(Data_forest_tmp__)){
          Enroll_date_data = data.frame(sort(unique(Data_forest_tmp__$Enroll_date)))
          colnames(Enroll_date_data) = c("Enroll_date_type")
          save(file = file.path(tempdir(), "Enroll_date_data.rda"), Enroll_date_data)
        }
        if(sum(Data_forest_tmp__$Best_effect == "SD" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        } else if(sum(Data_forest_tmp__$Best_effect != "SD" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        }
        if(sum(Data_forest_tmp__$Age == "Younger" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        } else if(sum(Data_forest_tmp__$Age == "Older" & Data_forest_tmp__$EP_treat == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        }
        if(sum(Data_forest_tmp__$Sex == "Male" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        } else if(sum(Data_forest_tmp__$Sex != "Male" & Data_forest_tmp__$EP_treat == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        }
        metastasis_columns <- c("Smoking_history", "Alcoholic_history", "Lymph_met", "Lung_met", "Brain_met", "Bone_met", "Liver_met")
        metastasis_columns = metastasis_columns[metastasis_columns %in% Data_forest_tmp__]

        # 削除対象の列名を取得（EP_treat == 1 の行において、"Yes" の数または "Yes" でない数が3未満なら削除）
        if(length(metastasis_columns > 0)){
          cols_to_drop <- metastasis_columns[
            sapply(metastasis_columns, function(col) {
              yes_count <- Data_forest_tmp__[EP_treat == 1, sum(get(col) == "Yes", na.rm = TRUE)]
              no_count  <- Data_forest_tmp__[EP_treat == 1, sum(get(col) != "Yes", na.rm = TRUE)]
              (yes_count < 3 || no_count < 3)
            })
          ]
        } else {
          cols_to_drop = NULL
        }

        # 削除対象の列がある場合、一括削除
        if (length(cols_to_drop) > 0) {
          Data_forest_tmp__[, (cols_to_drop) := NULL]
        }
        if(input$HER2 != "No"){
          if(sum(Data_forest_tmp__$HER2_IHC == "Positive" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          } else if(sum(Data_forest_tmp__$HER2_IHC != "Positive" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          }
        }
        if(input$MSI != "No"){
          if(sum(Data_forest_tmp__$MSI_PCR == "Positive" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          } else if(sum(Data_forest_tmp__$MSI_PCR != "Positive" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          }
        }
        if(input$MMR != "No"){
          if(sum(Data_forest_tmp__$MMR_IHC == "dMMR" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          } else if(sum(Data_forest_tmp__$MMR_IHC != "dMMR" & Data_forest_tmp__$EP_treat == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          }
        }
        Data_forest_tmp__$Histology = factor(Data_forest_tmp__$Histology)
        if(length(unique(Data_forest_tmp__$Histology))<2){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Histology)
        } else {
          Data_forest_tmp__$Histology = relevel(Data_forest_tmp__$Histology,
                                                ref=names(sort(table(Data_forest_tmp__$Histology),decreasing = T))[[1]])
        }
        Data_forest_tmp__ = Data_forest_tmp__ %>%
          dplyr::select(-Histology_dummy)

        if(length(colnames(Data_forest_tmp__)) > 1){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_treat"]
          for(i in length(Factor_names_tmp__):1){
            if(all(as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]]) ==
                   as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]][1]))){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Factor_names_tmp__[i])
            }
          }
        }
        Data_forest_tmp__univariant = Data_forest_tmp__
        if(length(colnames(Data_forest_tmp__)) > 2){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_treat"]
          for(i in 1:(length(Factor_names_tmp__) - 1)){
            for(j in (i+1):length(Factor_names_tmp__)){
              if(Factor_names_tmp__[i] %in% colnames(Data_forest_tmp__) &
                 Factor_names_tmp__[j] %in% colnames(Data_forest_tmp__)){
                if(abs(cor(as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[i]]])),
                           as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[j]]])))) > 0.95){
                  Data_forest_tmp__ = Data_forest_tmp__ %>%
                    dplyr::select(-Factor_names_tmp__[j])
                }
              }
            }
          }
        }
        incProgress(1 / 13)

        Data_forest_tmp_7 = data.frame(Data_forest_tmp__)
        Data_forest_tmp_7__univariant = data.frame(Data_forest_tmp__univariant)
        Data_forest_tmp_7 = Data_forest_tmp_7[,colnames(Data_forest_tmp_7) %in%
                                                c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_treat","Panel",
                                                  "HER2_IHC", "MSI_PCR", "MMR_IHC","time_diagnosis_enroll", "Smoking_history", "Alcoholic_history"), drop=FALSE]
        Data_forest_tmp_7__univariant = Data_forest_tmp_7__univariant[,colnames(Data_forest_tmp_7__univariant) %in%
                                                                        c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_treat","Panel",
                                                                          "HER2_IHC", "MSI_PCR", "MMR_IHC","time_diagnosis_enroll", "Smoking_history", "Alcoholic_history", "Enroll_date"), drop=FALSE]
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant = Data_forest_tmp_7__univariant

        Data_forest_tmp__ = Data_forest %>%
          dplyr::select(-time_enroll_final, -censor, -EP_treat, -ID)
        Disease_tmp = unique(sort(Data_forest_tmp__$Histology))
        for(i in 1:length(Disease_tmp)){
          Data_forest_tmp_2_ = Data_forest_tmp__ %>% dplyr::filter(Histology == Disease_tmp[i])
          if(sum(Data_forest_tmp_2_$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::mutate(Histology = case_when(
                Histology == Disease_tmp[i] ~ Histology_dummy,
                TRUE ~ Histology
              )
              )
          }
        }
        for(i in length(gene_names):1){
          kk = 14+col_added+i
          if(sum(Data_forest_tmp__[[kk]] == "mut(+)" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          } else if(sum(Data_forest_tmp__[[kk]] != "mut(+)" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          }
        }

        if(sum(Data_forest_tmp__$time_diagnosis_enroll == ">1-year" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        } else if(sum(Data_forest_tmp__$time_diagnosis_enroll != ">1-year" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        }

        if(sum(Data_forest_tmp__$PS == "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS != "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS == "2_4" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "1"] = "1_4"
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "2_4"] = "1_4"
        }
        if((sum(Data_forest_tmp__$Panel == "NCC OncoPanel" & Data_forest_tmp__$EP_option == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "FoundationOne CDx" & Data_forest_tmp__$EP_option == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "GenMineTOP" & Data_forest_tmp__$EP_option == 1) < 3)){
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "NCC OncoPanel"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "FoundationOne CDx"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "GenMineTOP"] = "Solid"
          if(length(unique(Data_forest_tmp__$Panel))>1){
            Data_forest_tmp__$Panel = factor(Data_forest_tmp__$Panel)
            Data_forest_tmp__$Panel = relevel(Data_forest_tmp__$Panel, ref=names(sort(table(Data_forest_tmp__$Panel),decreasing = T))[[1]])
          }
          if(sum(Data_forest_tmp__$Panel == "Solid" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          } else if(sum(Data_forest_tmp__$Panel != "Solid" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          }
        }
        if("Panel" %in% colnames(Data_forest_tmp__)){
          Panel_data_rec = data.frame(sort(unique(Data_forest_tmp__$Panel)))
          colnames(Panel_data_rec) = c("Panel_type")
          save(file = file.path(tempdir(), "Panel_data_rec.rda"), Panel_data_rec)
        }

        if(sum(Data_forest_tmp__$Enroll_date == "2019" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019"] = "2019/20"
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2020"] = "2019/20"
          if(sum(Data_forest_tmp__$Enroll_date == "2019/20" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019/20"] = "2019-21"
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2021"] = "2019-21"
            if(sum(Data_forest_tmp__$Enroll_date == "2019-21" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            } else if(sum(Data_forest_tmp__$Enroll_date != "2019-21" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            }
          } else if(sum(Data_forest_tmp__$Enroll_date != "2019/20" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Enroll_date)
          }
        } else if(sum(Data_forest_tmp__$Enroll_date != "2019" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Enroll_date)
        }
        if("Enroll_date" %in% colnames(Data_forest_tmp__)){
          Enroll_date_data_rec = data.frame(sort(unique(Data_forest_tmp__$Enroll_date)))
          colnames(Enroll_date_data_rec) = c("Enroll_date_type")
          save(file = file.path(tempdir(), "Enroll_date_data_rec.rda"), Enroll_date_data_rec)
        }
        if(sum(Data_forest_tmp__$Lines == "1" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines != "1" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines == "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "0"] = "0~1"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "0~1"
        } else if(sum(Data_forest_tmp__$Lines == "2~" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "2~"] = "1~"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "1~"
        }
        if(sum(Data_forest_tmp__$Best_effect == "SD" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        } else if(sum(Data_forest_tmp__$Best_effect != "SD" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        }
        if(sum(Data_forest_tmp__$Age == "Younger" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        } else if(sum(Data_forest_tmp__$Age == "Older" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        }
        if(sum(Data_forest_tmp__$Sex == "Male" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        } else if(sum(Data_forest_tmp__$Sex != "Male" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        }
        metastasis_columns <- c("Smoking_history", "Alcoholic_history", "Lymph_met", "Lung_met", "Brain_met", "Bone_met", "Liver_met")
        metastasis_columns = metastasis_columns[metastasis_columns %in% Data_forest_tmp__]

        # 削除対象の列名を取得（EP_option == 1 の行において、"Yes" の数または "Yes" でない数が3未満なら削除）
        if(length(metastasis_columns > 0)){
          cols_to_drop <- metastasis_columns[
            sapply(metastasis_columns, function(col) {
              yes_count <- Data_forest_tmp__[EP_option == 1, sum(get(col) == "Yes", na.rm = TRUE)]
              no_count  <- Data_forest_tmp__[EP_option == 1, sum(get(col) != "Yes", na.rm = TRUE)]
              (yes_count < 3 || no_count < 3)
            })
          ]
        } else {
          cols_to_drop = NULL
        }

        # 削除対象の列がある場合、一括削除
        if (length(cols_to_drop) > 0) {
          Data_forest_tmp__[, (cols_to_drop) := NULL]
        }
        if(input$HER2 != "No"){
          if(sum(Data_forest_tmp__$HER2_IHC == "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          } else if(sum(Data_forest_tmp__$HER2_IHC != "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          }
        }
        if(input$MSI != "No"){
          if(sum(Data_forest_tmp__$MSI_PCR == "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          } else if(sum(Data_forest_tmp__$MSI_PCR != "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          }
        }
        if(input$MMR != "No"){
          if(sum(Data_forest_tmp__$MMR_IHC == "dMMR" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          } else if(sum(Data_forest_tmp__$MMR_IHC != "dMMR" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          }
        }
        Data_forest_tmp__$Histology = factor(Data_forest_tmp__$Histology)
        if(length(unique(Data_forest_tmp__$Histology))<2){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Histology)
        } else {
          Data_forest_tmp__$Histology = relevel(Data_forest_tmp__$Histology,
                                                ref=names(sort(table(Data_forest_tmp__$Histology),decreasing = T))[[1]])
        }
        Data_forest_tmp__ = Data_forest_tmp__ %>%
          dplyr::select(-Histology_dummy)

        if(length(colnames(Data_forest_tmp__)) > 1){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_option"]
          for(i in length(Factor_names_tmp__):1){
            if(all(as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]]) ==
                   as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]][1]))){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Factor_names_tmp__[i])
            }
          }
        }
        Data_forest_tmp__univariant = Data_forest_tmp__
        if(length(colnames(Data_forest_tmp__)) > 2){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_option"]
          for(i in 1:(length(Factor_names_tmp__) - 1)){
            for(j in (i+1):length(Factor_names_tmp__)){
              if(Factor_names_tmp__[i] %in% colnames(Data_forest_tmp__) &
                 Factor_names_tmp__[j] %in% colnames(Data_forest_tmp__)){
                if(abs(cor(as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[i]]])),
                           as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[j]]])))) > 0.95){
                  Data_forest_tmp__ = Data_forest_tmp__ %>%
                    dplyr::select(-Factor_names_tmp__[j])
                }
              }
            }
          }
        }
        Data_forest_tmp_7_rec = data.frame(Data_forest_tmp__)
        Data_forest_tmp_7__univariant_rec = data.frame(Data_forest_tmp__univariant)
        Data_forest_tmp_7_rec = Data_forest_tmp_7_rec[,colnames(Data_forest_tmp_7_rec) %in%
                                                        c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_option","Panel",
                                                          "HER2_IHC", "MSI_PCR", "MMR_IHC", "Smoking_history", "Alcoholic_history", "time_diagnosis_enroll"), drop=FALSE]
        Data_forest_tmp_7__univariant_rec = Data_forest_tmp_7__univariant_rec[,colnames(Data_forest_tmp_7__univariant_rec) %in%
                                                                                c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_option","Panel",
                                                                                  "HER2_IHC", "MSI_PCR", "MMR_IHC", "Smoking_history", "Alcoholic_history", "time_diagnosis_enroll", "Enroll_date"), drop=FALSE]
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant_rec = Data_forest_tmp_7__univariant_rec

        incProgress(1 / 13)

        Data_forest_tmp__ = Data_forest %>%
          dplyr::mutate(EP_option = case_when(
            ID %in% ID_evidence_ABC ~ 1,
            TRUE ~ 0
          )) %>%
          dplyr::select(-time_enroll_final, -censor, -EP_treat, -ID)
        Disease_tmp = unique(sort(Data_forest_tmp__$Histology))
        for(i in 1:length(Disease_tmp)){
          Data_forest_tmp_2_ = Data_forest_tmp__ %>% dplyr::filter(Histology == Disease_tmp[i])
          if(sum(Data_forest_tmp_2_$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::mutate(Histology = case_when(
                Histology == Disease_tmp[i] ~ Histology_dummy,
                TRUE ~ Histology
              )
              )
          }
        }
        for(i in length(gene_names):1){
          kk = 14+col_added+i
          if(sum(Data_forest_tmp__[[kk]] == "mut(+)" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          } else if(sum(Data_forest_tmp__[[kk]] != "mut(+)" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__[,-..kk]
          }
        }
        if(sum(Data_forest_tmp__$PS == "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS != "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-PS)
        } else if(sum(Data_forest_tmp__$PS == "2_4" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "1"] = "1_4"
          Data_forest_tmp__$PS[Data_forest_tmp__$PS == "2_4"] = "1_4"
        }
        if(sum(Data_forest_tmp__$time_diagnosis_enroll == ">1-year" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        } else if(sum(Data_forest_tmp__$time_diagnosis_enroll != ">1-year" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-time_diagnosis_enroll)
        }
        if((sum(Data_forest_tmp__$Panel == "NCC OncoPanel" & Data_forest_tmp__$EP_option == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "FoundationOne CDx" & Data_forest_tmp__$EP_option == 1) < 3) |
           (sum(Data_forest_tmp__$Panel == "GenMineTOP" & Data_forest_tmp__$EP_option == 1) < 3)){
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "NCC OncoPanel"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "FoundationOne CDx"] = "Solid"
          Data_forest_tmp__$Panel[Data_forest_tmp__$Panel == "GenMineTOP"] = "Solid"
          if(length(unique(Data_forest_tmp__$Panel))>1){
            Data_forest_tmp__$Panel = factor(Data_forest_tmp__$Panel)
            Data_forest_tmp__$Panel = relevel(Data_forest_tmp__$Panel, ref=names(sort(table(Data_forest_tmp__$Panel),decreasing = T))[[1]])
          }
          if(sum(Data_forest_tmp__$Panel == "Solid" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          } else if(sum(Data_forest_tmp__$Panel != "Solid" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Panel)
          }
        }
        if("Panel" %in% colnames(Data_forest_tmp__)){
          Panel_data_GMT = data.frame(sort(unique(Data_forest_tmp__$Panel)))
          colnames(Panel_data_GMT) = c("Panel_type")
          save(file = file.path(tempdir(), "Panel_data_GMT.rda"), Panel_data_GMT)
        }
        if(sum(Data_forest_tmp__$Enroll_date == "2019" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019"] = "2019/20"
          Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2020"] = "2019/20"
          if(sum(Data_forest_tmp__$Enroll_date == "2019/20" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2019/20"] = "2019-21"
            Data_forest_tmp__$Enroll_date[Data_forest_tmp__$Enroll_date == "2021"] = "2019-21"
            if(sum(Data_forest_tmp__$Enroll_date == "2019-21" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            } else if(sum(Data_forest_tmp__$Enroll_date != "2019-21" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Enroll_date)
            }
          } else if(sum(Data_forest_tmp__$Enroll_date != "2019/20" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-Enroll_date)
          }
        } else if(sum(Data_forest_tmp__$Enroll_date != "2019" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Enroll_date)
        }
        if("Enroll_date" %in% colnames(Data_forest_tmp__)){
          Enroll_date_data_GMT = data.frame(sort(unique(Data_forest_tmp__$Enroll_date)))
          colnames(Enroll_date_data_GMT) = c("Enroll_date_type")
          save(file = file.path(tempdir(), "Enroll_date_data_GMT.rda"), Enroll_date_data_GMT)
        }

        if(sum(Data_forest_tmp__$Lines == "1" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines != "1" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Lines)
        } else if(sum(Data_forest_tmp__$Lines == "0" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "0"] = "0~1"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "0~1"
        } else if(sum(Data_forest_tmp__$Lines == "2~" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "2~"] = "1~"
          Data_forest_tmp__$Lines[Data_forest_tmp__$Lines == "1"] = "1~"
        }
        if(sum(Data_forest_tmp__$Best_effect == "SD" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        } else if(sum(Data_forest_tmp__$Best_effect != "SD" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Best_effect)
        }
        if(sum(Data_forest_tmp__$Age == "Younger" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        } else if(sum(Data_forest_tmp__$Age == "Older" & Data_forest_tmp__$EP_option == 1) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Age)
        }
        if(sum(Data_forest_tmp__$Sex == "Male" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        } else if(sum(Data_forest_tmp__$Sex != "Male" & Data_forest_tmp__$EP_option == 1,na.rm = T) < 3){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Sex)
        }
        metastasis_columns <- c("Smoking_history", "Alcoholic_history", "Lymph_met", "Lung_met", "Brain_met", "Bone_met", "Liver_met")
        metastasis_columns = metastasis_columns[metastasis_columns %in% Data_forest_tmp__]

        # 削除対象の列名を取得（EP_option == 1 の行において、"Yes" の数または "Yes" でない数が3未満なら削除）
        if(length(metastasis_columns > 0)){
          cols_to_drop <- metastasis_columns[
            sapply(metastasis_columns, function(col) {
              yes_count <- Data_forest_tmp__[EP_option == 1, sum(get(col) == "Yes", na.rm = TRUE)]
              no_count  <- Data_forest_tmp__[EP_option == 1, sum(get(col) != "Yes", na.rm = TRUE)]
              (yes_count < 3 || no_count < 3)
            })
          ]
        } else {
          cols_to_drop = NULL
        }

        # 削除対象の列がある場合、一括削除
        if (length(cols_to_drop) > 0) {
          Data_forest_tmp__[, (cols_to_drop) := NULL]
        }
        if(input$HER2 != "No"){
          if(sum(Data_forest_tmp__$HER2_IHC == "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          } else if(sum(Data_forest_tmp__$HER2_IHC != "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-HER2_IHC)
          }
        }
        if(input$MSI != "No"){
          if(sum(Data_forest_tmp__$MSI_PCR == "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          } else if(sum(Data_forest_tmp__$MSI_PCR != "Positive" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MSI_PCR)
          }
        }
        if(input$MMR != "No"){
          if(sum(Data_forest_tmp__$MMR_IHC == "dMMR" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          } else if(sum(Data_forest_tmp__$MMR_IHC != "dMMR" & Data_forest_tmp__$EP_option == 1) < 3){
            Data_forest_tmp__ = Data_forest_tmp__ %>%
              dplyr::select(-MMR_IHC)
          }
        }
        Data_forest_tmp__$Histology = factor(Data_forest_tmp__$Histology)
        if(length(unique(Data_forest_tmp__$Histology))<2){
          Data_forest_tmp__ = Data_forest_tmp__ %>%
            dplyr::select(-Histology)
        } else {
          Data_forest_tmp__$Histology = relevel(Data_forest_tmp__$Histology,
                                                ref=names(sort(table(Data_forest_tmp__$Histology),decreasing = T))[[1]])
        }
        Data_forest_tmp__ = Data_forest_tmp__ %>%
          dplyr::select(-Histology_dummy)

        if(length(colnames(Data_forest_tmp__)) > 1){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_option"]
          for(i in length(Factor_names_tmp__):1){
            if(all(as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]]) ==
                   as.character(Data_forest_tmp__[[Factor_names_tmp__[i]]][1]))){
              Data_forest_tmp__ = Data_forest_tmp__ %>%
                dplyr::select(-Factor_names_tmp__[i])
            }
          }
        }
        Data_forest_tmp__univariant = Data_forest_tmp__
        if(length(colnames(Data_forest_tmp__)) > 2){
          Factor_names_tmp__ = colnames(Data_forest_tmp__)
          Factor_names_tmp__ = Factor_names_tmp__[Factor_names_tmp__ != "EP_option"]
          for(i in 1:(length(Factor_names_tmp__) - 1)){
            for(j in (i+1):length(Factor_names_tmp__)){
              if(Factor_names_tmp__[i] %in% colnames(Data_forest_tmp__) &
                 Factor_names_tmp__[j] %in% colnames(Data_forest_tmp__)){
                if(abs(cor(as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[i]]])),
                           as.numeric(as.factor(Data_forest_tmp__[[Factor_names_tmp__[j]]])))) > 0.95){
                  Data_forest_tmp__ = Data_forest_tmp__ %>%
                    dplyr::select(-Factor_names_tmp__[j])
                }
              }
            }
          }
        }
        incProgress(1 / 13)

        Data_forest_tmp_7_GMT = data.frame(Data_forest_tmp__)
        Data_forest_tmp_7__univariant_GMT = data.frame(Data_forest_tmp__univariant)
        Data_forest_tmp_7_GMT = Data_forest_tmp_7_GMT[,colnames(Data_forest_tmp_7_GMT) %in%
                                                        c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_option","Panel",
                                                          "HER2_IHC", "MSI_PCR", "MMR_IHC","Smoking_history","Alcoholic_history","time_diagnosis_enroll"), drop=FALSE]
        Data_forest_tmp_7__univariant_GMT = Data_forest_tmp_7__univariant_GMT[,colnames(Data_forest_tmp_7__univariant_GMT) %in%
                                                                                c("Age","Lymph_met","Brain_met","Lung_met","Bone_met","Liver_met","Sex","Histology","PS","Lines","Best_effect","EP_option","Panel",
                                                                                  "HER2_IHC", "MSI_PCR", "MMR_IHC","Alcoholic_history","time_diagnosis_enroll", "Enroll_date"), drop=FALSE]
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant_GMT = Data_forest_tmp_7__univariant_GMT
        Data_forest_tmp_ = Data_forest  %>%
          dplyr::filter(EP_option == 1) %>%
          dplyr::select(-time_enroll_final, -censor, -EP_option, -ID)
        Disease_tmp = unique(sort(Data_forest_tmp_$Histology))
        for(i in 1:length(Disease_tmp)){
          Data_forest_tmp_2 = Data_forest_tmp_ %>% dplyr::filter(Histology == Disease_tmp[i])
          if(sum(Data_forest_tmp_2$EP_treat == 1) < 3){
            Data_forest_tmp_ = Data_forest_tmp_ %>%
              dplyr::mutate(Histology = case_when(
                Histology == Disease_tmp[i] ~ Histology_dummy,
                TRUE ~ Histology
              )
              )
          }
        }
        for(i in length(gene_names):1){
          kk = 14+col_added+i
          if(sum(Data_forest_tmp_[[kk]] == "mut(+)" & Data_forest_tmp_$EP_treat == 1) < 3){
            Data_forest_tmp_ = Data_forest_tmp_[,-..kk]
          } else if(sum(Data_forest_tmp_[[kk]] != "mut(+)" & Data_forest_tmp_$EP_treat == 1) < 3){
            Data_forest_tmp_ = Data_forest_tmp_[,-..kk]
          }
        }
        Data_forest_tmp_ = optimize_data_datatable(Data_forest_tmp_, input)
        Data_forest_tmp_$Histology = factor(Data_forest_tmp_$Histology)
        if(length(unique(Data_forest_tmp_$Histology))<2){
          Data_forest_tmp_ = Data_forest_tmp_ %>%
            dplyr::select(-Histology)
        } else {
          Data_forest_tmp_$Histology = relevel(Data_forest_tmp_$Histology,
                                               ref=names(sort(table(Data_forest_tmp_$Histology),decreasing = T))[[1]])
        }

        Data_forest_tmp_ = Data_forest_tmp_ %>%
          dplyr::select(-Histology_dummy)
        if(length(colnames(Data_forest_tmp_)) > 1){
          Factor_names_tmp_ = colnames(Data_forest_tmp_)
          Factor_names_tmp_ = Factor_names_tmp_[Factor_names_tmp_ != "EP_treat"]
          for(i in length(Factor_names_tmp_):1){
            if(all(as.character(Data_forest_tmp_[[Factor_names_tmp_[i]]]) ==
                   as.character(Data_forest_tmp_[[Factor_names_tmp_[i]]][1]))){
              Data_forest_tmp_ = Data_forest_tmp_ %>%
                dplyr::select(-Factor_names_tmp_[i])
            }
          }
        }

        Data_forest_tmp_univariant = Data_forest_tmp_

        if(length(colnames(Data_forest_tmp_)) > 2){
          Factor_names_tmp_ = colnames(Data_forest_tmp_)
          Factor_names_tmp_ = Factor_names_tmp_[Factor_names_tmp_ != "EP_treat"]
          for(i in 1:(length(Factor_names_tmp_) - 1)){
            for(j in (i+1):length(Factor_names_tmp_)){
              if(Factor_names_tmp_[i] %in% colnames(Data_forest_tmp_) &
                 Factor_names_tmp_[j] %in% colnames(Data_forest_tmp_)){
                if(abs(cor(as.numeric(as.factor(Data_forest_tmp_[[Factor_names_tmp_[i]]])),
                           as.numeric(as.factor(Data_forest_tmp_[[Factor_names_tmp_[j]]])))) > 0.95){
                  Data_forest_tmp_ = Data_forest_tmp_ %>%
                    dplyr::select(-Factor_names_tmp_[j])
                }
              }
            }
          }
        }
        incProgress(1 / 13)

        Data_forest_tmp_6 = Data_forest_tmp_
        Data_forest_tmp_6_univariant = Data_forest_tmp_univariant
        OUTPUT_DATA$table_prediction_Data_forest_tmp_6 = Data_forest_tmp_6
        OUTPUT_DATA$table_prediction_Data_forest_tmp_6_univariant = Data_forest_tmp_6_univariant
        Data_forest_tmp_8 = Data_forest_tmp_7
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7 = Data_forest_tmp_7
        if("Lines" %in% colnames(Data_forest_tmp_8)){
          Data_forest_tmp_8 = Data_forest_tmp_8 %>% select(-Lines)
        }

        Data_ML = Data_forest_tmp_8
        if("Histology" %in% colnames(Data_forest_tmp_8)){
          Organ_data = data.frame(Data_ML$Histology, Data_forest_tmp_8$Histology) %>% dplyr::distinct()
        } else {
          Organ_data = data.frame("not significant", "not significant")
        }
        colnames(Organ_data) = c("Organ_type_raw", "Organ_type")
        save(file = file.path(tempdir(), "Organ_data.rda"), Organ_data)
        Data_forest_tmp_8 = data.frame( lapply(Data_forest_tmp_8, as.factor) )
        OUTPUT_DATA$table_prediction_Data_forest_tmp_8 = Data_forest_tmp_8
        if(length(colnames(Data_forest_tmp_8)) > 1){
          if("Histology" %in% colnames(Data_forest_tmp_8)){
            Data_forest_tmp_8$Histology <- relevel(Data_forest_tmp_8$Histology, ref=names(sort(table(Data_forest_tmp_8$Histology),decreasing = T))[[1]])
          }
          final_mv_regxlevels <- (prepare_data(as.numeric(as.character(Data_forest_tmp_8$EP_treat)),
                                               as.matrix(data.frame(lapply(Data_forest_tmp_8[, setdiff(names(Data_forest_tmp_8), "EP_treat")],
                                                                           function(x) as.numeric(as.factor(x))))), type = "logistic") %>% stepwise(aic))$model
          OUTPUT_DATA$table_prediction_final_mv_regxlevels = final_mv_regxlevels
          if(length(final_mv_regxlevels) != 0){
            Data_forest_tmp_treat = Data_forest_tmp_8 %>% dplyr::select(c("EP_treat", final_mv_regxlevels))
            OUTPUT_DATA$table_prediction_Data_forest_tmp_treat = Data_forest_tmp_treat
            if(length(colnames(Data_forest_tmp_treat)) > 1){
              dd <- datadist(Data_forest_tmp_treat)
              options(datadist=dd)
              formula_treat = as.formula(paste0("EP_treat ~ `",paste(colnames(Data_forest_tmp_treat)[colnames(Data_forest_tmp_treat)!="EP_treat"], collapse="` + `"),"`"))
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "m2.qs"))){
                m2 = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "m2.qs"))
              } else {
                m2 <- lrm(formula_treat,
                          data = Data_forest_tmp_treat,x=TRUE,y=TRUE,penalty = Penal, tol=Toler)
              }
              BootNo = BootNoSet(Data_forest_tmp_treat)
              OUTPUT_DATA$table_prediction_BootNo <<- BootNo
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "val.qs"))){
                val = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "val.qs"))
              } else {
                # error in rms package
                val <- try(rms::validate(m2, B=BootNo), silent = FALSE)
                if (class(val) == "try-error") {
                  val <- try(rms::validate(m2, B=BootNo), silent = FALSE)
                  if (class(val) == "try-error") {
                    return(NULL)
                  }
                }
                # do not delete
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, val, file = file.path(tempdir(), "val.qs"))
                }
              }
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "cal.qs"))){
                cal = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "cal.qs"))
              } else {
                cal <- try(rms::calibrate(m2, B=BootNo), silent = FALSE)
                if (class(cal) == "try-error") {
                  cal = try(rms::calibrate(m2, B=BootNo), silent = FALSE)
                  if (class(cal) == "try-error") {
                    return(NULL)
                  }
                }
                OUTPUT_DATA$table_prediction_cal <<- cal
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, cal, file = file.path(tempdir(), "cal.qs"))
                }
              }
              options(datadist=NULL)
              log_nomogram <- nomogram(m2, fun = plogis, lp=F, funlabel="Matched therapy")
              OUTPUT_DATA$table_prediction_log_nomogram = log_nomogram
              OUTPUT_DATA$table_prediction_m2_pvalue <<- m2$stats["P"][[1]]
              OUTPUT_DATA$table_prediction_m2_val_1 <<- val[9,5]
              OUTPUT_DATA$table_prediction_m2_val_2 <<- val[2,5]
              OUTPUT_DATA$table_prediction_m2_val_3 <<- val[1,5]
              OUTPUT_DATA$table_prediction_m2_val_4 <<- val[3,5]
              OUTPUT_DATA$table_prediction_m2_val_5 <<- val[4,5]
              m2$x <- NULL
              m2$y <- NULL
              m2$linear.predictors <- NULL
              m2$fitted.values <- NULL
              m2$residuals <- NULL
              m2$weights <- NULL
              m2$na.action <- NULL
              QS_SAVE(nthreads = PARALLEL, m2, file = file.path(tempdir(), "m2.qs"))
              rm(m2)
              rm(val)
            }
          }
        }
        incProgress(1 / 13)

        Data_ML = Data_forest_tmp_7_rec
        Data_forest_tmp_8_rec = Data_forest_tmp_7_rec
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7_rec = Data_forest_tmp_7_rec
        if("Lines" %in% colnames(Data_forest_tmp_8_rec)){
          Data_forest_tmp_8_rec = Data_forest_tmp_8_rec %>% select(-Lines)
        }
        if("Histology" %in% colnames(Data_forest_tmp_8_rec)){
          Organ_data_rec = data.frame(Data_ML$Histology, Data_forest_tmp_8_rec$Histology) %>% dplyr::distinct()
        } else {
          Organ_data_rec = data.frame("not significant", "not significant")
        }
        colnames(Organ_data_rec) = c("Organ_type_raw", "Organ_type")
        save(file = file.path(tempdir(), "Organ_data_rec.rda"), Organ_data_rec)
        Data_forest_tmp_8_rec = data.frame( lapply(Data_forest_tmp_8_rec, as.factor) )
        OUTPUT_DATA$table_prediction_Data_forest_tmp_8_rec = Data_forest_tmp_8_rec
        if(length(colnames(Data_forest_tmp_8_rec)) > 1){
          if("Histology" %in% colnames(Data_forest_tmp_8_rec)){
            Data_forest_tmp_8_rec$Histology <- relevel(Data_forest_tmp_8_rec$Histology, ref=names(sort(table(Data_forest_tmp_8_rec$Histology),decreasing = T))[[1]])
          }
          final_mv_regxlevels_rec <- (prepare_data(as.numeric(as.character(Data_forest_tmp_8_rec$EP_option)),
                                                   as.matrix(data.frame(lapply(Data_forest_tmp_8_rec[, setdiff(names(Data_forest_tmp_8_rec), "EP_option")],
                                                                               function(x) as.numeric(as.factor(x))))), type = "logistic") %>% stepwise(aic))$model
          OUTPUT_DATA$table_prediction_final_mv_regxlevels_rec = final_mv_regxlevels_rec
          if(length(final_mv_regxlevels_rec) != 0){
            Data_forest_tmp_treat_rec = Data_forest_tmp_8_rec %>% dplyr::select(c("EP_option", final_mv_regxlevels_rec))
            OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_rec = Data_forest_tmp_treat_rec
            if(length(colnames(Data_forest_tmp_treat_rec)) > 1){
              dd <- datadist(Data_forest_tmp_treat_rec)
              options(datadist=dd)
              formula_treat_rec = as.formula(paste0("EP_option ~ `",paste(colnames(Data_forest_tmp_treat_rec)[colnames(Data_forest_tmp_treat_rec)!="EP_option"], collapse="` + `"),"`"))
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "m2_rec.qs"))){
                m2_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "m2_rec.qs"))
              } else {
                m2_rec <- lrm(formula_treat_rec,
                              data = Data_forest_tmp_treat_rec,x=TRUE,y=TRUE,penalty = Penal, tol=Toler)
              }
              BootNo = BootNoSet(Data_forest_tmp_treat_rec)
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "val_rec.qs"))){
                val_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "val_rec.qs"))
              } else {
                # error in rms package
                val_rec <- try(rms::validate(m2_rec, B=BootNo), silent = FALSE)
                if (class(val_rec) == "try-error") {
                  val_rec = try(rms::validate(m2_rec, B=BootNo), silent = FALSE)
                  if (class(val_rec) == "try-error") {
                    return(NULL)
                  }
                }
                # do not delete
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, val_rec, file = file.path(tempdir(), "val_rec.qs"))
                }
              }
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "cal_rec.qs"))){
                cal_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "cal_rec.qs"))
              } else {
                cal_rec <- try(rms::calibrate(m2_rec, B=BootNo), silent = FALSE)
                if (class(cal_rec) == "try-error") {
                  cal_rec = try(rms::calibrate(m2_rec, B=BootNo), silent = FALSE)
                  if (class(cal_rec) == "try-error") {
                    return(NULL)
                  }
                }
                OUTPUT_DATA$table_prediction_cal_rec <<- cal_rec
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, cal_rec, file = file.path(tempdir(), "cal_rec.qs"))
                }
              }
              options(datadist=NULL)
              log_nomogram_rec <- nomogram(m2_rec, fun = plogis, lp=F, funlabel="Treatment recommendation")
              OUTPUT_DATA$table_prediction_log_nomogram_rec = log_nomogram_rec
              OUTPUT_DATA$table_prediction_m2_rec_pvalue <<- m2_rec$stats["P"][[1]]
              OUTPUT_DATA$table_prediction_m2_rec_val_1 <<- val_rec[9,5]
              OUTPUT_DATA$table_prediction_m2_rec_val_2 <<- val_rec[2,5]
              OUTPUT_DATA$table_prediction_m2_rec_val_3 <<- val_rec[1,5]
              OUTPUT_DATA$table_prediction_m2_rec_val_4 <<- val_rec[3,5]
              OUTPUT_DATA$table_prediction_m2_rec_val_5 <<- val_rec[4,5]
              m2_rec$x <- NULL
              m2_rec$y <- NULL
              m2_rec$linear.predictors <- NULL
              m2_rec$fitted.values <- NULL
              m2_rec$residuals <- NULL
              m2_rec$weights <- NULL
              m2_rec$na.action <- NULL
              QS_SAVE(nthreads = PARALLEL, m2_rec, file = file.path(tempdir(), "m2_rec.qs"))
              rm(m2_rec)
              rm(val_rec)
            }
          }
        }
        incProgress(1 / 13)

        Data_ML = Data_forest_tmp_7_GMT
        Data_forest_tmp_8_GMT = Data_forest_tmp_7_GMT
        OUTPUT_DATA$table_prediction_Data_forest_tmp_7_GMT = Data_forest_tmp_7_GMT
        if("Lines" %in% colnames(Data_forest_tmp_8_GMT)){
          Data_forest_tmp_8_GMT = Data_forest_tmp_8_GMT %>% select(-Lines)
        }
        if("Histology" %in% colnames(Data_forest_tmp_8_GMT)){
          Organ_data_GMT = data.frame(Data_ML$Histology, Data_forest_tmp_8_GMT$Histology) %>% dplyr::distinct()
        } else {
          Organ_data_GMT = data.frame("not significant", "not significant")
        }
        colnames(Organ_data_GMT) = c("Organ_type_raw", "Organ_type")
        save(file = file.path(tempdir(), "Organ_data_GMT.rda"), Organ_data_GMT)
        Data_forest_tmp_8_GMT = data.frame( lapply(Data_forest_tmp_8_GMT, as.factor) )
        OUTPUT_DATA$table_prediction_Data_forest_tmp_8_GMT = Data_forest_tmp_8_GMT
        if(length(colnames(Data_forest_tmp_8_GMT)) > 1){
          if("Histology" %in% colnames(Data_forest_tmp_8_GMT)){
            Data_forest_tmp_8_GMT$Histology <- relevel(Data_forest_tmp_8_GMT$Histology, ref=names(sort(table(Data_forest_tmp_8_GMT$Histology),decreasing = T))[[1]])
          }
          final_mv_regxlevels_GMT <- (prepare_data(as.numeric(as.character(Data_forest_tmp_8_GMT$EP_option)),
                                                   as.matrix(data.frame(lapply(Data_forest_tmp_8_GMT[, setdiff(names(Data_forest_tmp_8_GMT), "EP_option")],
                                                                               function(x) as.numeric(as.factor(x))))), type = "logistic") %>% stepwise(aic))$model
          OUTPUT_DATA$table_prediction_final_mv_regxlevels_GMT = final_mv_regxlevels_GMT
          if(length(final_mv_regxlevels_GMT) != 0){
            Data_forest_tmp_treat_GMT = Data_forest_tmp_8_GMT %>% dplyr::select(c("EP_option", final_mv_regxlevels_GMT))
            OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_GMT = Data_forest_tmp_treat_GMT
            if(length(colnames(Data_forest_tmp_treat_GMT)) > 1){
              dd <- datadist(Data_forest_tmp_treat_GMT)
              options(datadist=dd)
              formula_treat_GMT = as.formula(paste0("EP_option ~ `",paste(colnames(Data_forest_tmp_treat_GMT)[colnames(Data_forest_tmp_treat_GMT)!="EP_option"], collapse="` + `"),"`"))
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "m2_GMT.qs"))){
                m2_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "m2_GMT.qs"))
              } else {
                m2_GMT <- lrm(formula_treat_GMT,
                              data = Data_forest_tmp_treat_GMT,x=TRUE,y=TRUE,penalty = Penal, tol=Toler)
              }
              BootNo = BootNoSet(Data_forest_tmp_treat_GMT)
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "val_GMT.qs"))){
                val_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "val_GMT.qs"))
              } else {
                # error in rms package
                val_GMT <- try(rms::validate(m2_GMT, B=BootNo), silent = FALSE)
                if (class(val_GMT) == "try-error") {
                  val_GMT = try(rms::validate(m2_GMT, B=BootNo), silent = FALSE)
                  if (class(val_GMT) == "try-error") {
                    return(NULL)
                  }
                }
                # do not delete
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, val_GMT, file = file.path(tempdir(), "val_GMT.qs"))
                }
              }
              if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "cal_GMT.qs"))){
                cal_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "cal_GMT.qs"))
              } else {
                cal_GMT <- try(rms::calibrate(m2_GMT, B=BootNo), silent = FALSE)
                if (class(cal_GMT) == "try-error") {
                  cal_GMT = try(rms::calibrate(m2_GMT, B=BootNo), silent = FALSE)
                  if (class(cal_GMT) == "try-error") {
                    return(NULL)
                  }
                }
                OUTPUT_DATA$table_prediction_cal_GMT <<- cal_GMT
                if(input$intermediate_file != "No"){
                  QS_SAVE(nthreads = PARALLEL, cal_GMT, file = file.path(tempdir(), "cal_GMT.qs"))
                }
              }
              options(datadist=NULL)
              log_nomogram_GMT <- nomogram(m2_GMT, fun = plogis, lp=F, funlabel=paste("Evidence level", paste(input$evidence_level, collapse = "/")))
              OUTPUT_DATA$table_prediction_log_nomogram_GMT = log_nomogram_GMT
              OUTPUT_DATA$table_prediction_m2_GMT_val_1 <<- val_GMT[9,5]
              OUTPUT_DATA$table_prediction_m2_GMT_val_2 <<- val_GMT[2,5]
              OUTPUT_DATA$table_prediction_m2_GMT_val_3 <<- val_GMT[1,5]
              OUTPUT_DATA$table_prediction_m2_GMT_val_4 <<- val_GMT[3,5]
              OUTPUT_DATA$table_prediction_m2_GMT_val_5 <<- val_GMT[4,5]
              OUTPUT_DATA$table_prediction_m2_GMT_pvalue <<- m2_GMT$stats["P"][[1]]
              m2_GMT$x <- NULL
              m2_GMT$y <- NULL
              m2_GMT$linear.predictors <- NULL
              m2_GMT$fitted.values <- NULL
              m2_GMT$residuals <- NULL
              m2_GMT$weights <- NULL
              m2_GMT$na.action <- NULL
              QS_SAVE(nthreads = PARALLEL, m2_GMT, file = file.path(tempdir(), "m2_GMT.qs"))
              rm(m2_GMT)
              rm(val_GMT)
            }
          }
        }
        incProgress(1 / 13)
        if(length(colnames(Data_forest_tmp_8)) > 1){
          Data_forest_tmp_8$patientid = 1:length(Data_forest$EP_treat)
          dca_thresholds = seq(0, 0.25, 0.01)
          set.seed(1212)
          FLAG = TRUE
          # create a 5-fold cross validation set, 1 repeat which is the base case, change to suit your use case
          if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "glm_param.qs"))){
            df_prediction = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "glm_param.qs"))
          } else {
            cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8, strata = "EP_treat", v = 5, repeats = 1)
            formula_treat = as.formula(paste0("EP_treat ~ `",paste(colnames(Data_forest_tmp_8)[!colnames(Data_forest_tmp_8) %in% c("EP_treat", "patientid")], collapse="` + `"),"`"))
            df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
            if (all(class(df_crossval_predictions) == "try-error")) {
              set.seed(12122)
              cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8, strata = "EP_treat", v = 5, repeats = 1)
              df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
              if (all(class(df_crossval_predictions) == "try-error")) {
                set.seed(12222)
                cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8, strata = "EP_treat", v = 5, repeats = 1)
                df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
                if (all(class(df_crossval_predictions) == "try-error")) {
                  FLAG = FALSE
                }
              }
            }
            incProgress(1 / 13)
            if(FLAG){
              df_prediction_id = c(rsample::assessment(df_crossval_predictions$splits[[1]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[2]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[3]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[4]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[5]])$patientid#,
                                   # rsample::assessment(df_crossval_predictions$splits[[6]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[7]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[8]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[9]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[10]])$patientid
              )
              df_prediction_score = c(predict(df_crossval_predictions$glm_analysis[[1]], rsample::assessment(df_crossval_predictions$splits[[1]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[2]], rsample::assessment(df_crossval_predictions$splits[[2]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[3]], rsample::assessment(df_crossval_predictions$splits[[3]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[4]], rsample::assessment(df_crossval_predictions$splits[[4]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[5]], rsample::assessment(df_crossval_predictions$splits[[5]]),type = "fitted")#,
                                      # predict(df_crossval_predictions$glm_analysis[[6]], rsample::assessment(df_crossval_predictions$splits[[6]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[7]], rsample::assessment(df_crossval_predictions$splits[[7]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[8]], rsample::assessment(df_crossval_predictions$splits[[8]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[9]], rsample::assessment(df_crossval_predictions$splits[[9]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[10]], rsample::assessment(df_crossval_predictions$splits[[10]]),type = "fitted")
              )
              rm(df_crossval_predictions)

              df_prediction = data.frame(df_prediction_id, df_prediction_score)
              colnames(df_prediction) = c("patientid", ".fitted")
              if(input$intermediate_file != "No"){
                QS_SAVE(nthreads = PARALLEL, df_prediction, file=file.path(tempdir(), "glm_param.qs"))
              }

              df_cv_pred <-
                Data_forest_tmp_8 %>%
                dplyr::left_join(
                  df_prediction,
                  by = 'patientid'
                )
              if("PS" %in% colnames(Data_forest_tmp_8)){
                PS_data = data.frame(sort(unique(Data_forest_tmp_8$PS)))
                colnames(PS_data) = c("PS_type")
                save(file = file.path(tempdir(), "PS_data.rda"), PS_data)
              }
              if("Lines" %in% colnames(Data_forest_tmp_8)){
                Lines_data = data.frame(sort(unique(Data_forest_tmp_8$Lines)))
                colnames(Lines_data) = c("Lines_type")
                save(file = file.path(tempdir(), "Lines_data.rda"), Lines_data)
              }
              if(input$machine_learning == "Yes"){
                set.seed(1212)
                Data_ML = data.table(Data_forest_tmp_8 %>%
                                       dplyr::select(-patientid))
                Data_ML$EP_treat = factor(Data_ML$EP_treat)

                Data_ML[, (names(Data_ML)) := lapply(.SD, function(x) if (is.character(x)) as.factor(x) else x)]
                factor_cols <- names(Data_ML)[sapply(Data_ML, is.factor)]
                factor_cols = factor_cols[factor_cols != "EP_treat"]
                # if (length(factor_cols) > 0) {
                #   Data_ML <- one_hot(Data_ML, cols = factor_cols, sparsifyNAs = FALSE, dropCols = TRUE)
                # }
                cls_met <- metric_set(roc_auc)
                Data_cv <- vfold_cv(v=5, Data_ML, strata = EP_treat)
                Data_train = vfold_cv(v=5,
                                      rsample::analysis(Data_cv$splits[[1]]) %>%
                                        dplyr::slice_sample(n = min(nrow(rsample::analysis(Data_cv$splits[[1]])), 20000)),
                                      strata = EP_treat)
                Data_parameter_train = rsample::analysis(Data_train$splits[[1]])
                Data_parameter_test = rsample::assessment(Data_train$splits[[1]])
                Data_parameter_cv = Data_parameter_train |>
                  rsample::vfold_cv(v = 5, strata = EP_treat)
                Data_recipe_train <- Data_parameter_train |>
                  recipe(EP_treat ~ .) |>
                  step_integer(all_nominal_predictors()) |>
                  step_nzv(all_predictors())# |>
                Data_recipe <- Data_ML |>
                  recipe(EP_treat ~ .) |>
                  step_integer(all_nominal_predictors()) |>
                  step_nzv(all_predictors())
                rf_param_space <- parameters(
                  mtry(range = c(1, 10)),
                  min_n(range = c(2, 50))
                )
                lgb_param_space <- parameters(
                  mtry(range = c(1, 10)),
                  tree_depth(range = c(2, 10)),
                  min_n(range = c(2, 50)),
                  num_leaves(range = c(5, 50))
                )
                cv_fit_pred <- function(recipe, spec, df_train, df_test, df_cv = NULL, iter = 500, init_grid_n = 8, param_space = NULL) {
                  if (is.null(df_cv)) {
                    wf <- workflow() |>
                      add_recipe(recipe) |>
                      add_model(spec)
                    params_grid <- NULL
                  } else {
                    cv_wf <- workflow() |>
                      add_recipe(recipe) |>
                      add_model(spec)
                    random_grid <- grid_random(param_space, size = init_grid_n)
                    # 初期グリッド（ランダム or グリッド）
                    initial_grid <- tune::tune_grid(
                      cv_wf,
                      resamples = df_cv,
                      grid = random_grid,
                      metrics = metric_set(roc_auc),
                      control = control_grid(allow_par = DOCKER, parallel_over = "everything",verbose = TRUE,save_pred = TRUE)
                    )

                    # ベイズ最適化
                    params_grid <- tune::tune_bayes(
                      cv_wf,
                      resamples = df_cv,
                      param_info = param_space,
                      initial = initial_grid,
                      iter = iter,
                      metrics = metric_set(roc_auc),
                      control = tune::control_bayes(
                        no_improve = NO_IMPROVE,          # 改善がないと停止
                        parallel_over = "everything",
                        allow_par = DOCKER,
                        verbose = TRUE,
                        verbose_iter = TRUE,
                        save_pred = TRUE
                      )
                    )

                    best_params <- params_grid |> tune::select_best(metric = "roc_auc")

                    wf <- cv_wf |> tune::finalize_workflow(best_params)
                  }

                  wf_fit <- wf |> fit(data = df_train)
                  pred <- wf_fit |> augment(new_data = df_test)

                  return(list(
                    params_grid = params_grid,
                    wf_fit = wf_fit,
                    pred = pred
                  ))
                }
                if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_param.qs"))){
                  lgb_parameter = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_param.qs"))
                } else {
                  lgb_spec <-
                    parsnip::boost_tree(mode = "classification") %>%
                    parsnip::set_args(
                      tree_depth = tune(), min_n = tune(), mtry = tune(), trees = 2000, learn_rate = 0.05, lambda_l1 = 0.9) |>
                    set_engine("lightgbm", num_leaves = tune(),
                               nthread = PARALLEL)
                  lgb_result <-
                    Data_recipe_train |>
                    cv_fit_pred(lgb_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = lgb_param_space)
                  lgb_parameter = (lgb_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                  if(input$intermediate_file != "No"){
                    QS_SAVE(nthreads = PARALLEL, lgb_parameter, file=file.path(tempdir(), "lgb_param.qs"))
                  }
                }
                if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_param.qs"))){
                  rf_parameter = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_param.qs"))
                } else {
                  rf_spec <-
                    parsnip::rand_forest() |>
                    parsnip::set_args(mtry = tune(), min_n = tune(), trees = 1000) |>
                    parsnip::set_engine(
                      'ranger',
                      num.threads = PARALLEL
                    ) |>
                    parsnip::set_mode('classification')
                  rf_result <-
                    Data_recipe_train |>
                    cv_fit_pred(rf_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = rf_param_space)
                  rf_parameter = (rf_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                  if(input$intermediate_file != "No"){
                    QS_SAVE(nthreads = PARALLEL, rf_parameter, file=file.path(tempdir(), "rf_param.qs"))
                  }
                }
                set.seed(1212)
                lgb_model <-
                  parsnip::boost_tree(mode = "classification") %>%
                  set_engine("lightgbm", num_leaves = lgb_parameter$num_leaves,
                             importance = TRUE,  # 変数重要度を有効化
                             lambda_l1 = .9,
                             min_n = lgb_parameter$min_n,
                             tree_depth = lgb_parameter$tree_depth,
                             learn_rate = 0.05,
                             trees = 2000,
                             mtry = lgb_parameter$mtry,
                             nthread = PARALLEL)
                lgb_wflow <-
                  workflow() %>%
                  add_model(lgb_model) %>%
                  add_recipe(Data_recipe)
                if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_m2.qs"))){
                  lgb_m2 = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_m2.qs"))
                } else {
                  lgb_m2 = lgb_wflow %>% fit(Data_ML)
                  lgb_m2 <- butcher:: butcher(lgb_m2)
                  QS_SAVE(nthreads = PARALLEL, lgb_m2, file = file.path(tempdir(), "lgb_m2.qs"))
                }
                rf_model <-
                  parsnip::rand_forest() |>
                  parsnip::set_args(trees = 2000, min_n = rf_parameter$min_n, mtry = rf_parameter$mtry) |>
                  parsnip::set_engine(
                    'ranger', keep.inbag = TRUE, importance = "impurity",
                    num.threads = PARALLEL
                  ) |>
                  parsnip::set_mode('classification')
                rf_wflow <-
                  workflow() %>%
                  add_model(rf_model) %>%
                  add_recipe(Data_recipe)
                if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_m2.qs"))){
                  rf_m2 = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_m2.qs"))
                } else {
                  rf_m2 = rf_wflow %>% fit(Data_ML)
                  rf_m2 <- butcher:: butcher(rf_m2)
                  QS_SAVE(nthreads = PARALLEL, rf_m2, file = file.path(tempdir(), "rf_m2.qs"))
                }
                X <- Data_ML %>% select(-EP_treat)
                X <- if (nrow(X) > 10000) {
                  X[sample(nrow(X), 10000), , drop = FALSE]
                } else {
                  X
                }
                pred_wrapper <- function(object, newdata) {
                  predict(object, new_data = newdata, type = "prob") %>% dplyr::pull(2)
                }
                set.seed(123)
                if(input$importance == "SHAP value (heavy analysis)"){
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_lgb.qs"))){
                    shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_lgb.qs"))
                  } else {
                    shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                     # ncores = PARALLEL,
                                                     object = lgb_m2,
                                                     X = X[seq_len(min(nrow(X), 5000)),],
                                                     pred_wrapper = pred_wrapper,
                                                     nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                    )
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_lgb.qs"))
                    }
                  }
                  rm(lgb_m2)
                  # 各特徴量ごとの平均絶対値を計算
                  mean_abs_shap <- colMeans(abs(shap_values))
                  # データフレームに変換
                  importance_df <- data.frame(
                    feature = names(mean_abs_shap),
                    mean_abs_shap = mean_abs_shap
                  )
                  importance_df = importance_df %>%
                    dplyr::arrange(mean_abs_shap)
                  shap_long <- as.data.frame(shap_values) %>%
                    mutate(id = row_number()) %>%   # サンプルIDを追加
                    pivot_longer(
                      cols = -id,
                      names_to = "feature",
                      values_to = "shap_value"
                    )
                  # ② Xの元データもlong形式に変換（行番号を付与）
                  X_long <- X %>%
                    mutate(id = row_number()) %>%
                    mutate(across(where(is.factor),
                                  ~ case_when(
                                    . %in% c("Yes", "Male") ~ 1,
                                    . %in% c("No", "Female") ~ 0,
                                    TRUE ~ as.numeric(as.character(.))
                                  ))) %>%
                    pivot_longer(
                      cols = -id,
                      names_to = "feature",
                      values_to = "orig_value"
                    ) %>%
                    group_by(feature) %>%
                    group_by(feature) %>%
                    mutate(norm_value = case_when(
                      max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                      TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                    )) %>%
                    ungroup()
                  # ③ shap_longとX_longをidとfeatureで結合
                  shap_long <- shap_long %>%
                    left_join(X_long, by = c("id", "feature"))
                  shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                  importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                  shap_long = shap_long %>%
                    filter(feature %in% importance_df$feature)
                  # 棒グラフでプロット（横向き）
                  gSHAP1_lgb <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                    geom_bar(stat = "identity") +
                    coord_flip() +
                    labs(title = "LightGBM, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                    scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                    theme_minimal() +
                    theme(legend.title = element_blank(),  # 凡例タイトルを削除
                          legend.position = "bottom")  # 凡例を下に配置
                  gSHAP2_lgb = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                    geom_quasirandom(alpha = 0.2, width = 0.2) +
                    scale_color_gradient(low = "blue", high = "red",
                                         limits = c(0, 1),
                                         breaks = c(0, 1),
                                         labels = c("Low/No/Female", "High/Yes/Male")) +
                    theme_minimal() +
                    labs(title = "SHAP Value Distribution", color = "Feature Value")
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_rf.qs"))){
                    shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_rf.qs"))
                  } else {
                    shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                     # ncores = PARALLEL,
                                                     object = rf_m2,
                                                     X = X[seq_len(min(nrow(X), 5000)),],
                                                     pred_wrapper = pred_wrapper,
                                                     nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                    )
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_rf.qs"))
                    }
                  }
                  rm(rf_m2)
                  # 各特徴量のSHAP値の符号（正/負）を判別し、色を付ける
                  shap_values_with_sign <- shap_values
                  shap_values_with_sign[] <- ifelse(shap_values_with_sign >= 0, "Positive", "Negative")
                  # 各特徴量ごとの平均絶対値を計算
                  mean_abs_shap <- colMeans(abs(shap_values))
                  # データフレームに変換
                  importance_df <- data.frame(
                    feature = names(mean_abs_shap),
                    mean_abs_shap = mean_abs_shap,
                    sign = shap_values_with_sign[1,]  # ここでSHAP値の符号を代表として使う
                  )
                  importance_df = importance_df %>%
                    dplyr::arrange(mean_abs_shap)
                  shap_long <- as.data.frame(shap_values) %>%
                    mutate(id = row_number()) %>%   # サンプルIDを追加
                    pivot_longer(
                      cols = -id,
                      names_to = "feature",
                      values_to = "shap_value"
                    )
                  # ② Xの元データもlong形式に変換（行番号を付与）
                  X_long <- X %>%
                    mutate(id = row_number()) %>%
                    mutate(across(where(is.factor),
                                  ~ case_when(
                                    . %in% c("Yes", "Male") ~ 1,
                                    . %in% c("No", "Female") ~ 0,
                                    TRUE ~ as.numeric(as.character(.))
                                  ))) %>%
                    pivot_longer(
                      cols = -id,
                      names_to = "feature",
                      values_to = "orig_value"
                    ) %>%
                    group_by(feature) %>%
                    group_by(feature) %>%
                    mutate(norm_value = case_when(
                      max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                      TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                    )) %>%
                    ungroup()
                  # ③ shap_longとX_longをidとfeatureで結合
                  shap_long <- shap_long %>%
                    left_join(X_long, by = c("id", "feature"))
                  shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                  importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                  shap_long = shap_long %>%
                    filter(feature %in% importance_df$feature)
                  # 棒グラフでプロット（横向き）
                  gSHAP1_rf <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                    geom_bar(stat = "identity") +
                    coord_flip() +
                    labs(title = "RandomForest, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                    scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                    theme_minimal() +
                    theme(legend.title = element_blank(),  # 凡例タイトルを削除
                          legend.position = "bottom")  # 凡例を下に配置
                  gSHAP2_rf = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                    geom_quasirandom(alpha = 0.2, width = 0.2) +
                    scale_color_gradient(low = "blue", high = "red",
                                         limits = c(0, 1),
                                         breaks = c(0, 1),
                                         labels = c("Low/No/Female", "High/Yes/Male")) +
                    theme_minimal() +
                    labs(title = "SHAP Value Distribution", color = "Feature Value")
                } else {
                  importances <- map(Data_cv$splits, function(split) {
                    lgb_fit <- fit(lgb_wflow, data = analysis(split))
                    lgb_booster <- extract_fit_engine(lgb_fit)
                    imp <- lgb.importance(lgb_booster)
                    imp %>% select(Feature, Gain)  # 必要な列だけ
                  })
                  # 2. importance を平均化（split ごとの dataframe を結合 → group_by → summarise）
                  lgb_vi <- bind_rows(importances) %>%
                    group_by(Feature) %>%
                    summarise(importance = mean(Gain, na.rm = TRUE)) %>%
                    arrange(desc(importance))
                  # 3. 可視化
                  gSHAP1_lgb <- ggplot(lgb_vi[seq_len(min(30,nrow(lgb_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                    geom_bar(stat = "identity", fill = "steelblue") +
                    coord_flip() +
                    labs(
                      title = "Top 30 Variable Importance from LightGBM (CV Averaged)",
                      x = "Feature",
                      y = "Importance"
                    ) +
                    theme_minimal()
                  gSHAP2_lgb = gg_empty()
                  rf_importances <- map(Data_cv$splits, function(split) {
                    rf_fit <- fit(rf_wflow, data = analysis(split))
                    rf_model_fitted <- extract_fit_engine(rf_fit)
                    # 変数重要度を tibble 形式で取得
                    tibble(
                      Feature = names(rf_model_fitted$variable.importance),
                      Importance = rf_model_fitted$variable.importance
                    )
                  })
                  # 重要度を平均化
                  rf_vi <- bind_rows(rf_importances) %>%
                    group_by(Feature) %>%
                    summarise(importance = mean(Importance, na.rm = TRUE)) %>%
                    arrange(desc(importance))
                  # 3. 可視化
                  gSHAP1_rf <- ggplot(rf_vi[seq_len(min(30,nrow(rf_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                    geom_bar(stat = "identity", fill = "steelblue") +
                    coord_flip() +
                    labs(
                      title = "Top 30 Variable Importance from Random Forest (CV Averaged)",
                      x = "Feature",
                      y = "Importance"
                    ) +
                    theme_minimal()
                  gSHAP2_rf = gg_empty()
                }
                keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

                set.seed(1212)
                rf_res <- rf_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                rf_res_all = rf_res %>%
                  pull(.predictions) %>%
                  bind_rows() %>%
                  dplyr::arrange(.row)
                # df_cv_pred$rf  = rf_res_all$.pred_1
                ROC_rf <- roc(EP_treat ~ .pred_1, data = rf_res_all, ci = TRUE)
                gROC_rf = ggroc(ROC_rf,
                                size = 1, #サイズ
                                legacy.axes = TRUE) +
                  geom_abline(color = "dark grey", size = 0.5) +
                  theme_classic() +
                  labs(
                    title =
                      paste0("Random forest, AUC: ", format_p(ROC_rf$auc[1], digits = 3), " (95%CI:",format_p(ROC_rf$ci[1], digits = 3), "-", format_p(ROC_rf$ci[3], digits = 3), ")"),
                    subtitle = paste0("mtry/min_n=",rf_parameter$mtry,"/",rf_parameter$min_n)
                  )

                set.seed(1212)
                lgb_res <- lgb_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                lgb_res_all = lgb_res %>%
                  pull(.predictions) %>%
                  bind_rows() %>%
                  dplyr::arrange(.row)
                # df_cv_pred$lgb  = lgb_res_all$.pred_1
                ROC_lgb <- roc(EP_treat ~ .pred_1, data = lgb_res_all, ci = TRUE)
                gROC_lgb = ggroc(ROC_lgb,
                                 size = 1, #サイズ
                                 legacy.axes = TRUE) +
                  geom_abline(color = "dark grey", size = 0.5) +
                  theme_classic() +
                  labs(
                    title =
                      paste0("LightGBM, AUC: ", format_p(ROC_lgb$auc[1], digits = 3), " (95%CI:",format_p(ROC_lgb$ci[1], digits = 3), "-", format_p(ROC_lgb$ci[3], digits = 3), ")"),
                    subtitle = paste0("mtry/min_n/tree_depth/num_leaves=",lgb_parameter$mtry,"/",lgb_parameter$min_n,"/",lgb_parameter$tree_depth,"/",lgb_parameter$num_leaves)
                  )

                g1 = collect_predictions(lgb_res) %>%
                  ggplot(aes(.pred_1)) +
                  geom_histogram(col = "white", bins = 40) +
                  facet_wrap(~ EP_treat, ncol = 1) +
                  geom_rug(col = "blue", alpha = 1 / 2) +
                  labs(x = "LightGBM: Probability estimate of treatment reach")
                g2 = lgb_res %>% cal_plot_windowed(truth = EP_treat, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                  ggtitle("LightGBM: Calibration curve")
                set.seed(1212)
                iso_val <- cal_validate_logistic(lgb_res, metrics = cls_met,
                                                 save_pred = TRUE, times = 25)
                # collect_metrics(iso_val)
                g3= collect_predictions(iso_val) %>%
                  filter(.type == "calibrated") %>%
                  cal_plot_windowed(truth = EP_treat, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                  ggtitle("LightGBM: Logistic calibration via GAM")

                df_cv_pred$lgb  = (iso_val %>%
                                     pull(.predictions_cal) %>%
                                     bind_rows() %>%
                                     dplyr::arrange(.row))$.pred_1
                g4 = collect_predictions(rf_res) %>%
                  ggplot(aes(.pred_1)) +
                  geom_histogram(col = "white", bins = 40) +
                  facet_wrap(~ EP_treat, ncol = 1) +
                  geom_rug(col = "blue", alpha = 1 / 2) +
                  labs(x = "Random forest: Probability estimate of treatment reach")
                g5 = rf_res %>% cal_plot_windowed(truth = EP_treat, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                  ggtitle("Random forest: Calibration curve")

                iso_val <- cal_validate_logistic(rf_res, metrics = cls_met,
                                                 save_pred = TRUE, times = 25)
                collect_metrics(iso_val)
                g6= collect_predictions(iso_val) %>%
                  filter(.type == "calibrated") %>%
                  cal_plot_windowed(truth = EP_treat, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                  ggtitle("Random forest: Logistic calibration via GAM")
                df_cv_pred$rf  = (iso_val %>%
                                    pull(.predictions_cal) %>%
                                    bind_rows() %>%
                                    dplyr::arrange(.row))$.pred_1
                suffix <- ""
                plot_names <- c(
                  "g1", "g2", "g3", "gROC_lgb", "gSHAP1_lgb", "gSHAP2_lgb",
                  "g4", "g5", "g6", "gROC_rf", "gSHAP1_rf", "gSHAP2_rf"
                )
                plot_var_names <- paste0(plot_names, suffix)
                plots_list <- mget(plot_var_names)
                # 5. lapply() を使って、OUTPUT_DATA に動的に保存
                lapply(1:length(plot_names), function(i) {
                  output_key <- paste0("figure_DCA_pre_CGP_Machine_learning_", plot_names[i], suffix)
                  OUTPUT_DATA[[output_key]] <- plots_list[[paste0(plot_names[i], suffix)]]
                })
              }
              OUTPUT_DATA$table_prediction_df_cv_pred = df_cv_pred
            }
          }
        }
        incProgress(1 / 13)
        if(length(colnames(Data_forest_tmp_8_rec)) > 1){
          Data_forest_tmp_8_rec$patientid = 1:length(Data_forest_tmp_8_rec$EP_option)
          dca_thresholds = seq(0, 1, 0.01)
          set.seed(1212)
          FLAG = TRUE
          # create a 5-fold cross validation set, 1 repeat which is the base case, change to suit your use case
          if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "glm_param_rec.qs"))){
            df_prediction_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "glm_param_rec.qs"))
          } else {
            cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_rec, strata = "EP_option", v = 5, repeats = 1)
            formula_treat = as.formula(paste0("EP_option ~ `",paste(colnames(Data_forest_tmp_8_rec)[!colnames(Data_forest_tmp_8_rec) %in% c("EP_option", "patientid")], collapse="` + `"),"`"))
            df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
            if (all(class(df_crossval_predictions) == "try-error")) {
              set.seed(12122)
              cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_rec, strata = "EP_option", v = 5, repeats = 1)
              df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
              if (all(class(df_crossval_predictions) == "try-error")) {
                set.seed(12222)
                cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_rec, strata = "EP_option", v = 5, repeats = 1)
                df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
                if (all(class(df_crossval_predictions) == "try-error")) {
                  FLAG = FALSE
                }
              }}
            incProgress(1 / 13)
            if(FLAG){
              df_prediction_id = c(rsample::assessment(df_crossval_predictions$splits[[1]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[2]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[3]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[4]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[5]])$patientid#,
                                   # rsample::assessment(df_crossval_predictions$splits[[6]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[7]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[8]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[9]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[10]])$patientid
              )
              df_prediction_score = c(predict(df_crossval_predictions$glm_analysis[[1]], rsample::assessment(df_crossval_predictions$splits[[1]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[2]], rsample::assessment(df_crossval_predictions$splits[[2]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[3]], rsample::assessment(df_crossval_predictions$splits[[3]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[4]], rsample::assessment(df_crossval_predictions$splits[[4]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[5]], rsample::assessment(df_crossval_predictions$splits[[5]]),type = "fitted")#,
                                      # predict(df_crossval_predictions$glm_analysis[[6]], rsample::assessment(df_crossval_predictions$splits[[6]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[7]], rsample::assessment(df_crossval_predictions$splits[[7]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[8]], rsample::assessment(df_crossval_predictions$splits[[8]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[9]], rsample::assessment(df_crossval_predictions$splits[[9]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[10]], rsample::assessment(df_crossval_predictions$splits[[10]]),type = "fitted")
              )
              rm(df_crossval_predictions)

              df_prediction_rec = data.frame(df_prediction_id, df_prediction_score)
              colnames(df_prediction_rec) = c("patientid", ".fitted")
              if(input$intermediate_file != "No"){
                QS_SAVE(nthreads = PARALLEL, df_prediction_rec, file=file.path(tempdir(), "glm_param_rec.qs"))
              }

              df_cv_pred_rec <-
                Data_forest_tmp_8_rec %>%
                dplyr::left_join(
                  df_prediction_rec,
                  by = 'patientid'
                )
              if("PS" %in% colnames(Data_forest_tmp_8_rec)){
                PS_data_rec = data.frame(sort(unique(Data_forest_tmp_8_rec$PS)))
                colnames(PS_data_rec) = c("PS_type")
                save(file = file.path(tempdir(), "PS_data_rec.rda"), PS_data_rec)
              }
              if("Lines" %in% colnames(Data_forest_tmp_8_rec)){
                Lines_data_rec = data.frame(sort(unique(Data_forest_tmp_8_rec$Lines)))
                colnames(Lines_data_rec) = c("Lines_type")
                save(file = file.path(tempdir(), "Lines_data_rec.rda"), Lines_data_rec)
              }
              calculation_env <- new.env()
              with(calculation_env, {
                if(input$machine_learning == "Yes"){
                  set.seed(1212)
                  Data_ML = data.table(Data_forest_tmp_8_rec %>%
                                         dplyr::select(-patientid))
                  Data_ML$EP_option = factor(Data_ML$EP_option)

                  Data_ML[, (names(Data_ML)) := lapply(.SD, function(x) if (is.character(x)) as.factor(x) else x)]
                  factor_cols <- names(Data_ML)[sapply(Data_ML, is.factor)]
                  factor_cols = factor_cols[factor_cols != "EP_option"]
                  # if (length(factor_cols) > 0) {
                  #   Data_ML <- one_hot(Data_ML, cols = factor_cols, sparsifyNAs = FALSE, dropCols = TRUE)
                  # }
                  cls_met <- metric_set(roc_auc)
                  Data_cv <- vfold_cv(v=5, Data_ML, strata = EP_option)
                  Data_train = vfold_cv(v=5,
                                        rsample::analysis(Data_cv$splits[[1]]) %>%
                                          dplyr::slice_sample(n = min(nrow(rsample::analysis(Data_cv$splits[[1]])), 20000)),
                                        strata = EP_option)
                  Data_parameter_train = rsample::analysis(Data_train$splits[[1]])
                  Data_parameter_test = rsample::assessment(Data_train$splits[[1]])
                  Data_parameter_cv = Data_parameter_train |>
                    rsample::vfold_cv(v = 5, strata = EP_option)
                  Data_recipe_train <- Data_parameter_train |>
                    recipe(EP_option ~ .) |>
                    step_nzv(all_predictors()) |>
                    step_integer(all_nominal_predictors())
                  Data_recipe <- Data_ML |>
                    recipe(EP_option ~ .) |>
                    step_nzv(all_predictors()) |>
                    step_integer(all_nominal_predictors())
                  rf_param_space <- parameters(
                    mtry(range = c(1, 10)),
                    min_n(range = c(2, 50))
                  )
                  lgb_param_space <- parameters(
                    mtry(range = c(1, 10)),
                    tree_depth(range = c(2, 10)),
                    min_n(range = c(2, 50)),
                    num_leaves(range = c(5, 50))
                  )
                  cv_fit_pred <- function(recipe, spec, df_train, df_test, df_cv = NULL, iter = 500, init_grid_n = 8, param_space = NULL) {
                    if (is.null(df_cv)) {
                      wf <- workflow() |>
                        add_recipe(recipe) |>
                        add_model(spec)
                      params_grid <- NULL
                    } else {
                      cv_wf <- workflow() |>
                        add_recipe(recipe) |>
                        add_model(spec)
                      random_grid <- grid_random(param_space, size = init_grid_n)
                      # 初期グリッド（ランダム or グリッド）
                      initial_grid <- tune::tune_grid(
                        cv_wf,
                        resamples = df_cv,
                        grid = random_grid,
                        metrics = metric_set(roc_auc),
                        control = control_grid(allow_par = DOCKER,parallel_over = "everything",verbose = TRUE,save_pred = TRUE)
                      )

                      # ベイズ最適化
                      params_grid <- tune::tune_bayes(
                        cv_wf,
                        resamples = df_cv,
                        param_info = param_space,
                        initial = initial_grid,
                        iter = iter,
                        metrics = metric_set(roc_auc),
                        control = tune::control_bayes(
                          no_improve = NO_IMPROVE,          # 改善がないと停止
                          parallel_over = "everything",
                          allow_par = DOCKER,
                          verbose = TRUE,
                          verbose_iter = TRUE,
                          save_pred = TRUE
                        )
                      )

                      best_params <- params_grid |> tune::select_best(metric = "roc_auc")

                      wf <- cv_wf |> tune::finalize_workflow(best_params)
                    }

                    wf_fit <- wf |> fit(data = df_train)
                    pred <- wf_fit |> augment(new_data = df_test)

                    return(list(
                      params_grid = params_grid,
                      wf_fit = wf_fit,
                      pred = pred
                    ))
                  }
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_param_rec.qs"))){
                    lgb_parameter_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_param_rec.qs"))
                  } else {
                    lgb_spec <-
                      parsnip::boost_tree(mode = "classification") %>%
                      parsnip::set_args(
                        tree_depth = tune(), min_n = tune(), mtry = tune(), trees = 2000, learn_rate = 0.05, lambda_l1 = 0.9) |>
                      set_engine("lightgbm", num_leaves = tune(),
                                 nthread = PARALLEL)
                    lgb_result <-
                      Data_recipe_train |>
                      cv_fit_pred(lgb_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = lgb_param_space)
                    lgb_parameter_rec = (lgb_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, lgb_parameter_rec, file=file.path(tempdir(), "lgb_param_rec.qs"))
                    }
                  }
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_param_rec.qs"))){
                    rf_parameter_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_param_rec.qs"))
                  } else {
                    rf_spec <-
                      parsnip::rand_forest() |>
                      parsnip::set_args(mtry = tune(), min_n = tune(), trees = 1000) |>
                      parsnip::set_engine(
                        'ranger',
                        num.threads = PARALLEL
                      ) |>
                      parsnip::set_mode('classification')
                    rf_result <-
                      Data_recipe_train |>
                      cv_fit_pred(rf_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = rf_param_space)
                    rf_parameter_rec = (rf_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, rf_parameter_rec, file=file.path(tempdir(), "rf_param_rec.qs"))
                    }
                  }
                  set.seed(1212)
                  lgb_model <-
                    parsnip::boost_tree(mode = "classification") %>%
                    set_engine("lightgbm", num_leaves = lgb_parameter_rec$num_leaves,
                               importance = TRUE,  # 変数重要度を有効化
                               lambda_l1 = .9,
                               min_n = lgb_parameter_rec$min_n,
                               tree_depth = lgb_parameter_rec$tree_depth,
                               learn_rate = 0.05,
                               trees = 2000,
                               mtry = lgb_parameter_rec$mtry,
                               nthread = PARALLEL)
                  lgb_wflow <-
                    workflow() %>%
                    add_model(lgb_model) %>%
                    add_recipe(Data_recipe)
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_m2_rec.qs"))){
                    lgb_m2_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_m2_rec.qs"))
                  } else {
                    lgb_m2_rec = lgb_wflow %>% fit(Data_ML)
                    lgb_m2_rec <- butcher:: butcher(lgb_m2_rec)
                    QS_SAVE(nthreads = PARALLEL, lgb_m2_rec, file = file.path(tempdir(), "lgb_m2_rec.qs"))
                  }
                  rf_model <-
                    parsnip::rand_forest() |>
                    parsnip::set_args(trees = 2000, min_n = rf_parameter_rec$min_n, mtry = rf_parameter_rec$mtry) |>
                    parsnip::set_engine(
                      'ranger', keep.inbag = TRUE, importance = "impurity",
                      num.threads = PARALLEL
                    ) |>
                    parsnip::set_mode('classification')
                  rf_wflow <-
                    workflow() %>%
                    add_model(rf_model) %>%
                    add_recipe(Data_recipe)
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_m2_rec.qs"))){
                    rf_m2_rec = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_m2_rec.qs"))
                  } else {
                    rf_m2_rec = rf_wflow %>% fit(Data_ML)
                    rf_m2_rec <- butcher:: butcher(rf_m2_rec)
                    QS_SAVE(nthreads = PARALLEL, rf_m2_rec, file = file.path(tempdir(), "rf_m2_rec.qs"))
                  }
                  X <- Data_ML %>% select(-EP_option)
                  X <- if (nrow(X) > 10000) {
                    X[sample(nrow(X), 10000), , drop = FALSE]
                  } else {
                    X
                  }
                  pred_wrapper <- function(object, newdata) {
                    predict(object, new_data = newdata, type = "prob") %>% dplyr::pull(2)
                  }
                  set.seed(123)
                  if(input$importance == "SHAP value (heavy analysis)"){
                    if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_lgb_rec.qs"))){
                      shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_lgb_rec.qs"))
                    } else {
                      shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                       # ncores = PARALLEL,
                                                       object = lgb_m2_rec,
                                                       X = X[seq_len(min(nrow(X), 5000)),],
                                                       pred_wrapper = pred_wrapper,
                                                       nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                      )
                      if(input$intermediate_file != "No"){
                        QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_lgb_rec.qs"))
                      }
                    }
                    rm(lgb_m2_rec)
                    # 各特徴量のSHAP値の符号（正/負）を判別し、色を付ける
                    shap_values_with_sign <- shap_values
                    shap_values_with_sign[] <- ifelse(shap_values_with_sign >= 0, "Positive", "Negative")
                    # 各特徴量ごとの平均絶対値を計算
                    mean_abs_shap <- colMeans(abs(shap_values))
                    # データフレームに変換
                    importance_df <- data.frame(
                      feature = names(mean_abs_shap),
                      mean_abs_shap = mean_abs_shap,
                      sign = shap_values_with_sign[1,]  # ここでSHAP値の符号を代表として使う
                    )
                    importance_df = importance_df %>%
                      dplyr::arrange(mean_abs_shap)
                    shap_long <- as.data.frame(shap_values) %>%
                      mutate(id = row_number()) %>%   # サンプルIDを追加
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "shap_value"
                      )
                    # ② Xの元データもlong形式に変換（行番号を付与）
                    X_long <- X %>%
                      mutate(id = row_number()) %>%
                      mutate(across(where(is.factor),
                                    ~ case_when(
                                      . %in% c("Yes", "Male") ~ 1,
                                      . %in% c("No", "Female") ~ 0,
                                      TRUE ~ as.numeric(as.character(.))
                                    ))) %>%
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "orig_value"
                      ) %>%
                      group_by(feature) %>%
                      group_by(feature) %>%
                      mutate(norm_value = case_when(
                        max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                        TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                      )) %>%
                      ungroup()
                    # ③ shap_longとX_longをidとfeatureで結合
                    shap_long <- shap_long %>%
                      left_join(X_long, by = c("id", "feature"))
                    shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                    importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                    shap_long = shap_long %>%
                      filter(feature %in% importance_df$feature)
                    # 棒グラフでプロット（横向き）
                    gSHAP1_lgb_rec <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                      geom_bar(stat = "identity") +
                      coord_flip() +
                      labs(title = "LightGBM, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                      scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                      theme_minimal() +
                      theme(legend.title = element_blank(),  # 凡例タイトルを削除
                            legend.position = "bottom")  # 凡例を下に配置
                    gSHAP2_lgb_rec = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                      geom_quasirandom(alpha = 0.2, width = 0.2) +
                      scale_color_gradient(low = "blue", high = "red",
                                           limits = c(0, 1),
                                           breaks = c(0, 1),
                                           labels = c("Low/No/Female", "High/Yes/Male")) +
                      theme_minimal() +
                      labs(title = "SHAP Value Distribution", color = "Feature Value")
                    if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_rf_rec.qs"))){
                      shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_rf_rec.qs"))
                    } else {
                      shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                       # ncores = PARALLEL,
                                                       object = rf_m2_rec,
                                                       X = X[seq_len(min(nrow(X), 5000)),],
                                                       pred_wrapper = pred_wrapper,
                                                       nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                      )
                      if(input$intermediate_file != "No"){
                        QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_rf_rec.qs"))
                      }
                    }
                    rm(rf_m2_rec)
                    # 各特徴量のSHAP値の符号（正/負）を判別し、色を付ける
                    shap_values_with_sign <- shap_values
                    shap_values_with_sign[] <- ifelse(shap_values_with_sign >= 0, "Positive", "Negative")
                    # 各特徴量ごとの平均絶対値を計算
                    mean_abs_shap <- colMeans(abs(shap_values))
                    # データフレームに変換
                    importance_df <- data.frame(
                      feature = names(mean_abs_shap),
                      mean_abs_shap = mean_abs_shap,
                      sign = shap_values_with_sign[1,]  # ここでSHAP値の符号を代表として使う
                    )
                    importance_df = importance_df %>%
                      dplyr::arrange(mean_abs_shap)
                    shap_long <- as.data.frame(shap_values) %>%
                      mutate(id = row_number()) %>%   # サンプルIDを追加
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "shap_value"
                      )
                    # ② Xの元データもlong形式に変換（行番号を付与）
                    X_long <- X %>%
                      mutate(id = row_number()) %>%
                      mutate(across(where(is.factor),
                                    ~ case_when(
                                      . %in% c("Yes", "Male") ~ 1,
                                      . %in% c("No", "Female") ~ 0,
                                      TRUE ~ as.numeric(as.character(.))
                                    ))) %>%
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "orig_value"
                      ) %>%
                      group_by(feature) %>%
                      group_by(feature) %>%
                      mutate(norm_value = case_when(
                        max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                        TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                      )) %>%
                      ungroup()
                    # ③ shap_longとX_longをidとfeatureで結合
                    shap_long <- shap_long %>%
                      left_join(X_long, by = c("id", "feature"))
                    shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                    importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                    shap_long = shap_long %>%
                      filter(feature %in% importance_df$feature)
                    # 棒グラフでプロット（横向き）
                    gSHAP1_rf_rec <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                      geom_bar(stat = "identity") +
                      coord_flip() +
                      labs(title = "RandomForest, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                      scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                      theme_minimal() +
                      theme(legend.title = element_blank(),  # 凡例タイトルを削除
                            legend.position = "bottom")  # 凡例を下に配置
                    gSHAP2_rf_rec = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                      geom_quasirandom(alpha = 0.2, width = 0.2) +
                      scale_color_gradient(low = "blue", high = "red",
                                           limits = c(0, 1),
                                           breaks = c(0, 1),
                                           labels = c("Low/No/Female", "High/Yes/Male")) +
                      theme_minimal() +
                      labs(title = "SHAP Value Distribution", color = "Feature Value")
                  } else {
                    importances <- map(Data_cv$splits, function(split) {
                      lgb_fit <- fit(lgb_wflow, data = analysis(split))
                      lgb_booster <- extract_fit_engine(lgb_fit)
                      imp <- lgb.importance(lgb_booster)
                      imp %>% select(Feature, Gain)  # 必要な列だけ
                    })
                    # 2. importance を平均化（split ごとの dataframe を結合 → group_by → summarise）
                    lgb_vi <- bind_rows(importances) %>%
                      group_by(Feature) %>%
                      summarise(importance = mean(Gain, na.rm = TRUE)) %>%
                      arrange(desc(importance))
                    # 3. 可視化
                    gSHAP1_lgb_rec <- ggplot(lgb_vi[seq_len(min(30,nrow(lgb_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      labs(
                        title = "Top 30 Variable Importance from LightGBM (CV Averaged)",
                        x = "Feature",
                        y = "Importance"
                      ) +
                      theme_minimal()
                    gSHAP2_lgb_rec = gg_empty()
                    rf_importances <- map(Data_cv$splits, function(split) {
                      rf_fit <- fit(rf_wflow, data = analysis(split))
                      rf_model_fitted <- extract_fit_engine(rf_fit)
                      # 変数重要度を tibble 形式で取得
                      tibble(
                        Feature = names(rf_model_fitted$variable.importance),
                        Importance = rf_model_fitted$variable.importance
                      )
                    })
                    # 重要度を平均化
                    rf_vi <- bind_rows(rf_importances) %>%
                      group_by(Feature) %>%
                      summarise(importance = mean(Importance, na.rm = TRUE)) %>%
                      arrange(desc(importance))
                    # 3. 可視化
                    gSHAP1_rf_rec <- ggplot(rf_vi[seq_len(min(30,nrow(rf_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      labs(
                        title = "Top 30 Variable Importance from Random Forest (CV Averaged)",
                        x = "Feature",
                        y = "Importance"
                      ) +
                      theme_minimal()
                    gSHAP2_rf_rec = gg_empty()
                  }
                  keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

                  set.seed(1212)
                  rf_res <- rf_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                  rf_res_all = rf_res %>%
                    pull(.predictions) %>%
                    bind_rows() %>%
                    dplyr::arrange(.row)
                  ROC_rf_rec <- roc(EP_option ~ .pred_1, data = rf_res_all, ci = TRUE)
                  gROC_rf_rec = ggroc(ROC_rf_rec,
                                      size = 1, #サイズ
                                      legacy.axes = TRUE) +
                    geom_abline(color = "dark grey", size = 0.5) +
                    theme_classic() +
                    labs(
                      title =
                        paste0("Random forest, AUC: ", format_p(ROC_rf_rec$auc[1], digits = 3), " (95%CI:",format_p(ROC_rf_rec$ci[1], digits = 3), "-", format_p(ROC_rf_rec$ci[3], digits = 3), ")"),
                      subtitle = paste0("mtry/min_n=",rf_parameter_rec$mtry,"/",rf_parameter_rec$min_n)
                    )

                  set.seed(1212)
                  lgb_res <- lgb_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                  lgb_res_all = lgb_res %>%
                    pull(.predictions) %>%
                    bind_rows() %>%
                    dplyr::arrange(.row)
                  ROC_lgb_rec <- roc(EP_option ~ .pred_1, data = lgb_res_all, ci = TRUE)
                  gROC_lgb_rec = ggroc(ROC_lgb_rec,
                                       size = 1, #サイズ
                                       legacy.axes = TRUE) +
                    geom_abline(color = "dark grey", size = 0.5) +
                    theme_classic() +
                    labs(
                      title =
                        paste0("LightGBM, AUC: ", format_p(ROC_lgb_rec$auc[1], digits = 3), " (95%CI:",format_p(ROC_lgb_rec$ci[1], digits = 3), "-", format_p(ROC_lgb_rec$ci[3], digits = 3), ")"),
                      subtitle =   paste0("mtry/min_n/tree_depth/num_leaves=",lgb_parameter_rec$mtry,"/",lgb_parameter_rec$min_n,"/",lgb_parameter_rec$tree_depth,"/",lgb_parameter_rec$num_leaves)
                    )

                  g1_rec = collect_predictions(lgb_res) %>%
                    ggplot(aes(.pred_1)) +
                    geom_histogram(col = "white", bins = 40) +
                    facet_wrap(~ EP_option, ncol = 1) +
                    geom_rug(col = "blue", alpha = 1 / 2) +
                    labs(x = "LightGBM: Probability estimate of treatment recommendation")
                  g2_rec = lgb_res %>% cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("LightGBM: Calibration curve")
                  set.seed(1212)
                  iso_val <- cal_validate_logistic(lgb_res, metrics = cls_met,
                                                   save_pred = TRUE, times = 25)
                  # collect_metrics(iso_val)
                  g3_rec = collect_predictions(iso_val) %>%
                    filter(.type == "calibrated") %>%
                    cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("LightGBM: Logistic calibration via GAM")

                  df_cv_pred_rec$lgb  = (iso_val %>%
                                           pull(.predictions_cal) %>%
                                           bind_rows() %>%
                                           dplyr::arrange(.row))$.pred_1
                  g4_rec = collect_predictions(rf_res) %>%
                    ggplot(aes(.pred_1)) +
                    geom_histogram(col = "white", bins = 40) +
                    facet_wrap(~ EP_option, ncol = 1) +
                    geom_rug(col = "blue", alpha = 1 / 2) +
                    labs(x = "Random forest: Probability estimate of treatment recommendation")
                  g5_rec = rf_res %>% cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("Random forest: Calibration curve")

                  iso_val <- cal_validate_logistic(rf_res, metrics = cls_met,
                                                   save_pred = TRUE, times = 25)
                  collect_metrics(iso_val)
                  g6_rec = collect_predictions(iso_val) %>%
                    filter(.type == "calibrated") %>%
                    cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("Random forest: Logistic calibration via GAM")
                  df_cv_pred_rec$rf  = (iso_val %>%
                                          pull(.predictions_cal) %>%
                                          bind_rows() %>%
                                          dplyr::arrange(.row))$.pred_1
                  suffix <- "_rec"
                  plot_names <- c(
                    "g1", "g2", "g3", "gROC_lgb", "gSHAP1_lgb", "gSHAP2_lgb",
                    "g4", "g5", "g6", "gROC_rf", "gSHAP1_rf", "gSHAP2_rf"
                  )
                  plot_var_names <- paste0(plot_names, suffix)
                  plots_list <- mget(plot_var_names)
                  # 5. lapply() を使って、OUTPUT_DATA に動的に保存
                  lapply(1:length(plot_names), function(i) {
                    output_key <- paste0("figure_DCA_pre_CGP_Machine_learning_", plot_names[i], suffix)
                    OUTPUT_DATA[[output_key]] <- plots_list[[paste0(plot_names[i], suffix)]]
                  })
                }
                OUTPUT_DATA$table_prediction_df_cv_pred_rec = df_cv_pred_rec
              })
              rm(calculation_env)
              gc()
            }
          }
        }
        incProgress(1 / 13)
        if(length(colnames(Data_forest_tmp_8_GMT)) > 1){
          Data_forest_tmp_8_GMT$patientid = 1:length(Data_forest_tmp_8_GMT$EP_option)
          dca_thresholds = seq(0, 1, 0.01)
          set.seed(1212)
          FLAG = TRUE
          # create a 5-fold cross validation set, 1 repeat which is the base case, change to suit your use case
          if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "glm_param_GMT.qs"))){
            df_prediction_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "glm_param_GMT.qs"))
          } else {
            cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_GMT, strata = "EP_option", v = 5, repeats = 1)
            df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
            if (all(class(df_crossval_predictions) == "try-error")) {
              set.seed(12122)
              cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_GMT, strata = "EP_option", v = 5, repeats = 1)
              df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
              if (all(class(df_crossval_predictions) == "try-error")) {
                set.seed(12222)
                cross_validation_samples <- rsample::vfold_cv(Data_forest_tmp_8_GMT, strata = "EP_option", v = 5, repeats = 1)
                df_crossval_predictions <- run_crossval(cross_validation_samples, formula_treat, Penal, Toler)
                if (all(class(df_crossval_predictions) == "try-error")) {
                  FLAG = FALSE
                }
              }
            }
            incProgress(1 / 13)
            if(FLAG){
              df_prediction_id = c(rsample::assessment(df_crossval_predictions$splits[[1]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[2]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[3]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[4]])$patientid,
                                   rsample::assessment(df_crossval_predictions$splits[[5]])$patientid#,
                                   # rsample::assessment(df_crossval_predictions$splits[[6]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[7]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[8]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[9]])$patientid,
                                   # rsample::assessment(df_crossval_predictions$splits[[10]])$patientid
              )
              df_prediction_score = c(predict(df_crossval_predictions$glm_analysis[[1]], rsample::assessment(df_crossval_predictions$splits[[1]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[2]], rsample::assessment(df_crossval_predictions$splits[[2]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[3]], rsample::assessment(df_crossval_predictions$splits[[3]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[4]], rsample::assessment(df_crossval_predictions$splits[[4]]),type = "fitted"),
                                      predict(df_crossval_predictions$glm_analysis[[5]], rsample::assessment(df_crossval_predictions$splits[[5]]),type = "fitted")#,
                                      # predict(df_crossval_predictions$glm_analysis[[6]], rsample::assessment(df_crossval_predictions$splits[[6]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[7]], rsample::assessment(df_crossval_predictions$splits[[7]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[8]], rsample::assessment(df_crossval_predictions$splits[[8]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[9]], rsample::assessment(df_crossval_predictions$splits[[9]]),type = "fitted"),
                                      # predict(df_crossval_predictions$glm_analysis[[10]], rsample::assessment(df_crossval_predictions$splits[[10]]),type = "fitted")
              )
              rm(df_crossval_predictions)

              df_prediction_GMT = data.frame(df_prediction_id, df_prediction_score)
              colnames(df_prediction_GMT) = c("patientid", ".fitted")
              if(input$intermediate_file != "No"){
                QS_SAVE(nthreads = PARALLEL, df_prediction_GMT, file=file.path(tempdir(), "glm_param_GMT.qs"))
              }

              df_cv_pred_GMT <-
                Data_forest_tmp_8_GMT %>%
                dplyr::left_join(
                  df_prediction_GMT,
                  by = 'patientid'
                )
              if("PS" %in% colnames(Data_forest_tmp_8_GMT)){
                PS_data_GMT = data.frame(sort(unique(Data_forest_tmp_8_GMT$PS)))
                colnames(PS_data_GMT) = c("PS_type")
                save(file = file.path(tempdir(), "PS_data_GMT.rda"), PS_data_GMT)
              }
              if("Lines" %in% colnames(Data_forest_tmp_8_GMT)){
                Lines_data_GMT = data.frame(sort(unique(Data_forest_tmp_8_GMT$Lines)))
                colnames(Lines_data_GMT) = c("Lines_type")
                save(file = file.path(tempdir(), "Lines_data_GMT.rda"), Lines_data_GMT)
              }
              calculation_env <- new.env()
              with(calculation_env, {
                if(input$machine_learning == "Yes"){
                  set.seed(1212)
                  Data_ML = data.table(Data_forest_tmp_8_GMT %>%
                                         dplyr::select(-patientid))
                  Data_ML$EP_option = factor(Data_ML$EP_option)

                  Data_ML[, (names(Data_ML)) := lapply(.SD, function(x) if (is.character(x)) as.factor(x) else x)]
                  factor_cols <- names(Data_ML)[sapply(Data_ML, is.factor)]
                  factor_cols = factor_cols[factor_cols != "EP_option"]
                  # if (length(factor_cols) > 0) {
                  #   Data_ML <- one_hot(Data_ML, cols = factor_cols, sparsifyNAs = FALSE, dropCols = TRUE)
                  # }
                  cls_met <- metric_set(roc_auc)
                  Data_cv <- vfold_cv(v=5, Data_ML, strata = EP_option)
                  Data_train = vfold_cv(v=5,
                                        rsample::analysis(Data_cv$splits[[1]]) %>%
                                          dplyr::slice_sample(n = min(nrow(rsample::analysis(Data_cv$splits[[1]])), 20000)),
                                        strata = EP_option)
                  Data_parameter_train = rsample::analysis(Data_train$splits[[1]])
                  Data_parameter_test = rsample::assessment(Data_train$splits[[1]])
                  Data_parameter_cv = Data_parameter_train |>
                    rsample::vfold_cv(v = 5, strata = EP_option)
                  Data_recipe_train <- Data_parameter_train |>
                    recipe(EP_option ~ .) |>
                    step_nzv(all_predictors()) |>
                    step_integer(all_nominal_predictors())
                  Data_recipe <- Data_ML |>
                    recipe(EP_option ~ .) |>
                    step_nzv(all_predictors()) |>
                    step_integer(all_nominal_predictors())
                  rf_param_space <- parameters(
                    mtry(range = c(1, 10)),
                    min_n(range = c(2, 50))
                  )
                  lgb_param_space <- parameters(
                    mtry(range = c(1, 10)),
                    tree_depth(range = c(2, 10)),
                    min_n(range = c(2, 50)),
                    num_leaves(range = c(5, 50))
                  )
                  cv_fit_pred <- function(recipe, spec, df_train, df_test, df_cv = NULL, iter = 500, init_grid_n = 8, param_space = NULL) {
                    if (is.null(df_cv)) {
                      wf <- workflow() |>
                        add_recipe(recipe) |>
                        add_model(spec)
                      params_grid <- NULL
                    } else {
                      cv_wf <- workflow() |>
                        add_recipe(recipe) |>
                        add_model(spec)
                      random_grid <- grid_random(param_space, size = init_grid_n)
                      # 初期グリッド（ランダム or グリッド）
                      initial_grid <- tune::tune_grid(
                        cv_wf,
                        resamples = df_cv,
                        grid = random_grid,
                        metrics = metric_set(roc_auc),
                        control = control_grid(allow_par = DOCKER, parallel_over = "everything", verbose = TRUE, save_pred = TRUE)
                      )

                      # ベイズ最適化
                      params_grid <- tune::tune_bayes(
                        cv_wf,
                        resamples = df_cv,
                        param_info = param_space,
                        initial = initial_grid,
                        iter = iter,
                        metrics = metric_set(roc_auc),
                        control = tune::control_bayes(
                          no_improve = NO_IMPROVE,          # 改善がないと停止
                          parallel_over = "everything",
                          allow_par = DOCKER,
                          verbose = TRUE,
                          verbose_iter = TRUE,
                          save_pred = TRUE
                        )
                      )

                      best_params <- params_grid |> tune::select_best(metric = "roc_auc")

                      wf <- cv_wf |> tune::finalize_workflow(best_params)
                    }

                    wf_fit <- wf |> fit(data = df_train)
                    pred <- wf_fit |> augment(new_data = df_test)

                    return(list(
                      params_grid = params_grid,
                      wf_fit = wf_fit,
                      pred = pred
                    ))
                  }
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_param_GMT.qs"))){
                    lgb_parameter_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_param_GMT.qs"))
                  } else {
                    lgb_spec <-
                      parsnip::boost_tree(mode = "classification") %>%
                      parsnip::set_args(
                        tree_depth = tune(), min_n = tune(), mtry = tune(), trees = 2000, learn_rate = 0.05, lambda_l1 = 0.9) |>
                      set_engine("lightgbm", num_leaves = tune(),
                                 nthread = PARALLEL)
                    lgb_result <-
                      Data_recipe_train |>
                      cv_fit_pred(lgb_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = lgb_param_space)
                    lgb_parameter_GMT = (lgb_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, lgb_parameter_GMT, file=file.path(tempdir(), "lgb_param_GMT.qs"))
                    }
                  }
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_param_GMT.qs"))){
                    rf_parameter_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_param_GMT.qs"))
                  } else {
                    rf_spec <-
                      parsnip::rand_forest() |>
                      parsnip::set_args(mtry = tune(), min_n = tune(), trees = 1000) |>
                      parsnip::set_engine(
                        'ranger',
                        num.threads = PARALLEL
                      ) |>
                      parsnip::set_mode('classification')
                    rf_result <-
                      Data_recipe_train |>
                      cv_fit_pred(rf_spec, Data_parameter_train, Data_parameter_test, Data_parameter_cv, param_space = rf_param_space)
                    rf_parameter_GMT = (rf_result$params_grid %>% show_best(metric = "roc_auc"))[1,]
                    if(input$intermediate_file != "No"){
                      QS_SAVE(nthreads = PARALLEL, rf_parameter_GMT, file=file.path(tempdir(), "rf_param_GMT.qs"))
                    }
                  }
                  set.seed(1212)
                  lgb_model <-
                    parsnip::boost_tree(mode = "classification") %>%
                    set_engine("lightgbm", num_leaves = lgb_parameter_GMT$num_leaves,
                               importance = TRUE,  # 変数重要度を有効化
                               lambda_l1 = .9,
                               min_n = lgb_parameter_GMT$min_n,
                               tree_depth = lgb_parameter_GMT$tree_depth,
                               learn_rate = 0.05,
                               trees = 2000,
                               mtry = lgb_parameter_GMT$mtry,
                               nthread = PARALLEL)
                  lgb_wflow <-
                    workflow() %>%
                    add_model(lgb_model) %>%
                    add_recipe(Data_recipe)
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "lgb_m2_GMT.qs"))){
                    lgb_m2_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "lgb_m2_GMT.qs"))
                  } else {
                    lgb_m2_GMT = lgb_wflow %>% fit(Data_ML)
                    lgb_m2_GMT <- butcher:: butcher(lgb_m2_GMT)
                    QS_SAVE(nthreads = PARALLEL, lgb_m2_GMT, file = file.path(tempdir(), "lgb_m2_GMT.qs"))
                  }
                  rf_model <-
                    parsnip::rand_forest() |>
                    parsnip::set_args(trees = 2000, min_n = rf_parameter_GMT$min_n, mtry = rf_parameter_GMT$mtry) |>
                    parsnip::set_engine(
                      'ranger', keep.inbag = TRUE, importance = "impurity",
                      num.threads = PARALLEL
                    ) |>
                    parsnip::set_mode('classification')
                  rf_wflow <-
                    workflow() %>%
                    add_model(rf_model) %>%
                    add_recipe(Data_recipe)
                  if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "rf_m2_GMT.qs"))){
                    rf_m2_GMT = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "rf_m2_GMT.qs"))
                  } else {
                    rf_m2_GMT = rf_wflow %>% fit(Data_ML)
                    rf_m2_GMT <- butcher:: butcher(rf_m2_GMT)
                    QS_SAVE(nthreads = PARALLEL, rf_m2_GMT, file = file.path(tempdir(), "rf_m2_GMT.qs"))
                  }
                  X <- Data_ML %>% select(-EP_option)
                  X <- if (nrow(X) > 10000) {
                    X[sample(nrow(X), 10000), , drop = FALSE]
                  } else {
                    X
                  }
                  pred_wrapper <- function(object, newdata) {
                    predict(object, new_data = newdata, type = "prob") %>% dplyr::pull(2)
                  }
                  set.seed(123)
                  if(input$importance == "SHAP value (heavy analysis)"){
                    if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_lgb_GMT.qs"))){
                      shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_lgb_GMT.qs"))
                    } else {
                      shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                       # ncores = PARALLEL,
                                                       object = lgb_m2_GMT,
                                                       X = X[seq_len(min(nrow(X), 5000)),],
                                                       pred_wrapper = pred_wrapper,
                                                       nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                      )
                      if(input$intermediate_file != "No"){
                        QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_lgb_GMT.qs"))
                      }
                    }
                    rm(lgb_m2_GMT)
                    # 各特徴量のSHAP値の符号（正/負）を判別し、色を付ける
                    shap_values_with_sign <- shap_values
                    shap_values_with_sign[] <- ifelse(shap_values_with_sign >= 0, "Positive", "Negative")
                    # 各特徴量ごとの平均絶対値を計算
                    mean_abs_shap <- colMeans(abs(shap_values))
                    # データフレームに変換
                    importance_df <- data.frame(
                      feature = names(mean_abs_shap),
                      mean_abs_shap = mean_abs_shap,
                      sign = shap_values_with_sign[1,]  # ここでSHAP値の符号を代表として使う
                    )
                    importance_df = importance_df %>%
                      dplyr::arrange(mean_abs_shap)
                    shap_long <- as.data.frame(shap_values) %>%
                      mutate(id = row_number()) %>%   # サンプルIDを追加
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "shap_value"
                      )
                    # ② Xの元データもlong形式に変換（行番号を付与）
                    X_long <- X %>%
                      mutate(id = row_number()) %>%
                      mutate(across(where(is.factor),
                                    ~ case_when(
                                      . %in% c("Yes", "Male") ~ 1,
                                      . %in% c("No", "Female") ~ 0,
                                      TRUE ~ as.numeric(as.character(.))
                                    ))) %>%
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "orig_value"
                      ) %>%
                      group_by(feature) %>%
                      group_by(feature) %>%
                      mutate(norm_value = case_when(
                        max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                        TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                      )) %>%
                      ungroup()
                    # ③ shap_longとX_longをidとfeatureで結合
                    shap_long <- shap_long %>%
                      left_join(X_long, by = c("id", "feature"))
                    shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                    importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                    shap_long = shap_long %>%
                      filter(feature %in% importance_df$feature)
                    # 棒グラフでプロット（横向き）
                    gSHAP1_lgb_GMT <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                      geom_bar(stat = "identity") +
                      coord_flip() +
                      labs(title = "LightGBM, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                      scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                      theme_minimal() +
                      theme(legend.title = element_blank(),  # 凡例タイトルを削除
                            legend.position = "bottom")  # 凡例を下に配置
                    gSHAP2_lgb_GMT = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                      geom_quasirandom(alpha = 0.2, width = 0.2) +
                      scale_color_gradient(low = "blue", high = "red",
                                           limits = c(0, 1),
                                           breaks = c(0, 1),
                                           labels = c("Low/No/Female", "High/Yes/Male")) +
                      theme_minimal() +
                      labs(title = "SHAP Value Distribution", color = "Feature Value")
                    if(input$new_analysis =="No, use the previous dataset" & input$intermediate_file != "No" & file.exists(file.path(tempdir(), "shap_values_rf_GMT.qs"))){
                      shap_values = QS_READ(nthreads = PARALLEL, file=file.path(tempdir(), "shap_values_rf_GMT.qs"))
                    } else {
                      shap_values <- fastshap::explain(parallel = ifelse(nrow(Data_ML) < 10000, TRUE, FALSE),
                                                       # ncores = PARALLEL,
                                                       object = rf_m2_GMT,
                                                       X = X[seq_len(min(nrow(X), 5000)),],
                                                       pred_wrapper = pred_wrapper,
                                                       nsim = 3  # モンテカルロシミュレーション回数（計算精度と時間のバランスで調整）
                      )
                      if(input$intermediate_file != "No"){
                        QS_SAVE(nthreads = PARALLEL, shap_values, file = file.path(tempdir(), "shap_values_rf_GMT.qs"))
                      }
                    }
                    rm(rf_m2_GMT)
                    # 各特徴量のSHAP値の符号（正/負）を判別し、色を付ける
                    shap_values_with_sign <- shap_values
                    shap_values_with_sign[] <- ifelse(shap_values_with_sign >= 0, "Positive", "Negative")
                    # 各特徴量ごとの平均絶対値を計算
                    mean_abs_shap <- colMeans(abs(shap_values))
                    # データフレームに変換
                    importance_df <- data.frame(
                      feature = names(mean_abs_shap),
                      mean_abs_shap = mean_abs_shap,
                      sign = shap_values_with_sign[1,]  # ここでSHAP値の符号を代表として使う
                    )
                    importance_df = importance_df %>%
                      dplyr::arrange(mean_abs_shap)
                    shap_long <- as.data.frame(shap_values) %>%
                      mutate(id = row_number()) %>%   # サンプルIDを追加
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "shap_value"
                      )
                    # ② Xの元データもlong形式に変換（行番号を付与）
                    X_long <- X %>%
                      mutate(id = row_number()) %>%
                      mutate(across(where(is.factor),
                                    ~ case_when(
                                      . %in% c("Yes", "Male") ~ 1,
                                      . %in% c("No", "Female") ~ 0,
                                      TRUE ~ as.numeric(as.character(.))
                                    ))) %>%
                      pivot_longer(
                        cols = -id,
                        names_to = "feature",
                        values_to = "orig_value"
                      ) %>%
                      group_by(feature) %>%
                      group_by(feature) %>%
                      mutate(norm_value = case_when(
                        max(orig_value) == min(orig_value) ~ orig_value,  # すべて同じならそのまま
                        TRUE ~ (orig_value - min(orig_value)) / (max(orig_value) - min(orig_value))
                      )) %>%
                      ungroup()
                    # ③ shap_longとX_longをidとfeatureで結合
                    shap_long <- shap_long %>%
                      left_join(X_long, by = c("id", "feature"))
                    shap_long$feature = factor(shap_long$feature, levels = importance_df$feature)
                    importance_df = importance_df[max(1, nrow(importance_df)-29):nrow(importance_df),]
                    shap_long = shap_long %>%
                      filter(feature %in% importance_df$feature)
                    # 棒グラフでプロット（横向き）
                    gSHAP1_rf_GMT <- ggplot(importance_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
                      geom_bar(stat = "identity") +
                      coord_flip() +
                      labs(title = "RandomForest, Mean |SHAP value| of Features", x = "Features", y = "Mean |SHAP value|") +
                      scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +  # 色の指定
                      theme_minimal() +
                      theme(legend.title = element_blank(),  # 凡例タイトルを削除
                            legend.position = "bottom")  # 凡例を下に配置
                    gSHAP2_rf_GMT = ggplot(shap_long, aes(x = shap_value, y = feature, color = norm_value)) +
                      geom_quasirandom(alpha = 0.2, width = 0.2) +
                      scale_color_gradient(low = "blue", high = "red",
                                           limits = c(0, 1),
                                           breaks = c(0, 1),
                                           labels = c("Low/No/Female", "High/Yes/Male")) +
                      theme_minimal() +
                      labs(title = "SHAP Value Distribution", color = "Feature Value")
                  } else {
                    importances <- map(Data_cv$splits, function(split) {
                      lgb_fit <- fit(lgb_wflow, data = analysis(split))
                      lgb_booster <- extract_fit_engine(lgb_fit)
                      imp <- lgb.importance(lgb_booster)
                      imp %>% select(Feature, Gain)  # 必要な列だけ
                    })
                    # 2. importance を平均化（split ごとの dataframe を結合 → group_by → summarise）
                    lgb_vi <- bind_rows(importances) %>%
                      group_by(Feature) %>%
                      summarise(importance = mean(Gain, na.rm = TRUE)) %>%
                      arrange(desc(importance))
                    # 3. 可視化
                    gSHAP1_lgb_GMT <- ggplot(lgb_vi[seq_len(min(30,nrow(lgb_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      labs(
                        title = "Top 30 Variable Importance from LightGBM (CV Averaged)",
                        x = "Feature",
                        y = "Importance"
                      ) +
                      theme_minimal()
                    gSHAP2_lgb_GMT = gg_empty()
                    rf_importances <- map(Data_cv$splits, function(split) {
                      rf_fit <- fit(rf_wflow, data = analysis(split))
                      rf_model_fitted <- extract_fit_engine(rf_fit)
                      # 変数重要度を tibble 形式で取得
                      tibble(
                        Feature = names(rf_model_fitted$variable.importance),
                        Importance = rf_model_fitted$variable.importance
                      )
                    })
                    # 重要度を平均化
                    rf_vi <- bind_rows(rf_importances) %>%
                      group_by(Feature) %>%
                      summarise(importance = mean(Importance, na.rm = TRUE)) %>%
                      arrange(desc(importance))
                    # 3. 可視化
                    gSHAP1_rf_GMT <- ggplot(rf_vi[seq_len(min(30,nrow(rf_vi))),], aes(x = reorder(Feature, importance), y = importance)) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      labs(
                        title = "Top 30 Variable Importance from Random Forest (CV Averaged)",
                        x = "Feature",
                        y = "Importance"
                      ) +
                      theme_minimal()
                    gSHAP2_rf_GMT = gg_empty()
                  }
                  keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

                  set.seed(1212)
                  rf_res <- rf_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                  rf_res_all = rf_res %>%
                    pull(.predictions) %>%
                    bind_rows() %>%
                    dplyr::arrange(.row)
                  # df_cv_pred_GMT$rf  = rf_res_all$.pred_1
                  ROC_rf_GMT <- roc(EP_option ~ .pred_1, data = rf_res_all, ci = TRUE)
                  gROC_rf_GMT = ggroc(ROC_rf_GMT,
                                      size = 1, #サイズ
                                      legacy.axes = TRUE) +
                    geom_abline(color = "dark grey", size = 0.5) +
                    theme_classic() +
                    labs(
                      title =
                        paste0("Random forest, AUC: ", format_p(ROC_rf_GMT$auc[1], digits = 3), " (95%CI:",format_p(ROC_rf_GMT$ci[1], digits = 3), "-", format_p(ROC_rf_GMT$ci[3], digits = 3), ")"),
                      subtitle = paste0("mtry/min_n=",rf_parameter_GMT$mtry,"/",rf_parameter_GMT$min_n)
                    )

                  set.seed(1212)
                  lgb_res <- lgb_wflow %>% fit_resamples(resamples = Data_cv, metrics = cls_met, control = keep_pred)
                  lgb_res_all = lgb_res %>%
                    pull(.predictions) %>%
                    bind_rows() %>%
                    dplyr::arrange(.row)
                  # df_cv_pred_rec$lgb  = lgb_res_all$.pred_1
                  ROC_lgb_GMT <- roc(EP_option ~ .pred_1, data = lgb_res_all, ci = TRUE)
                  gROC_lgb_GMT = ggroc(ROC_lgb_GMT,
                                       size = 1, #サイズ
                                       legacy.axes = TRUE) +
                    geom_abline(color = "dark grey", size = 0.5) +
                    theme_classic() +
                    labs(
                      title =
                        paste0("LightGBM, AUC: ", format_p(ROC_lgb_GMT$auc[1], digits = 3), " (95%CI:",format_p(ROC_lgb_GMT$ci[1], digits = 3), "-", format_p(ROC_lgb_GMT$ci[3], digits = 3), ")"),
                      subtitle = paste0("mtry/min_n/tree_depth/num_leaves=",lgb_parameter_GMT$mtry,"/",lgb_parameter_GMT$min_n,"/",lgb_parameter_GMT$tree_depth,"/",lgb_parameter_GMT$num_leaves)
                    )

                  g1_GMT = collect_predictions(lgb_res) %>%
                    ggplot(aes(.pred_1)) +
                    geom_histogram(col = "white", bins = 40) +
                    facet_wrap(~ EP_option, ncol = 1) +
                    geom_rug(col = "blue", alpha = 1 / 2) +
                    labs(x = paste("LightGBM: Probability estimate of evidence level", paste(input$evidence_level, collapse = "/")))
                  g2_GMT = lgb_res %>% cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("LightGBM: Calibration curve")
                  set.seed(1212)
                  iso_val <- cal_validate_logistic(lgb_res, metrics = cls_met,
                                                   save_pred = TRUE, times = 25)
                  # collect_metrics(iso_val)
                  g3_GMT = collect_predictions(iso_val) %>%
                    filter(.type == "calibrated") %>%
                    cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("LightGBM: Logistic calibration via GAM")

                  df_cv_pred_GMT$lgb  = (iso_val %>%
                                           pull(.predictions_cal) %>%
                                           bind_rows() %>%
                                           dplyr::arrange(.row))$.pred_1
                  g4_GMT = collect_predictions(rf_res) %>%
                    ggplot(aes(.pred_1)) +
                    geom_histogram(col = "white", bins = 40) +
                    facet_wrap(~ EP_option, ncol = 1) +
                    geom_rug(col = "blue", alpha = 1 / 2) +
                    labs(x = paste("Random forest: Probability estimate of evidence level", paste(input$evidence_level, collapse = "/")))
                  g5_GMT = rf_res %>% cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("Random forest: Calibration curve")

                  iso_val <- cal_validate_logistic(rf_res, metrics = cls_met,
                                                   save_pred = TRUE, times = 25)
                  collect_metrics(iso_val)
                  g6_GMT = collect_predictions(iso_val) %>%
                    filter(.type == "calibrated") %>%
                    cal_plot_windowed(truth = EP_option, event_level = 'second', estimate = .pred_1, step_size = 0.025) +
                    ggtitle("Random forest: Logistic calibration via GAM")
                  df_cv_pred_GMT$rf  = (iso_val %>%
                                          pull(.predictions_cal) %>%
                                          bind_rows() %>%
                                          dplyr::arrange(.row))$.pred_1
                  suffix <- "_GMT"
                  plot_names <- c(
                    "g1", "g2", "g3", "gROC_lgb", "gSHAP1_lgb", "gSHAP2_lgb",
                    "g4", "g5", "g6", "gROC_rf", "gSHAP1_rf", "gSHAP2_rf"
                  )
                  plot_var_names <- paste0(plot_names, suffix)
                  plots_list <- mget(plot_var_names)
                  # 5. lapply() を使って、OUTPUT_DATA に動的に保存
                  lapply(1:length(plot_names), function(i) {
                    output_key <- paste0("figure_DCA_pre_CGP_Machine_learning_", plot_names[i], suffix)
                    OUTPUT_DATA[[output_key]] <- plots_list[[paste0(plot_names[i], suffix)]]
                  })
                }
                OUTPUT_DATA$table_prediction_df_cv_pred_GMT = df_cv_pred_GMT
              })
              rm(calculation_env)
              gc()
            }
          }
        }
      }
      incProgress(1 / 13)
    })
  })
  rm(analysis_env)
  gc()
}



# データを切り替える reactive
selected_df <- reactive({
  req(input$nomogram_type_4)
  if (input$nomogram_type_4 == "EP_option") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_rec)
    return(OUTPUT_DATA$table_prediction_df_cv_pred_rec)
  } else if (input$nomogram_type_4 == "_GMT") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_GMT)
    return(OUTPUT_DATA$table_prediction_df_cv_pred_GMT)
  } else {
    req(OUTPUT_DATA$table_prediction_df_cv_pred)
    return(OUTPUT_DATA$table_prediction_df_cv_pred)
  }
})

output$table_prediction <- DT::renderDataTable(server = FALSE,{
  req(selected_df())
  create_datatable_with_confirm(selected_df())
})

output$figure_DCA_pre_CGP_Machine_learning = renderPlot({
  req(input$nomogram_type_4, input$machine_learning)
  req(input$machine_learning == "Yes")
  suffix <- if (input$nomogram_type_4 == "EP_option") {
    "_rec"
  } else if (input$nomogram_type_4 == "_GMT") {
    "_GMT"
  } else {
    "" # "EP_treat" の場合
  }
  plot_names <- c(
    "g1", "g2", "g3", "gROC_lgb", "gSHAP1_lgb", "gSHAP2_lgb",
    "g4", "g5", "g6", "gROC_rf", "gSHAP1_rf", "gSHAP2_rf"
  )
  full_plot_names <- paste0("figure_DCA_pre_CGP_Machine_learning_", plot_names, suffix)
  plot_list <- lapply(full_plot_names, function(name) {
    OUTPUT_DATA[[name]]
  })
  req(all(!sapply(plot_list, is.null)))
  do.call(grid.arrange, c(plot_list, nrow = 4, ncol = 3))
})

output$figure_DCA_pre_CGP_1 = renderPlot({
  # 必要な入力が揃っているか確認
  req(input$nomogram_type_3, input$machine_learning)
  # 1. nomogram_type_3 に応じて data, formula, xlim を動的に設定
  if (input$nomogram_type_3 == "EP_treat") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred
    formula_target <- EP_treat ~ .fitted
    xlim_value <- c(0, 0.25)
    dca_thresholds = seq(0, 0.25, 0.01)
  } else if (input$nomogram_type_3 == "EP_option") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_rec)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_rec
    formula_target <- EP_option ~ .fitted
    xlim_value <- c(0, 1)
    dca_thresholds = seq(0, 1, 0.01)
  } else if (input$nomogram_type_3 == "_GMT") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_GMT)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_GMT
    formula_target <- EP_option ~ .fitted
    xlim_value <- c(0, 1)
    dca_thresholds = seq(0, 1, 0.01)
  }
  if (input$machine_learning == "Yes") {
    formula_full <- update(formula_target, . ~ . + rf + lgb)
    labels_list <- list(
      .fitted = "5-fold cross-validated logistic regression model",
      rf = "5-fold cross-validated random forest model",
      lgb = "5-fold cross-validated lightGBM model"
    )
  } else {
    formula_full <- formula_target
    labels_list <- list(
      .fitted = "5-fold cross-validated logistic regression model"
    )
  }
  dcurves::dca(
    data = data_to_use,
    formula = formula_full,
    thresholds = dca_thresholds,
    label = labels_list
  )$dca |>
    ggplot(aes(x = threshold, y = net_benefit, color = label)) +
    stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2, fullrange = TRUE) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Threshold probability", y = "Net Benefit", color = "") +
    coord_cartesian(xlim = xlim_value, ylim = c(-0.01, NA)) +
    theme_bw()
})

output$figure_DCA_pre_CGP_ROC = renderPlot({
  # 必要な入力が揃っているか確認
  req(input$nomogram_type_4, input$machine_learning)

  # 1. nomogram_type_4 に応じてデータ、目的変数、タイトルを動的に設定
  if (input$nomogram_type_4 == "EP_treat") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred)
    df_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred
    target_var <- "EP_treat"
    x_label_nomogram <- "Logistic regression: Probability estimate of treatment reach"
    roc_title_base <- "ROC curves for treatment reach"
  } else if (input$nomogram_type_4 == "EP_option") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_rec)
    df_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_rec
    target_var <- "EP_option"
    x_label_nomogram <- "Logistic regression: Probability estimate of treatment recommendation"
    roc_title_base <- "ROC curves for treatment recommendation"
  } else if (input$nomogram_type_4 == "_GMT") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_GMT)
    df_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_GMT
    target_var <- "EP_option"
    x_label_nomogram <- paste("Logistic regression: Probability estimate of evidence level", paste(input$evidence_level, collapse = "/"))
    roc_title_base <- paste("ROC curves for evidence level", paste(input$evidence_level, collapse = "/"))
  }
  gNomogram <- df_to_use %>%
    ggplot(aes(.fitted)) +
    theme_classic() +
    geom_histogram(col = "white", bins = 40) +
    facet_wrap(as.formula(paste0("~", target_var)), ncol = 1) +
    geom_rug(col = "blue", alpha = 1 / 2) +
    labs(x = x_label_nomogram)
  ROC_base <- roc(as.formula(paste0(target_var, " ~ .fitted")), data = df_to_use, ci = TRUE)
  gROC2 <- ggroc(ROC_base, size = 1, legacy.axes = TRUE) +
    geom_abline(color = "dark grey", size = 0.5) +
    theme_classic() +
    labs(
      title = paste0("Logistic regression, AUC: ", format_p(ROC_base$auc[1], digits = 3),
                     " (95%CI:", format_p(ROC_base$ci[1], digits = 3), "-", format_p(ROC_base$ci[3], digits = 3), ")")
    )
  if (input$machine_learning == "Yes") {
    formula_full <- as.formula(paste0(target_var, " ~ .fitted + lgb + rf"))
    labels_list <- c(.fitted = "Logistic regression", lgb = "Light GBM", rf = "Random forest")
  } else {
    formula_full <- as.formula(paste0(target_var, " ~ .fitted"))
    labels_list <- c(.fitted = "Logistic regression")
  }
  ROC_full <- roc(formula_full, data = df_to_use, ci = TRUE)
  gROC3 <- ggroc(ROC_full, size = 1, legacy.axes = TRUE) +
    geom_abline(color = "dark grey", size = 0.5) +
    theme_classic() +
    scale_color_hue(name = "Prediction methods", labels = labels_list) +
    labs(title = roc_title_base)
  grid.arrange(gNomogram, gROC2, gROC3, nrow = 3)
})

output$table_DCA_pre_CGP_2 = render_gt({
  req(input$nomogram_type_3, input$machine_learning)
  if (input$nomogram_type_3 == "EP_treat") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred
    formula_target <- EP_treat ~ .fitted
    dca_thresholds = seq(0, 0.25, 0.01)
  } else if (input$nomogram_type_3 == "EP_option") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_rec)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_rec
    formula_target <- EP_option ~ .fitted
    dca_thresholds = seq(0, 1, 0.01)
  } else if (input$nomogram_type_3 == "_GMT") {
    req(OUTPUT_DATA$table_prediction_df_cv_pred_GMT)
    data_to_use <- OUTPUT_DATA$table_prediction_df_cv_pred_GMT
    formula_target <- EP_option ~ .fitted
    dca_thresholds = seq(0, 1, 0.01)
  }
  if (input$machine_learning == "Yes") {
    formula_full <- update(formula_target, . ~ . + rf + lgb)
    labels_list <- list(
      .fitted = "5-fold cross-validated logistic regression model",
      rf = "5-fold cross-validated random forest model",
      lgb = "5-fold cross-validated lightGBM model"
    )
  } else {
    formula_full <- formula_target
    labels_list <- list(
      .fitted = "5-fold cross-validated logistic regression model"
    )
  }
  dcurves::dca(
    data = data_to_use,
    formula = formula_full,
    thresholds = dca_thresholds,
    label = labels_list
  ) %>%
    net_intervention_avoided() %>%
    as_tibble() %>%
    gt::gt() %>%
    gt::fmt_percent(columns = threshold, decimals = 0) %>%
    gt::fmt(columns = net_benefit, fns = function(x) style_sigfig(x, digits = 3)) %>%
    gt::cols_label(
      label = "Strategy",
      threshold = "Decision Threshold",
      net_benefit = "Net Benefit"
    ) %>%
    gt::cols_align("left", columns = label)
})

output$histology_nomogram_pre_CGP = renderText({
  # 必要な入力が揃っているか確認
  req(input$nomogram_type)

  # nomogram_type に応じて動的にノモグラムオブジェクトを選択
  if (input$nomogram_type == "EP_treat") {
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_8)) > 1)
    req(length(OUTPUT_DATA$table_prediction_final_mv_regxlevels) != 0)
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_treat)) > 1)
    nomogram_obj <- OUTPUT_DATA$table_prediction_log_nomogram
  } else if (input$nomogram_type == "EP_option") {
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_8_rec)) > 1)
    req(length(OUTPUT_DATA$table_prediction_final_mv_regxlevels_rec) != 0)
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_rec)) > 1)
    nomogram_obj <- OUTPUT_DATA$table_prediction_log_nomogram_rec
  } else if (input$nomogram_type == "_GMT") {
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_8_GMT)) > 1)
    req(length(OUTPUT_DATA$table_prediction_final_mv_regxlevels_GMT) != 0)
    req(length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_GMT)) > 1)
    nomogram_obj <- OUTPUT_DATA$table_prediction_log_nomogram_GMT
  } else if (input$nomogram_type == "post") {
    req(OUTPUT_DATA$table_prediction_tmp_post_log_nomogram)
    nomogram_obj <- OUTPUT_DATA$table_prediction_tmp_post_log_nomogram
  } else {
    # どの条件にも当てはまらない場合は空の文字列を返す
    return("")
  }

  # ノモグラムオブジェクトの Histology 部分が存在するか確認
  shiny::validate(
    shiny::need(!is.null(nomogram_obj$Histology$points), "No Histology point data")
  )

  # 共通のテキスト生成ロジック
  paste0("Histology and points\n",
         paste(paste(nomogram_obj$Histology$Histology,
                     format_p(nomogram_obj$Histology$points, digits = 0),
                     sep = "\t"),
               collapse = "\n"))
})

output$figure_nomogram_treat_pre_CGP <- renderPlot({
  req(input$nomogram_type)
  if (input$nomogram_type == "post") {
    handle_post_analysis()
  } else if (input$nomogram_type %in% c("EP_treat", "EP_option", "_GMT")) {
    config <- get_nomogram_config(input$nomogram_type)

    if (is_data_valid(config$data, config$data_treat, config$xlevels)) {
      par(mfcol = c(2, 1))
      plot_calibration(config)
      plot(config$nomogram)
    } else {
      gg_empty()
    }
  } else {
    gg_empty()
  }
})

# ファイルの読み込みまたは計算のヘルパー関数
load_or_compute <- function(file_path, compute_func) {
  if (input$new_analysis == "No, use the previous dataset" &&
      input$intermediate_file != "No" &&
      file.exists(file_path)) {
    return(QS_READ(nthreads = PARALLEL, file = file_path))
  } else {
    result <- try(compute_func(), silent = FALSE)
    if (class(result) == "try-error") {
      result <- try(compute_func(), silent = FALSE)
      if (class(result) == "try-error") {
        return(NULL)
      }
    }
    if (input$intermediate_file != "No") {
      QS_SAVE(nthreads = PARALLEL, result, file = file_path)
    }
    return(result)
  }
}

# データが有効かチェックする関数
is_data_valid <- function(data, data_treat, xlevels) {
  return(!is.null(data) &&
           !is.null(data_treat) &&
           !is.null(xlevels) &&
           length(colnames(data)) > 1 &&
           length(xlevels) != 0 &&
           length(colnames(data_treat)) > 1)
}

# post解析用の特別な処理
handle_post_analysis <- function() {
  if (is.null(OUTPUT_DATA$table_prediction_Data_forest_tmp_6) ||
      length(colnames(OUTPUT_DATA$table_prediction_Data_forest_tmp_6)) <= 1) {
    return(gg_empty())
  }

  # データ前処理
  Data_forest_tmp_6 <- data.frame(lapply(OUTPUT_DATA$table_prediction_Data_forest_tmp_6, as.factor))
  if ("Histology" %in% colnames(Data_forest_tmp_6)) {
    Data_forest_tmp_6$Histology <- relevel(
      Data_forest_tmp_6$Histology,
      ref = names(sort(table(Data_forest_tmp_6$Histology), decreasing = TRUE))[[1]]
    )
  }
  Data_forest_tmp_6 <- Data_forest_tmp_6[, colnames(Data_forest_tmp_6) != "Lines"]

  # モデル構築
  final_mv_regxlevels <- (prepare_data(
    as.numeric(as.character(Data_forest_tmp_6$EP_treat)),
    as.matrix(data.frame(lapply(Data_forest_tmp_6[, setdiff(names(Data_forest_tmp_6), "EP_treat")],
                                function(x) as.numeric(as.factor(x))))),
    type = "logistic"
  ) %>% stepwise(aic))$model

  if (length(final_mv_regxlevels) == 0) {
    return(gg_empty())
  }

  # モデル学習と検証
  Data_forest_tmp_treat <- Data_forest_tmp_6 %>% dplyr::select(c("EP_treat", final_mv_regxlevels))
  dd <- datadist(Data_forest_tmp_treat)
  options(datadist = dd)

  # モデルの読み込みまたは作成
  model_file <- file.path(tempdir(), "m2_post.qs")
  if (input$new_analysis == "No, use the previous dataset" &&
      input$intermediate_file != "No" &&
      file.exists(model_file)) {
    m2 <- QS_READ(nthreads = PARALLEL, file = model_file)
  } else {
    m2 <- lrm(as.formula(paste0("EP_treat ~ ", paste(colnames(Data_forest_tmp_treat)[colnames(Data_forest_tmp_treat) != "EP_treat"], collapse = "+"))),
              data = Data_forest_tmp_treat, x = TRUE, y = TRUE, penalty = Penal)
    if (input$intermediate_file != "No") {
      QS_SAVE(nthreads = PARALLEL, m2, file = model_file)
    }
  }

  BootNo = BootNoSet(Data_forest_tmp_treat)

  # バリデーションとキャリブレーション
  val <- load_or_compute(file.path(tempdir(), "val_post.qs"), function() rms::validate(m2, B = BootNo))
  cal <- load_or_compute(file.path(tempdir(), "cal_post.qs"), function() rms::calibrate(m2, B = BootNo))

  log_nomogram <- nomogram(m2, fun = plogis, lp = FALSE, funlabel = "Matched therapy")
  options(datadist = NULL)

  # プロット
  par(mfcol = c(2, 1))
  plot(cal, lwd = 2, lty = 1,
       cex.lab = 1.2, cex.axis = 1, cex.main = 1.2, cex.sub = 0.6,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nomogram-Predicted Probability of recommended treatment receive",
       ylab = "Actual recommended treatment receive (proportion)",
       legend = FALSE)

  stats_text <- paste0(
    "Matched therapy receiving ratio, only patients with CTx recommendation, post CGP information", "\n",
    "Model likelihood ratio test P-value: ", format_p(m2$stats["P"][[1]], digits = 3), '\n',
    BootNo, "-time-bootstrapped-Brier score: ", format_p(val[9, 5], digits = 3), '\n',
    BootNo, "-time-bootstrapped-Nagelkerke R2 = ", format_p(val[2, 5], digits = 3), '\n',
    BootNo, "-time-bootstrapped-concordance index = ", format_p(((val[1, 5] + 1) / 2), digits = 3), "\n",
    "Intercept = ", format_p(val[3, 5], digits = 3), "\n",
    "Slope = ", format_p(val[4, 5], digits = 3)
  )

  text(x = 0.4, y = 0.0, adj = c(0, 0), stats_text)
  abline(0, 1, lty = 3, lwd = 2, col = "#224444")
  legend(x = 0.8, y = 0.2, legend = c("Apparent", "Bias-corrected", "Ideal"),
         lty = c(3, 1, 2), bty = "n")
  plot(log_nomogram)

  # 結果をグローバル変数に保存
  OUTPUT_DATA$table_prediction_tmp_post_log_nomogram = log_nomogram
}

# 標準的なキャリブレーションプロットを描画する関数
plot_calibration <- function(config) {
  plot(config$cal, lwd = 2, lty = 1,
       cex.lab = 1.2, cex.axis = 1, cex.main = 1.2, cex.sub = 0.6,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = config$xlabel,
       ylab = config$ylabel,
       legend = FALSE)

  abline(0, 1, lty = 3, lwd = 2, col = "#224444")

  # 統計情報テキスト
  stats_text <- paste0(
    config$title, "\n",
    "Model likelihood ratio test P-value: ", format_p(config$pvalue, digits = 3), '\n',
    OUTPUT_DATA$table_prediction_BootNo, "-time-bootstrapped-Brier score: ", format_p(config$val[[1]], digits = 3), '\n',
    OUTPUT_DATA$table_prediction_BootNo, "-time-bootstrapped-Nagelkerke R2 = ", format_p(config$val[[2]], digits = 3), '\n',
    OUTPUT_DATA$table_prediction_BootNo, "-time-bootstrapped-concordance index = ", format_p(((config$val[[3]] + 1) / 2), digits = 3), "\n",
    "Intercept = ", format_p(config$val[[4]], digits = 3), "\n",
    "Slope = ", format_p(config$val[[5]], digits = 3)
  )

  text(x = 0.4, y = 0.0, adj = c(0, 0), stats_text)
  legend(x = 0.8, y = 0.2, legend = c("Apparent", "Bias-corrected", "Ideal"),
         lty = c(3, 1, 2), bty = "n")
}

# 共通パラメータを定義する関数
get_nomogram_config <- function(type) {
  configs <- list(
    EP_treat = list(
      data = OUTPUT_DATA$table_prediction_Data_forest_tmp_8,
      data_treat = OUTPUT_DATA$table_prediction_Data_forest_tmp_treat,
      xlevels = OUTPUT_DATA$table_prediction_final_mv_regxlevels,
      cal = OUTPUT_DATA$table_prediction_cal,
      nomogram = OUTPUT_DATA$table_prediction_log_nomogram,
      pvalue = OUTPUT_DATA$table_prediction_m2_pvalue,
      val = list(OUTPUT_DATA$table_prediction_m2_val_1,
                 OUTPUT_DATA$table_prediction_m2_val_2,
                 OUTPUT_DATA$table_prediction_m2_val_3,
                 OUTPUT_DATA$table_prediction_m2_val_4,
                 OUTPUT_DATA$table_prediction_m2_val_5),
      xlabel = "Nomogram-Predicted Probability of recommended treatment receive",
      ylabel = "Actual recommended treatment receive (proportion)",
      title = "Matched therapy receiving ratio, pre CGP information (other than mutation)"
    ),
    EP_option = list(
      data = OUTPUT_DATA$table_prediction_Data_forest_tmp_8_rec,
      data_treat = OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_rec,
      xlevels = OUTPUT_DATA$table_prediction_final_mv_regxlevels_rec,
      cal = OUTPUT_DATA$table_prediction_cal_rec,
      nomogram = OUTPUT_DATA$table_prediction_log_nomogram_rec,
      pvalue = OUTPUT_DATA$table_prediction_m2_rec_pvalue,
      val = list(OUTPUT_DATA$table_prediction_m2_rec_val_1,
                 OUTPUT_DATA$table_prediction_m2_rec_val_2,
                 OUTPUT_DATA$table_prediction_m2_rec_val_3,
                 OUTPUT_DATA$table_prediction_m2_rec_val_4,
                 OUTPUT_DATA$table_prediction_m2_rec_val_5),
      xlabel = "Nomogram-Predicted Probability of treatment recommendation",
      ylabel = "Actual treatment recommendation (proportion)",
      title = "Matched therapy recommendation ratio, pre CGP information (other than mutation)"
    ),
    `_GMT` = list(
      data = OUTPUT_DATA$table_prediction_Data_forest_tmp_8_GMT,
      data_treat = OUTPUT_DATA$table_prediction_Data_forest_tmp_treat_GMT,
      xlevels = OUTPUT_DATA$table_prediction_final_mv_regxlevels_GMT,
      cal = OUTPUT_DATA$table_prediction_cal_GMT,
      nomogram = OUTPUT_DATA$table_prediction_log_nomogram_GMT,
      pvalue = OUTPUT_DATA$table_prediction_m2_GMT_pvalue,
      val = list(OUTPUT_DATA$table_prediction_m2_GMT_val_1,
                 OUTPUT_DATA$table_prediction_m2_GMT_val_2,
                 OUTPUT_DATA$table_prediction_m2_GMT_val_3,
                 OUTPUT_DATA$table_prediction_m2_GMT_val_4,
                 OUTPUT_DATA$table_prediction_m2_GMT_val_5),
      xlabel = paste("Nomogram-Predicted Probability of Evidence level", paste(input$evidence_level, collapse = "/")),
      ylabel = paste("Actual evidence level", paste(input$evidence_level, collapse = "/"), "(proportion)"),
      title = "Matched therapy recommendation ratio, pre CGP information (other than mutation)"
    )
  )
  return(configs[[type]])
}

output$figure_survival_CGP_7 <- render_gt({
  # nomogram_type_2 に応じて、必要なデータと変数を動的に設定
  req(input$nomogram_type_2)

  params <- switch(input$nomogram_type_2,
                   "EP_treat" = list(
                     data_uni = OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant,
                     data_mv = OUTPUT_DATA$table_prediction_Data_forest_tmp_7,
                     outcome_var = "EP_treat",
                     file_uni = file.path(tempdir(), "m2_uni.qs"),
                     file_mv = file.path(tempdir(), "m2.qs"),
                     caption_base = "Factors that led patients to receive the recommended treatment"
                   ),
                   "EP_option" = list(
                     data_uni = OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant_rec,
                     data_mv = OUTPUT_DATA$table_prediction_Data_forest_tmp_7_rec,
                     outcome_var = "EP_option",
                     file_uni = file.path(tempdir(), "m2_rec_uni.qs"),
                     file_mv = file.path(tempdir(), "m2_rec.qs"),
                     caption_base = "Factors that recommended treatment"
                   ),
                   "_GMT" = list(
                     data_uni = OUTPUT_DATA$table_prediction_Data_forest_tmp_7__univariant_GMT,
                     data_mv = OUTPUT_DATA$table_prediction_Data_forest_tmp_7_GMT,
                     outcome_var = "EP_option",
                     file_uni = file.path(tempdir(), "m2_GMT_uni.qs"),
                     file_mv = file.path(tempdir(), "m2_GMT.qs"),
                     caption_base = paste0("Factors for genome-matched therapy (Evidence level ",
                                           paste(input$evidence_level, collapse = "/"),
                                           ")")
                   ),
                   "post" = list(
                     data_uni = OUTPUT_DATA$table_prediction_Data_forest_tmp_6_univariant,
                     data_mv = OUTPUT_DATA$table_prediction_Data_forest_tmp_6,
                     outcome_var = "EP_treat",
                     file_uni = file.path(tempdir(), "m2_post_uni.qs"),
                     file_mv = file.path(tempdir(), "m2_post.qs"),
                     caption_base = "Factors that led patients with recommended treatment to receive the recommended treatment"
                   )
  )

  # データフレームの列数が1より大きいことを確認
  req(length(colnames(params$data_uni)) > 1)

  withProgress(message = sample(nietzsche)[1], {

    # 共通のデータ前処理
    # 全ての列を因子型に変換
    params$data_uni <- data.frame(lapply(params$data_uni, as.factor))
    params$data_mv <- data.frame(lapply(params$data_mv, as.factor))

    # Histology列が存在する場合、リファレンスレベルを設定
    if ("Histology" %in% colnames(params$data_mv)) {
      ref_level <- names(sort(table(params$data_mv$Histology), decreasing = TRUE))[[1]]
      params$data_mv$Histology <- relevel(params$data_mv$Histology, ref = ref_level)
    }
    if ("Histology" %in% colnames(params$data_uni)) {
      ref_level_uni <- names(sort(table(params$data_uni$Histology), decreasing = TRUE))[[1]]
      params$data_uni$Histology <- relevel(params$data_uni$Histology, ref = ref_level_uni)
    }

    # 列名のエイリアスを定義
    rename_cols <- function(df, input) {
      colnames_tmp <- colnames(df)
      colnames_tmp[colnames_tmp == "Age"] <- paste0("Age ~", input$mid_age, "/", input$mid_age + 1, "~")
      colnames_tmp[colnames_tmp == "Lymph_met"] <- "Lymphatic metastasis"
      colnames_tmp[colnames_tmp == "Liver_met"] <- "Liver metastasis"
      colnames_tmp[colnames_tmp == "Lung_met"] <- "Lung metastasis"
      colnames_tmp[colnames_tmp == "Bone_met"] <- "Bone metastasis"
      colnames_tmp[colnames_tmp == "Brain_met"] <- "Brain metastasis"
      colnames_tmp[colnames_tmp == "PS"] <- "Performance status"
      colnames_tmp[colnames_tmp == "Lines"] <- "CTx lines before CGP"
      colnames_tmp[colnames_tmp == "Best_effect"] <- "Best CTx effect before CGP"
      colnames(df) <- colnames_tmp
      return(df)
    }

    # 列名を変更
    params$data_uni <- rename_cols(params$data_uni, input)

    incProgress(1 / 3)

    # 単変量解析の実行または読み込み
    if (input$new_analysis == "No, use the previous dataset" && input$intermediate_file != "No" && file.exists(params$file_uni)) {
      univ_tab <- QS_READ(nthreads = PARALLEL, file = params$file_uni)
    } else {
      run_univariate_model <- function(var, data, outcome_var) {
        tryCatch({
          fml <- as.formula(paste(outcome_var, "~", paste("`", var, "`", sep = "")))
          model <- glm(fml, data = data, family = binomial)
          tbl_regression(model, exponentiate = TRUE, conf.int = FALSE) %>%
            add_global_p() %>%
            add_n(location = "level") %>%
            add_nevent(location = "level") %>%
            add_q() %>%
            bold_p() %>%
            apply_custom_format(columns = c("estimate", "p.value", "q.value")) %>%
            bold_labels() %>%
            tbl_butcher()
        }, error = function(e) {
          return(NULL) # エラー時はNULLを返す
        })
      }

      univ_vars <- setdiff(names(params$data_uni), params$outcome_var)
      num_cores <- min(1 + 3*DOCKER, PARALLEL)
      # 結果格納用リスト
      univ_models_list <- list()

      # ---------------------------------------------------------
      # 並列処理の実行トライ
      # ---------------------------------------------------------
      run_serial <- TRUE # デフォルトで直列実行フラグを立てておく

      if(DOCKER && num_cores > 1) {
        # 並列実行
        univ_models_parallel <- tryCatch({
          parallel::mclapply(univ_vars, run_univariate_model,
                             data = params$data_uni,
                             outcome_var = params$outcome_var,
                             mc.cores = num_cores)
        }, warning = function(w) {
          # "core did not deliver a result" 警告が出たらここに来る
          return(NULL)
        }, error = function(e) {
          return(NULL)
        })

        # 結果の検証: NULLが含まれている、またはリストの長さが足りない場合は失敗とみなす
        # parallel::mclapplyは失敗した要素に "try-error" クラスを入れることがある
        is_valid <- !is.null(univ_models_parallel) &&
          length(univ_models_parallel) == length(univ_vars) &&
          !any(sapply(univ_models_parallel, function(x) inherits(x, "try-error")))

        if (is_valid) {
          univ_models_list <- univ_models_parallel
          run_serial <- FALSE # 並列成功したので直列はやらない
        } else {
          # 失敗時のログ（コンソールに出力）
          message("Parallel processing failed (core crash or memory issue). Switching to serial processing.")
        }
      }

      # ---------------------------------------------------------
      # 直列処理 (Parallelがオフ または Parallelが失敗した場合に実行)
      # ---------------------------------------------------------
      if (run_serial) {
        univ_models_list <- lapply(univ_vars, run_univariate_model,
                                   data = params$data_uni,
                                   outcome_var = params$outcome_var)
      }

      names(univ_models_list) <- univ_vars
      univ_models_clean <- Filter(function(x) !is.null(x) && inherits(x, "gtsummary"), univ_models_list)

      if (length(univ_models_clean) == 0) {
        validate(need(FALSE, "No univariate models converged."))
      }

      univ_tab <- tbl_stack(univ_models_clean)

      if (input$intermediate_file != "No") {
        QS_SAVE(nthreads = PARALLEL, univ_tab, file = params$file_uni)
      }
    }

    incProgress(1 / 3)

    # 多変量解析のためのステップワイズ選択とデータ整形
    final_mv_regxlevels <- (prepare_data(as.numeric(as.character(params$data_mv[[params$outcome_var]])),
                                         as.matrix(data.frame(lapply(params$data_mv[, setdiff(names(params$data_mv), params$outcome_var)],
                                                                     function(x) as.numeric(as.factor(x))))), type = "logistic") %>% stepwise(aic))$model

    # 多変量解析に選択された変数でデータフレームをフィルタリング
    params$data_mv <- params$data_mv %>% select(c(final_mv_regxlevels, params$outcome_var))

    # 列名を変更
    params$data_mv <- rename_cols(params$data_mv, input)
    params$data_mv <- params$data_mv[, colnames(params$data_mv) != "CTx lines before CGP"]

    incProgress(1 / 3)

    # 多変量解析の実行または読み込み
    if (length(final_mv_regxlevels) != 0) {
      if (input$new_analysis == "No, use the previous dataset" && input$intermediate_file != "No" && file.exists(params$file_mv)) {
        mv_tab <- QS_READ(nthreads = PARALLEL, file = params$file_mv)
      } else {
        mv_tab <- glm(as.formula(paste(params$outcome_var, "~ .")),
                      data = params$data_mv,
                      family = binomial(link = "logit")) %>%
          tbl_regression(exponentiate = TRUE, conf.int = FALSE) %>%
          add_global_p() %>%
          add_q() %>%
          bold_p() %>%
          apply_custom_format(columns = c("estimate", "p.value", "q.value")) %>%
          bold_labels() %>%
          tbl_butcher()

        if (input$intermediate_file != "No") {
          QS_SAVE(nthreads = PARALLEL, mv_tab, file = params$file_mv)
        }
      }

      # 単変量と多変量のテーブルを結合してキャプションを設定
      tbl_merge(
        tbls = list(univ_tab, mv_tab),
        tab_spanner = c("**Univariate**", "**Multivariable**")
      ) %>%
        modify_caption(paste0(params$caption_base, " (Only factors with >2 events observed)")) %>%
        as_gt()
    } else {
      # 多変量解析で有意な因子がない場合
      univ_tab %>%
        modify_caption(paste0(params$caption_base, ", no significant factor in multivariable analysis (Only factors with >2 events observed)")) %>%
        as_gt()
    }
  })
})
