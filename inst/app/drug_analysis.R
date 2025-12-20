drug_analysis_logic <- function() {
  analysis_env <- new.env()
  with(analysis_env, {
    if(!is.null(input$drug)){
      withProgress(message = sample(nietzsche)[1], {
        Data_case_target = Data_case()
        if(input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed",
                                            "Only cases without mutations in the gene set are analyzed")){
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
                        Date_start_1L,
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
        if(nrow(Data_case_target)>0){
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
          if(input$drug_volcano != "Yes"){
            Data_cluster_ID_list = Data_cluster_ID() %>%
              dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster, driver_mutations)
            Data_case_target = left_join(Data_case_target,
                                         Data_cluster_ID_list,
                                         by = "C.CAT調査結果.基本項目.ハッシュID")
            Data_case_target$cluster[is.na(Data_case_target$cluster)] = max(Data_case_target$cluster, na.rm = T) + 1
            Data_cluster_ID_list$cluster = as.factor(Data_cluster_ID_list$cluster)
            Data_case_target$cluster = as.factor(Data_case_target$cluster)
          } else {
            Data_case_target$cluster = "1"
          }

          incProgress(1 / 13)

          Data_survival = Data_case_target
          incProgress(1 / 13)

          Data_drug = Data_drug_raw_rename()
          Data_drug = Data_drug %>%
            dplyr::filter(ID %in%
                            Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
          Data_drug$AE = ifelse(Data_drug$Adverse_effect == "AE (-)", 0, 1)
          significant_genes = NULL
          Data_drug <- Data_drug %>%
            add_count(Drug, name = "drug_count") %>%
            mutate(Drug = ifelse(drug_count < ifelse(!is.null(input$minimum_courses), input$minimum_courses, 10), "Others", Drug)) %>%
            select(-drug_count)
          Data_drug$治療ライン = as.character(Data_drug$CTx_line)
          if(input$ToT_TTF == "ToT"){
            Data_drug$Drug_length = Data_drug$ToT
            Data_drug$Drug_length_censor = Data_drug$ToT_censor
          } else if(input$ToT_TTF == "TTF"){
            Data_drug$Drug_length = Data_drug$TTF
            Data_drug$Drug_length_censor = Data_drug$TTF_censor
          }  else if(input$ToT_TTF == "TTAE"){
            Data_drug$Drug_length = Data_drug$TTAE
            Data_drug$Drug_length_censor = Data_drug$TTAE_censor
          } else {
            Data_drug$Drug_length = Data_drug$TTD
            Data_drug$Drug_length_censor = Data_drug$censor
          }
          Data_drug$ToT[Data_drug$ToT >= 10000] = 10000
          Data_drug$TTF[Data_drug$TTF >= 10000] = 10000
          Data_drug$TTAE[Data_drug$TTAE >= 10000] = 10000
          Data_drug$TTD[Data_drug$TTD >= 10000] = 10000

          Data_drug = Data_drug %>%
            dplyr::mutate(Early_AE = case_when(
              TTAE > 0 & TTAE <= 30 ~ "Early AE (+)",
              TTAE > 30 & TTAE <= 90 ~ "Middle AE (+)",
              TTAE > 90 ~ "Late AE (+)",
              AE == 1 ~ "Unknown AE (+)",
              TRUE ~ "AE (-)"
            )
            )
          Data_drug_RECIST = Data_drug
          Data_drug_original = Data_drug
          if(input$multiple_used_drug == "Only the first administration of the treatment was included in the analysis"){
            Data_drug_RECIST = Data_drug_RECIST %>% dplyr::distinct(ID, Drug, .keep_all = T)
            Data_drug_original = Data_drug_original %>% dplyr::distinct(ID, Drug, .keep_all = T)
          }
          OUTPUT_DATA$drug_analysis_Data_drug_RECIST = Data_drug_RECIST

          Data_drug = Data_drug %>% dplyr::filter(Drug_length>0 & !is.na(Drug_length) & is.finite(Drug_length))
          Data_drug = Data_drug %>%
            dplyr::arrange(ID, 治療ライン)

          if(input$drug_volcano != "Yes"){
            Data_drug$Drug_length_pre_censor = -1
            Data_drug$Drug_length_pre = -1
            optimized_datatable <- function(Data_drug) {
              dt <- as.data.table(Data_drug)
              setorder(dt, ID)  # IDでソート
              dt[, Drug_length_pre := shift(Drug_length, type = "lag"), by = ID]
              dt[, Drug_length_pre_censor := shift(Drug_length_censor, type = "lag"), by = ID]
              return(dt)
            }
            Data_drug = optimized_datatable(Data_drug)
            Data_drug$Drug_length_pre[is.na(Data_drug$Drug_length_pre)] = -1
            Data_drug$Drug_length_pre_censor[is.na(Data_drug$Drug_length_pre_censor)] = -1
            Data_drug = Data_drug %>%
              dplyr::filter(治療ライン %in% input$target_line)
            if(input$multiple_used_drug == "Only the first administration of the treatment was included in the analysis"){
              Data_drug = Data_drug %>% dplyr::distinct(ID, Drug, .keep_all = T)
            }
            Drug_summary =  Data_drug %>%
              dplyr::filter(Drug %in% input$drug) %>%
              dplyr::select(
                ID,
                Drug,
                治療ライン,
                最良総合効果,
                終了理由,
                Adverse_effect,
                Overall_AE_grage,
                all_of(names(classification_rules)),
                Drug_length,
                Drug_length_censor,
                Early_AE
              ) %>%
              dplyr::distinct()

            colnames(Drug_summary) =
              c("ID", "Drugs", "CTx line", "RECIST", "Reason to finish treatment",
                "Any adverse effect (G3-G5)",
                "Max AE grade",
                names(classification_rules),
                "Time on treatment", "Treatment finished or censored", "AE within 1-month, 3-month, or later")
            Drug_summary$`Reason to finish treatment`[is.na(Drug_summary$`Reason to finish treatment`)] = "Unknown"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "その他理由で中止"] = "Other reason"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "不明"] = "Unknown"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "副作用等で中止"] = "Side effect"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "本人希望により中止"] = "Patient will"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "無効中止"] = "Not effective"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "死亡中止"] = "Death"
            Drug_summary$`Reason to finish treatment`[Drug_summary$`Reason to finish treatment` == "計画通り終了"] = "As planned"
            Drug_summary$`Treatment finished or censored` = as.character(Drug_summary$`Treatment finished or censored`)
            Drug_summary$`Treatment finished or censored`[is.na(Drug_summary$`Treatment finished or censored`)] = "Unknown"
            Drug_summary$`Treatment finished or censored`[Drug_summary$`Treatment finished or censored` == "1"] = "Finished"
            Drug_summary$`Treatment finished or censored`[Drug_summary$`Treatment finished or censored` == "0"] = "Censored"

            if(!"LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
              Data_case_age_sex = Data_case_target %>%
                dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                              症例.基本情報.年齢,
                              YoungOld,
                              症例.基本情報.性別.名称.,
                              症例.背景情報.喫煙歴有無.名称.,
                              症例.背景情報.アルコール多飲有無.名称.
                ) %>% distinct()
              colnames(Data_case_age_sex) = c(
                "ID",
                "Age",
                paste0("Age <=", input$mid_age, " or older"),
                "Sex",
                "Smoking history",
                "Heavy alcohol"
              )
            } else{
              Data_case_age_sex = Data_case_target %>%
                dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                              症例.基本情報.年齢,
                              YoungOld,
                              症例.基本情報.性別.名称.,
                              症例.背景情報.喫煙歴有無.名称.,
                              症例.背景情報.アルコール多飲有無.名称.,
                              PD_L1
                ) %>% distinct()
              colnames(Data_case_age_sex) = c(
                "ID",
                "Age",
                paste0("Age <=", input$mid_age, " or older"),
                "Sex",
                "Smoking history",
                "Heavy alcohol",
                "PD-L1 (lung only)"
              )
            }
            Drug_summary = Drug_summary %>%
              dplyr::left_join(Data_case_age_sex, by = "ID")

            mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
            if(length(mut_gene_)==0){
              mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
              ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
              Drug_summary_line = Drug_summary %>% dplyr::mutate(
                separation_value = case_when(
                  ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
                  TRUE ~ paste0(mut_gene_, " mut(-)")
                )
              )
            } else{
              if(is.null(input$gene_group_2)){
                ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
                Drug_summary_line = Drug_summary %>% dplyr::mutate(
                  separation_value = case_when(
                    ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
                    TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
                  )
                )
              } else{
                ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
                ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
                Drug_summary_line = Drug_summary %>% dplyr::mutate(
                  separation_value = case_when(
                    ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                    !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                  )
                )
              }
            }
            ID_diagnosis_list = Data_case_target %>% dplyr::select(C.CAT調査結果.基本項目.ハッシュID,症例.基本情報.がん種.OncoTree.) %>% dplyr::distinct()
            Drug_summary_line$diagnosis = unlist(lapply(list(Drug_summary_line$ID), function(x) {
              as.vector(ID_diagnosis_list$症例.基本情報.がん種.OncoTree.[match(x, ID_diagnosis_list$C.CAT調査結果.基本項目.ハッシュID)])}))
            remove_cols <- names(Drug_summary_line) %in% names(classification_rules) &
              sapply(Drug_summary_line, is.logical)
            Drug_summary_line <- Drug_summary_line[, !remove_cols, with = FALSE]
            keep_cols <- names(Drug_summary_line)[names(Drug_summary_line) %in% names(classification_rules)]
            Drug_summary_line <- Drug_summary_line %>%
              mutate(across(all_of(keep_cols), as.factor))



            mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
            if(length(mut_gene_)==0){
              mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
              ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
              Drug_summary_line = Drug_summary_line %>% dplyr::mutate(
                mutation = case_when(
                  ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
                  TRUE ~ paste0(mut_gene_, " mut(-)")
                )
              )
            } else{
              if(is.null(input$gene_group_2)){
                ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
                Drug_summary_line = Drug_summary_line %>% dplyr::mutate(
                  mutation = case_when(
                    ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
                    TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
                  )
                )
              } else{
                ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
                ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
                Drug_summary_line = Drug_summary_line %>% dplyr::mutate(
                  mutation = case_when(
                    ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                    !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                  )
                )
              }
            }
            incProgress(1 / 13)


            Data_drug$serial = 1:length(Data_drug$ID)

            Data_drug_TTF = Data_drug
            gb = list()
            gc = list()
            kk = 1
            Data_Drug_length_compare = Data_drug %>%
              dplyr::filter(Drug_length_pre>0)
            if(nrow(Data_Drug_length_compare)>0){
              Data_Drug_length_compare$Drug_length = as.numeric(Data_Drug_length_compare$Drug_length)
              Data_Drug_length_compare$Drug_length_pre = as.numeric(Data_Drug_length_compare$Drug_length_pre)
              Data_Drug_length_compare_tmp = Data_Drug_length_compare %>% dplyr::filter(Drug %in% input$drug) %>%
                dplyr::arrange(Drug_length)
              if(nrow(Data_Drug_length_compare_tmp)>0){
                Data_Drug_length_compare_tmp = Data_Drug_length_compare_tmp %>%
                  dplyr::arrange(Drug_length_pre)
                Data_tmp_1 = Data_Drug_length_compare_tmp %>% dplyr::select(
                  Drug_length, Drug_length_censor)
                Data_tmp_2 = Data_Drug_length_compare_tmp %>% dplyr::select(
                  Drug_length_pre, Drug_length_pre_censor)
                colnames(Data_tmp_2) = colnames(Data_tmp_1)
                Data_tmp_1$line = "Treatment time"
                Data_tmp_2$line = "Treatment time, previous-line"
                Data_tmp = rbind(Data_tmp_1, Data_tmp_2)
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~line, data=Data_tmp,type = "kaplan-meier", conf.type = "log-log")
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~line,
                                  data=Data_tmp, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~line,
                                  data=Data_tmp, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~line,
                               data=Data_tmp)
                OUTPUT_DATA$drug_analysis_gb_previous = surv_curv_drug(survfit_t, Data_tmp, paste0("Treatment time, ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
              }
            } else {
              OUTPUT_DATA$drug_analysis_gb_previous = NULL
            }
            incProgress(1 / 13)

            if(nrow(Data_drug_TTF) > 0){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~1, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              OUTPUT_DATA$drug_analysis_gb_all_selected_line = surv_curv_drug(survfit_t, Data_drug_TTF, paste0("Treatment time, all drugs"), NULL, NULL, NULL)
            }
            mut_gene = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
            if(is.null(input$gene)){
              mut_gene = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
            }

            Data_drug_TTF = Data_drug_TTF %>% dplyr::mutate(mut = case_when(
              ID %in% unique((Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene))$Tumor_Sample_Barcode) ~ "+",
              TRUE ~ "-"))
            if(length(unique(Data_drug_TTF$mut)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~mut, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_TTF) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~mut,
                                  data=Data_drug_TTF, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~mut,
                                  data=Data_drug_TTF, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~mut,
                               data=Data_drug_TTF)
                OUTPUT_DATA$drug_analysis_gb_mutation = surv_curv_drug(survfit_t, Data_drug_TTF, paste("Treatment time, all drugs,", paste0(collapse = ";", mut_gene) , "mut"), diff_0, diff_1 ,diff_2)
              }
            }
            if(length(unique(Data_drug_TTF$AE)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~AE, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_TTF) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~AE,
                                  data=Data_drug_TTF, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~AE,
                                  data=Data_drug_TTF, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~AE,
                               data=Data_drug_TTF)
                OUTPUT_DATA$drug_analysis_gb_AE = surv_curv_drug(survfit_t, Data_drug_TTF, paste("Treatment time, all drugs, adverse effect (G3~G5)"), diff_0, diff_1 ,diff_2)
              }
            }
            if(length(unique(Data_drug_TTF$Early_AE)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~Early_AE, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_TTF) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~Early_AE,
                                  data=Data_drug_TTF, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~Early_AE,
                                  data=Data_drug_TTF, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~Early_AE,
                               data=Data_drug_TTF)
                OUTPUT_DATA$drug_analysis_gb_Early_AE = surv_curv_drug(survfit_t, Data_drug_TTF, paste("Treatment time, all drugs, adverse effect (G3~G5) within 1-month, 3-month, or later"), diff_0, diff_1 ,diff_2)
              }
            }
            Data_drug_TTF = Data_drug_TTF %>% dplyr::mutate(Regimen_type = case_when(
              Drug %in% input$drug ~ paste(input$drug, collapse = ";"),
              TRUE ~ "Other drugs"))
            if(length(unique(Data_drug_TTF$Regimen_type)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~Regimen_type, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_TTF) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~Regimen_type,
                                  data=Data_drug_TTF, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~Regimen_type,
                                  data=Data_drug_TTF, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~Regimen_type,
                               data=Data_drug_TTF)
                OUTPUT_DATA$drug_analysis_gb_regimen = surv_curv_drug(survfit_t, Data_drug_TTF, paste0("Treatment time, ", paste(input$drug, collapse = ";")), diff_0, diff_1 ,diff_2)
              }
            }

            `%||%` <- function(a, b) if (!is.null(a)) a else b

            Data_drug_TTF <- Data_drug_TTF %>%
              dplyr::mutate(Regimen = case_when(
                Drug %in% (input$drug_group_1 %||% character(0)) ~ input$drug_group_1_name %||% "Other drugs",
                Drug %in% (input$drug_group_2 %||% character(0)) ~ input$drug_group_2_name %||% "Other drugs",
                Drug %in% (input$drug_group_3 %||% character(0)) ~ input$drug_group_3_name %||% "Other drugs",
                Drug %in% (input$drug_group_4 %||% character(0)) ~ input$drug_group_4_name %||% "Other drugs",
                TRUE ~ "Other drugs"
              ))
            Data_drug_Drug_length_drug_select = Data_drug_TTF %>% dplyr::filter(Drug != "Other drugs")
            if(length(unique(Data_drug_Drug_length_drug_select$Regimen)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~Regimen, data=Data_drug_Drug_length_drug_select,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_Drug_length_drug_select) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~Regimen,
                                  data=Data_drug_Drug_length_drug_select, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~Regimen,
                                  data=Data_drug_Drug_length_drug_select, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~Regimen,
                               data=Data_drug_Drug_length_drug_select)
                OUTPUT_DATA$drug_analysis_gb_regimen_detail = surv_curv_drug(survfit_t, Data_drug_Drug_length_drug_select, "Treatment time, selected drugs", diff_0, diff_1, diff_2)
              }
            }

            if(nrow(Data_drug_TTF) > 0){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~mut+Regimen_type, data=Data_drug_TTF,type = "kaplan-meier", conf.type = "log-log")
              if(length(unique(Data_drug_TTF$mut)) > 1 & length(unique(Data_drug_TTF$Regimen_type)) > 1){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~mut+Regimen_type,
                                  data=Data_drug_TTF, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~mut+Regimen_type,
                                  data=Data_drug_TTF, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~mut+Regimen_type,
                               data=Data_drug_TTF)
                OUTPUT_DATA$drug_analysis_gb_mutation_regimen = surv_curv_drug(survfit_t, Data_drug_TTF, paste("Treatment time,", paste0(collapse = ";", mut_gene), "and", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
              }
            }

            # Survival analysis for mutation, all
            ID_mutation = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene))$Tumor_Sample_Barcode
            for(jj in 1:length(mut_gene)){
              ID_mutation = ID_mutation[ID_mutation %in% (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene[jj]))$Tumor_Sample_Barcode]
            }
            if(length(ID_mutation) > 0){
              Data_drug_TTF = Data_drug_TTF %>% dplyr::mutate(
                mut = case_when(
                  ID %in% ID_mutation ~ "All genes",
                  ID %in% (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene))$Tumor_Sample_Barcode ~ "Any genes",
                  TRUE ~ "No genes"
                )
              )
              Data_drug_Drug_length_mut = Data_drug_TTF %>% dplyr::filter(Regimen_type == paste(input$drug, collapse = ";"))
              if(length(unique(Data_drug_Drug_length_mut$mut)) > 1){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~mut, data=Data_drug_Drug_length_mut,type = "kaplan-meier", conf.type = "log-log")
                if(nrow(Data_drug_Drug_length_mut) != survfit_t[[1]][[1]]){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~mut,
                                    data=Data_drug_Drug_length_mut, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~mut,
                                    data=Data_drug_Drug_length_mut, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~mut,
                                 data=Data_drug_Drug_length_mut)
                  OUTPUT_DATA$drug_analysis_gb_mutation_detail = surv_curv_drug(survfit_t, Data_drug_Drug_length_mut, paste("Treatment time,", paste0(collapse = ";", mut_gene), "and", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
                }
              }
            }
            # Survival analysis for lymph node metastasis
            ID_mutation_lymph_yes = (Data_case_target %>% dplyr::filter(Lymph_met %in% c("Yes")))$C.CAT調査結果.基本項目.ハッシュID
            ID_mutation_lymph_no = (Data_case_target %>% dplyr::filter(Lymph_met %in% c("No")))$C.CAT調査結果.基本項目.ハッシュID
            Data_drug_TTF = Data_drug_TTF %>% dplyr::mutate(
              met = case_when(
                ID %in% ID_mutation_lymph_yes ~ "lymph met (+)",
                ID %in% ID_mutation_lymph_no ~ "lymph met (-)",
                TRUE ~ "Unknown"
              )
            )
            Data_drug_Drug_length_met = Data_drug_TTF %>%
              dplyr::filter(Regimen_type == paste(input$drug, collapse = ";")) %>%
              dplyr::filter(met != "Unknown")
            if(length(unique(Data_drug_Drug_length_met$met)) > 1){
              survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~met, data=Data_drug_Drug_length_met,type = "kaplan-meier", conf.type = "log-log")
              if(nrow(Data_drug_Drug_length_met) != survfit_t[[1]][[1]]){
                diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~met,
                                  data=Data_drug_Drug_length_met, rho=0)
                diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~met,
                                  data=Data_drug_Drug_length_met, rho=1)
                diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~met,
                               data=Data_drug_Drug_length_met)
                OUTPUT_DATA$drug_analysis_gb_lymph = surv_curv_drug(survfit_t, Data_drug_Drug_length_met, paste("Lymph node metastasis, treatment time with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
                kk = kk + 1
              }
            }
            Data_drug_TTF$diagnosis = ""
            for(target_diseases in unique(Data_case_target$症例.基本情報.がん種.OncoTree.)){
              Data_drug_TTF = Data_drug_TTF %>% dplyr::mutate(diagnosis = case_when(
                ID %in% unique((Data_case_target %>% dplyr::filter(症例.基本情報.がん種.OncoTree. == target_diseases))$C.CAT調査結果.基本項目.ハッシュID) ~ target_diseases,
                TRUE ~ diagnosis))
            }
            # # 副作用の列名
            # # 1. 想定される副作用列名（全体）
            # side_effect_cols_all <- names(classification_rules)
            # # 実在する列を選び、0/1変換
            # SE_data_all <- Drug_summary %>%
            #   dplyr::filter(`Any adverse effect (G3-G5)` != "AE (-)") %>%
            #   select(any_of(side_effect_cols_all)) %>%
            #   mutate(across(everything(), ~ ifelse(!is.na(.), 1, 0)))
            # # 各列の合計が1以上の列だけを選択
            # valid_cols <- names(SE_data_all)[colSums(SE_data_all, na.rm = TRUE) >= 1]
            # if(length(valid_cols)>1){
            #   # フィルタ済みの副作用データ
            #   SE_data <- SE_data_all %>% select(all_of(valid_cols))
            #   # 2. 相関係数と p 値の計算
            #   cor_result <- cor(SE_data, method = "kendall", use = "pairwise.complete.obs")
            #   p_mat <- matrix(NA, ncol = ncol(SE_data), nrow = ncol(SE_data))
            #   colnames(p_mat) <- rownames(p_mat) <- colnames(cor_result)
            #   # 3. p値の行列を作成
            #   for (i in 1:ncol(SE_data)) {
            #     for (j in 1:ncol(SE_data)) {
            #       test <- cor.test(SE_data[[i]], SE_data[[j]],method = "kendall")
            #       p_mat[i, j] <- test$p.value
            #     }
            #   }
            #   # 4. データ整形
            #   cor_melt <- melt(cor_result)
            #   p_melt <- melt(p_mat)
            #   heat_data <- cor_melt %>%
            #     rename(Cor = value) %>%
            #     mutate(p_value = p_melt$value) %>%
            #     mutate(label = sprintf("%.2f\n(p=%.3f)", Cor, p_value))
            #   heat_data$Var1 <- factor(heat_data$Var1, levels = unique(heat_data$Var1))
            #   heat_data$Var2 <- factor(heat_data$Var2, levels = (unique(heat_data$Var2)))
            #   heat_data <- heat_data[
            #     as.numeric(heat_data$Var2) >= as.numeric(heat_data$Var1),
            #   ]
            #   heat_data$Var2 <- factor(heat_data$Var2, levels = rev(unique(heat_data$Var2)))
            #   # 5. ヒートマッププロット
            #   OUTPUT_DATA$drug_analysis_gb_AE_heatmap = ggplot(heat_data, aes(x = Var1, y = Var2, fill = Cor)) +
            #     geom_tile(color = "white") +
            #     geom_text(aes(label = label), size = 3) +
            #     scale_fill_gradient2(
            #       low = "blue", high = "red", mid = "white", midpoint = 0,
            #       limit = c(-1, 1), name = "correlation"
            #     ) +
            #     theme_minimal()+theme(panel.grid=element_blank()) +
            #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            #     labs(
            #       x = "Adverse effect 1", y = "Adverse effect 2",
            #       title = "Concurrent side effects correlation matrix",
            #       subtitle = "Only patients with adverse effect (G3~5)"
            #     )
            # }



            incProgress(1 / 13)

            candidate_genes = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                                            paste0(input$special_gene, "_NOS"))))
            candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
            OUTPUT_DATA$drug_analysis_candidate_genes = candidate_genes[!is.na(candidate_genes)]
            Top_gene = unique(names(sort(table(Data_report_raw()$Hugo_Symbol), decreasing = T)))
            OUTPUT_DATA$drug_analysis_Top_gene = Top_gene[Top_gene %in% candidate_genes]
            candidate_drugs = sort(unique(c(Data_drug_TTF$Drug)))
            OUTPUT_DATA$drug_analysis_candidate_drugs = candidate_drugs[!is.na(candidate_drugs)]
            OUTPUT_DATA$drug_analysis_candidate_Overall_AE_grage = sort(unique(Data_drug_TTF$Overall_AE_grage))
            OUTPUT_DATA$drug_analysis_candidate_lines = sort(unique(Data_drug_TTF$CTx_line))
            OUTPUT_DATA$drug_analysis_candidate_Age = sort(unique(Data_case_target$YoungOld))
            OUTPUT_DATA$drug_analysis_candidate_RECIST = c("CR"=5,"PR"=4,"SD"=3,"PD"=2,"NE"=1)
            OUTPUT_DATA$drug_analysis_candidate_Sex = sort(unique(Data_case_target$症例.基本情報.性別.名称.))
            OUTPUT_DATA$drug_analysis_candidate_Histology = sort(unique(Data_case_target$症例.基本情報.がん種.OncoTree.))
            OUTPUT_DATA$drug_analysis_candidate_cluster = sort(unique(Data_case_target$cluster))
            OUTPUT_DATA$drug_analysis_candidate_meta = c('Lymph_met','Brain_met','Lung_met','Bone_met','Liver_met')
            OUTPUT_DATA$drug_analysis_candidate_AE = names(classification_rules)


            OUTPUT_DATA$drug_analysis_Data_drug_TTF_comparison = Data_drug_TTF
            incProgress(1 / 13)

            k_3 = 1
            if(length(unique(Data_drug_TTF$diagnosis)) >= 1){
              diagnosis_list = (Data_drug_TTF %>%
                                  dplyr::distinct(ID, .keep_all = T))$diagnosis
              diagnosis_list = names(sort(table(diagnosis_list),decreasing = T))[1:7]
              Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
                diagnosis = case_when(
                  diagnosis %in% diagnosis_list ~ diagnosis,
                  TRUE ~ "Other histology"
                )
              )
              Data_drug_Drug_length_tmp2 = Data_drug_Drug_length_tmp
              Data_drug_Drug_length_tmp2$diagnosis = " ALL"
              Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(
                diagnosis %in% diagnosis_list)
              Data_drug_Drug_length_tmp2 = rbind(Data_drug_Drug_length_tmp, Data_drug_Drug_length_tmp2)
              Data_drug_Drug_length_tmp$diagnosis = as.factor(Data_drug_Drug_length_tmp$diagnosis)
              Data_drug_Drug_length_tmp2$diagnosis = as.factor(Data_drug_Drug_length_tmp2$diagnosis)
              Data_drug_Drug_length_tmp$diagnosis = relevel(Data_drug_Drug_length_tmp$diagnosis, ref=names(sort(table(Data_drug_Drug_length_tmp$diagnosis),decreasing = T))[[1]])
              Data_drug_Drug_length_tmp2$diagnosis = relevel(Data_drug_Drug_length_tmp2$diagnosis, ref=names(sort(table(Data_drug_Drug_length_tmp2$diagnosis),decreasing = T))[[1]])
              if(nrow(Data_drug_Drug_length_tmp) > 0){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~diagnosis, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
                if(length(unique(Data_drug_Drug_length_tmp$diagnosis)) > 1){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                 data=Data_drug_Drug_length_tmp)
                  OUTPUT_DATA$drug_analysis_gc_diagnosis = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, diagnosis, all drugs"), diff_0, diff_1, diff_2)
                } else{
                  OUTPUT_DATA$drug_analysis_gc_diagnosis = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, diagnosis, all drugs"), NULL, NULL, NULL)
                }
              }
              if(nrow(Data_drug_Drug_length_tmp2) > 0){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~diagnosis, data=Data_drug_Drug_length_tmp2,type = "kaplan-meier", conf.type = "log-log")
                if(length(unique(Data_drug_Drug_length_tmp2$diagnosis)) > 1){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp2, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp2, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                 data=Data_drug_Drug_length_tmp2)
                  OUTPUT_DATA$drug_analysis_gc_diagnosis_all = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp2, paste0("Treatment time, diagnosis, all drugs"), diff_0, diff_1, diff_2)
                  k_3 = k_3 + 1
                } else{
                  OUTPUT_DATA$drug_analysis_gc_diagnosis_all = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp2, paste0("Treatment time, diagnosis, all drugs"), NULL, NULL,NULL)
                  k_3 = k_3 + 1
                }
              }

              Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(Drug %in% input$drug)
              Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::filter(
                diagnosis %in% diagnosis_list)
              if(length(unique(Data_drug_Drug_length_tmp$diagnosis)) >= 1){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~diagnosis, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
                if(length(unique(Data_drug_Drug_length_tmp$diagnosis)) > 1){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                    data=Data_drug_Drug_length_tmp, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~diagnosis,
                                 data=Data_drug_Drug_length_tmp)
                  OUTPUT_DATA$drug_analysis_gc_regimen = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, diagnosis, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
                  k_3 = k_3 + 1
                } else {
                  OUTPUT_DATA$drug_analysis_gc_regimen = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, diagnosis, treated with ", paste(input$drug, collapse = ";")), NULL, NULL)
                  k_3 = k_3 + 1
                }
              }
              Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(Drug %in% input$drug)
              if(length(unique(Data_drug_Drug_length_tmp$AE)) >= 1){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~AE, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
                if(length(unique(Data_drug_Drug_length_tmp$AE)) > 1){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~AE,
                                    data=Data_drug_Drug_length_tmp, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~AE,
                                    data=Data_drug_Drug_length_tmp, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~AE,
                                 data=Data_drug_Drug_length_tmp)
                  OUTPUT_DATA$drug_analysis_gc_AE = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, adverse effect (G3-G5), treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
                } else {
                  OUTPUT_DATA$drug_analysis_gc_AE = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, adverse effect (G3-G5), treated with ", paste(input$drug, collapse = ";")), NULL, NULL)
                }
              }
              if(length(unique(Data_drug_Drug_length_tmp$Early_AE)) >= 1){
                survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~Early_AE, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
                if(length(unique(Data_drug_Drug_length_tmp$Early_AE)) > 1){
                  diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~Early_AE,
                                    data=Data_drug_Drug_length_tmp, rho=0)
                  diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~Early_AE,
                                    data=Data_drug_Drug_length_tmp, rho=1)
                  diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~Early_AE,
                                 data=Data_drug_Drug_length_tmp)
                  OUTPUT_DATA$drug_analysis_gc_EAE = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, adverse effect within 1-month, 3-month, or later, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
                } else {
                  OUTPUT_DATA$drug_analysis_gc_EAE = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, adverse effect within 1-month, 3-month, or later, treated with ", paste(input$drug, collapse = ";")), NULL, NULL)
                }
              }
              if(!is.null(input$gene_group_1)){
                if(!is.null(input$gene_group_2)){
                  Data_MAF_target_tmp = Data_MAF_target %>%
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
                  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                    mutation_type = case_when(
                      ID %in% ID_special_gene_mutation_1 &
                        ID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut, ", paste0(paste(input$gene_group_2, collapse="/"), " mut")),
                      ID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut, ", paste0(paste(input$gene_group_2, collapse="/"), " WT")),
                      ID %in% ID_special_gene_mutation_2 ~ paste0(paste(input$gene_group_1, collapse="/"), " WT, ", paste0(paste(input$gene_group_2, collapse="/"), " mut")),
                      TRUE ~ paste0("No mut in ",  paste(input$gene_group_1, collapse="/"), "/", paste(input$gene_group_2, collapse="/"))
                    )
                  )
                } else {
                  Data_MAF_target_tmp = Data_MAF_target %>%
                    dplyr::filter(Tumor_Sample_Barcode %in%
                                    Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
                    dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)

                  ID_special_gene_mutation_1 = (Data_MAF_target_tmp %>%
                                                  dplyr::filter(Hugo_Symbol %in% input$gene_group_1 | str_detect(Hugo_Symbol, paste(paste0(input$gene_group_1, "_"),collapse ="|")) &
                                                                  amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
                  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                    mutation_type = case_when(
                      ID %in% ID_special_gene_mutation_1 ~ paste0(paste(input$gene_group_1, collapse="/"), " mut"),
                      TRUE ~ paste0("No mut in ", paste(input$gene_group_1, collapse="/"))
                    )
                  )
                }
              } else if(!input$special_gene == ""){
                Data_MAF_target_tmp = Data_MAF_target %>%
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
                Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                  mutation_type = case_when(
                    ID %in% ID_special_gene_mutation_1 ~ paste0(input$special_gene_mutation_1_name, " in ", input$special_gene),
                    ID %in% ID_special_gene_mutation_2 ~ paste0(input$special_gene_mutation_2_name, " in ", input$special_gene),
                    ID %in% ID_special_gene ~ paste0("Other mut in ", input$special_gene),
                    TRUE ~ paste0("No mut in ", input$special_gene)
                  )
                )
              } else {
                Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                  mutation_type = "No mutation pattern")
              }
            }

            Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(Drug %in% input$drug)
            oncogenic_genes = Data_drug_Drug_length_tmp %>%
              dplyr::select(ID, diagnosis) %>%
              dplyr::distinct() %>%
              dplyr::count(diagnosis) %>%
              dplyr::arrange(-n)
            colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
            if(length(unique(Data_drug_TTF$diagnosis)) >= 1 & nrow(Data_drug_Drug_length_tmp) > 0 & nrow(oncogenic_genes) > 0){
              oncogenic_genes = oncogenic_genes[1:min(nrow(oncogenic_genes), 10),]
              gene_table = data.frame(oncogenic_genes$gene_mutation)
              colnames(gene_table) = c("Gene")
              gene_table$positive_patients = oncogenic_genes$all_patients
              gene_table$positive_freq = format_p(100 * gene_table$positive_patients / length(unique(Data_drug_Drug_length_tmp$ID)), digits = 1)
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
              withProgress(message = "analysis for diaagnosis", {
                for(i in 1:nrow(oncogenic_genes)){
                  ID_mutation = (Data_drug_TTF %>% dplyr::filter(diagnosis %in% oncogenic_genes$gene_mutation[i]))$ID
                  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                    gene_mut = case_when(
                      ID %in% ID_mutation ~ 1,
                      TRUE ~ 0
                    )
                  )
                  traditional_fit = survfit(Surv(event = Drug_length_censor,
                                                 time = Drug_length) ~ gene_mut,
                                            data = Data_drug_Drug_length_tmp)
                  if(nrow(Data_drug_Drug_length_tmp) != traditional_fit[[1]][[1]]){
                    tau0 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 0))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau1 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 1))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau = floor(min(tau0, tau1, input$RMST_drug, na.rm = T) * 10) / 10
                    verify <- survRM2::rmst2(
                      time = Data_drug_Drug_length_tmp$Drug_length,
                      status = Data_drug_Drug_length_tmp$Drug_length_censor,
                      arm = Data_drug_Drug_length_tmp$gene_mut,
                      tau = tau * 365.25 / 12,
                      alpha = 0.05
                    )

                    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=0)
                    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=1)
                    mut_gene = oncogenic_genes$gene_mutation[i]
                    tmp = data.frame(summary(traditional_fit)$table)
                    gene_table$positive_median[i] = round(tmp$median[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_LL[i] = round(tmp$X0.95LCL[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_UL[i] = round(tmp$X0.95UCL[2] / 365.25 * 12, digits = 1)
                    gene_table$negative_median[i] = round(tmp$median[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp$X0.95LCL[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp$X0.95UCL[1] / 365.25 * 12, digits = 1)
                    gene_table$diff_median[i] = round(verify$unadjusted.result[1] / 365.25 * 12, digits = 2)
                    gene_table$diff_UL[i] = round(verify$unadjusted.result[7] / 365.25 * 12, digits = 2)
                    gene_table$diff_LL[i] = round(verify$unadjusted.result[4] / 365.25 * 12, digits = 2)
                  } else{
                    tmp = summary(traditional_fit)$table
                    legends = paste0(format_p(tmp[[7]] / 365.25 * 12, digits = 1), " (",
                                     format_p(tmp[[8]] / 365.25 * 12, digits = 1), "-",
                                     format_p(tmp[[9]] / 365.25 * 12, digits = 1),")")
                    gene_table$negative_median[i] = round(tmp[[7]] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp[[8]] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp[[9]] / 365.25 * 12, digits = 1)
                  }
                  incProgress(1/nrow(oncogenic_genes))
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
                labs(x=paste0("Median treatment time and difference in restricted mean survival time in ", input$RMST_drug, " months with ", paste(input$drug, collapse = ";"), " (months)"), y="Diagnosis") +
                coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                                xlim=c( - max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5,
                                        max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)) +
                annotate("text",
                         x = (-max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) - 0.5)/2,
                         y = length(Gene_arrange) + 1,
                         label = "diagnosis=short survival") +
                annotate("text",
                         x = (max(abs(gene_table$diff_LL), abs(gene_table$diff_UL)) + 0.5)/2,
                         y = length(Gene_arrange) + 1,
                         label = "diagnosis=long survival") +
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
                    Gene = "Dignosis",
                    estimate_lab_1 = "ToT, diagnosis",
                    estimate_lab_2 = "ToT, others",
                    estimate_lab_3 = "RMST difference",
                    patients = "Patients"
                  ), gene_table)
              Gene_arrange = gene_table$Gene
              gene_table = gene_table %>% dplyr::mutate(
                Gene = factor(gene_table$Gene, levels = Gene_arrange))

              p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
              p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                           fontface = ifelse(gene_table$Gene == "Dignosis", "bold", "bold.italic"))
              p_left <- p_left + geom_text(aes(x = 1.8, label = estimate_lab_1), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_1 == "ToT, diagnosis", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 3.2, label = estimate_lab_2), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_2 == "ToT, others", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 4.6, label = estimate_lab_3), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_3 == "RMST difference", "bold", "plain"))
              p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 7.2))

              # right side of plot - pvalues
              p_right <- gene_table  |> ggplot() +
                geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                          fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
                theme_void()

              # final plot arrangement
              layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 7.2),
                          patchwork::area(t = 1, l = 7.2, b = 30, r = 12.5),
                          patchwork::area(t = 0, l = 12.5, b = 30, r = 14))
              OUTPUT_DATA$drug_analysis_has = p_left + p_mid + p_right + plot_layout(design = layout)
              OUTPUT_DATA$drug_analysis_gene_table_3 = gene_table
            }

            incProgress(1 / 13)

            # analysis for common oncogenic mutations
            Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(Drug %in% input$drug)
            hs = list()
            hsc = list()

            oncogenic_genes = Data_MAF_target %>%
              dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
              dplyr::filter(Tumor_Sample_Barcode %in% Data_drug_Drug_length_tmp$ID) %>%
              dplyr::distinct() %>%
              dplyr::count(Hugo_Symbol) %>%
              dplyr::arrange(-n)
            colnames(oncogenic_genes) = c("gene_mutation", "all_patients")
            oncogenic_genes = rbind(oncogenic_genes %>% dplyr::filter(gene_mutation %in% input$gene),
                                    oncogenic_genes %>% dplyr::filter(!gene_mutation %in% input$gene))
            if(nrow(Data_drug_Drug_length_tmp) > 0 & nrow(oncogenic_genes)>0){
              oncogenic_genes = oncogenic_genes[1:min(nrow(oncogenic_genes), input$gene_no),]
              gene_table = data.frame(oncogenic_genes$gene_mutation)
              colnames(gene_table) = c("Gene")
              gene_table$positive_patients = oncogenic_genes$all_patients
              gene_table$positive_freq = format_p(100 * gene_table$positive_patients / length(unique(Data_drug_Drug_length_tmp$ID)), digits = 1)
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
              withProgress(message = "analysis for common oncogenic mutations", {
                for(i in 1:nrow(oncogenic_genes)){
                  ID_mutation = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% oncogenic_genes$gene_mutation[i]))$Tumor_Sample_Barcode
                  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                    gene_mut = case_when(
                      ID %in% ID_mutation ~ 1,
                      TRUE ~ 0
                    )
                  )
                  traditional_fit = survfit(Surv(event = Drug_length_censor,
                                                 time = Drug_length) ~ gene_mut,
                                            data = Data_drug_Drug_length_tmp)
                  if(nrow(Data_drug_Drug_length_tmp) != traditional_fit[[1]][[1]]){
                    tau0 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 0))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau1 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 1))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau = floor(min(tau0, tau1, input$RMST_drug, na.rm = T) * 10) / 10
                    verify <- survRM2::rmst2(
                      time = Data_drug_Drug_length_tmp$Drug_length,
                      status = Data_drug_Drug_length_tmp$Drug_length_censor,
                      arm = Data_drug_Drug_length_tmp$gene_mut,
                      tau = tau * 365.25 / 12,
                      alpha = 0.05
                    )

                    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=0)
                    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=1)
                    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                   data=Data_drug_Drug_length_tmp)
                    mut_gene = paste(oncogenic_genes$gene_mutation[i],
                                     ", top ", i, " gene", sep="")
                    tmp = data.frame(summary(traditional_fit)$table)
                    gene_table$positive_median[i] = round(tmp$median[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_LL[i] = round(tmp$X0.95LCL[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_UL[i] = round(tmp$X0.95UCL[2] / 365.25 * 12, digits = 1)
                    gene_table$negative_median[i] = round(tmp$median[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp$X0.95LCL[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp$X0.95UCL[1] / 365.25 * 12, digits = 1)
                    gene_table$diff_median[i] = round(verify$unadjusted.result[1] / 365.25 * 12, digits = 2)
                    gene_table$diff_UL[i] = round(verify$unadjusted.result[7] / 365.25 * 12, digits = 2)
                    gene_table$diff_LL[i] = round(verify$unadjusted.result[4] / 365.25 * 12, digits = 2)
                    hs[[k]] = surv_curv_drug(traditional_fit, Data_drug_Drug_length_tmp, paste(oncogenic_genes$gene_mutation[i], "mut, treated with", paste(input$drug, collapse = ";")),
                                             diff_0, diff_1, diff_2)
                    k = k + 1
                  } else{
                    tmp = summary(traditional_fit)$table
                    legends = paste0(format_p(tmp[[7]] / 365.25 * 12, digits = 1), " (",
                                     format_p(tmp[[8]] / 365.25 * 12, digits = 1), "-",
                                     format_p(tmp[[9]] / 365.25 * 12, digits = 1),")")
                    gene_table$negative_median[i] = round(tmp[[7]] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp[[8]] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp[[9]] / 365.25 * 12, digits = 1)
                  }
                  incProgress(1/nrow(oncogenic_genes))
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
                labs(x=paste0("Median treatment time and difference in restricted mean survival time in ", input$RMST_drug, " months with ", paste(input$drug, collapse = ";"), " (months)"), y="Genes with oncogenic alteration") +
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
                    estimate_lab_1 = "ToT, mut (+)",
                    estimate_lab_2 = "ToT, mut (-)",
                    estimate_lab_3 = "RMST difference",
                    patients = "Patients"
                  ), gene_table)
              Gene_arrange = gene_table$Gene
              gene_table = gene_table %>% dplyr::mutate(
                Gene = factor(gene_table$Gene, levels = Gene_arrange))

              p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
              p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                           fontface = ifelse(gene_table$Gene == "Gene", "bold", "bold.italic"))
              p_left <- p_left + geom_text(aes(x = 2.3, label = estimate_lab_1), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_1 == "ToT, mut (+)", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 4.0, label = estimate_lab_2), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_2 == "ToT, mut (-)", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 5.0, label = estimate_lab_3), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_3 == "RMST difference", "bold", "plain"))
              p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 7.5))

              # right side of plot - pvalues
              p_right <- gene_table  |> ggplot() +
                geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                          fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
                theme_void()
              # final plot arrangement
              layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 7.5),
                          patchwork::area(t = 1, l = 7.5, b = 30, r = 12.5),
                          patchwork::area(t = 0, l = 12.5, b = 30, r = 14))
              OUTPUT_DATA$drug_analysis_h = p_left + p_mid + p_right + plot_layout(design = layout)

              for(kk in k:32){
                hs[[kk]] = ggsurvplot_empty()
              }
              OUTPUT_DATA$drug_analysis_gene_table_gene = gene_table
            }
            ####################
            # analysis for clusters
            Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::filter(Drug %in% input$drug)
            oncogenic_genes = Data_case_target %>%
              dplyr::select(C.CAT調査結果.基本項目.ハッシュID, cluster) %>%
              dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% Data_drug_Drug_length_tmp$ID) %>%
              dplyr::distinct() %>%
              dplyr::count(cluster) %>%
              dplyr::arrange(cluster)
            colnames(oncogenic_genes) = c("cluster", "all_patients")

            if(nrow(Data_drug_Drug_length_tmp) > 0 & nrow(oncogenic_genes) > 0){
              gene_table = data.frame(oncogenic_genes$cluster)
              colnames(gene_table) = c("Gene")
              gene_table$positive_patients = oncogenic_genes$all_patients
              gene_table$positive_freq = format_p(100 * gene_table$positive_patients / length(unique(Data_drug_Drug_length_tmp$ID)), digits = 1)
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
              withProgress(message = "analysis for mutational clusters", {
                for(i in 1:nrow(oncogenic_genes)){
                  ID_mutation = (Data_case_target %>% dplyr::filter(cluster == oncogenic_genes$cluster[i]))$C.CAT調査結果.基本項目.ハッシュID
                  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
                    gene_mut = case_when(
                      ID %in% ID_mutation ~ 1,
                      TRUE ~ 0
                    )
                  )
                  traditional_fit = survfit(Surv(event = Drug_length_censor,
                                                 time = Drug_length) ~ gene_mut,
                                            data = Data_drug_Drug_length_tmp)
                  if(nrow(Data_drug_Drug_length_tmp) != traditional_fit[[1]][[1]]){
                    tau0 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 0))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau1 = ((max((Data_drug_Drug_length_tmp %>% dplyr::filter(gene_mut == 1))$Drug_length, na.rm = T) - 1)/365.25*12)
                    tau = floor(min(tau0, tau1, input$RMST_drug, na.rm = T) * 10) / 10
                    verify <- survRM2::rmst2(
                      time = Data_drug_Drug_length_tmp$Drug_length,
                      status = Data_drug_Drug_length_tmp$Drug_length_censor,
                      arm = Data_drug_Drug_length_tmp$gene_mut,
                      tau = tau * 365.25 / 12,
                      alpha = 0.05
                    )

                    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=0)
                    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                      data=Data_drug_Drug_length_tmp, rho=1)
                    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~gene_mut,
                                   data=Data_drug_Drug_length_tmp)

                    mut_gene = paste(oncogenic_genes$cluster[i],
                                     ", top ", i, " gene", sep="")
                    tmp = data.frame(summary(traditional_fit)$table)
                    gene_table$positive_median[i] = round(tmp$median[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_LL[i] = round(tmp$X0.95LCL[2] / 365.25 * 12, digits = 1)
                    gene_table$positive_UL[i] = round(tmp$X0.95UCL[2] / 365.25 * 12, digits = 1)
                    gene_table$negative_median[i] = round(tmp$median[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp$X0.95LCL[1] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp$X0.95UCL[1] / 365.25 * 12, digits = 1)
                    data_tmp = tidy(diff_2, exponentiate=TRUE, conf.int=TRUE)
                    gene_table$diff_median[i] = round(1/data_tmp[[2]], digits = 3)
                    gene_table$diff_LL[i] = round(1/data_tmp[[7]], digits = 3)
                    gene_table$diff_UL[i] = round(1/data_tmp[[6]], digits = 3)
                    hsc[[k]] = surv_curv_drug(traditional_fit, Data_drug_Drug_length_tmp, paste0("Cluster ", oncogenic_genes$cluster[i], ", treated with", paste(input$drug, collapse = ";")),
                                              diff_0, diff_1, diff_2)
                    k = k + 1
                  } else{
                    tmp = summary(traditional_fit)$table
                    legends = paste0(format_p(tmp[[7]] / 365.25 * 12, digits = 1), " (",
                                     format_p(tmp[[8]] / 365.25 * 12, digits = 1), "-",
                                     format_p(tmp[[9]] / 365.25 * 12, digits = 1),")")
                    gene_table$negative_median[i] = round(tmp[[7]] / 365.25 * 12, digits = 1)
                    gene_table$negative_LL[i] = round(tmp[[8]] / 365.25 * 12, digits = 1)
                    gene_table$negative_UL[i] = round(tmp[[9]] / 365.25 * 12, digits = 1)
                  }
                  incProgress(1/nrow(oncogenic_genes))
                }
              })
              Gene_arrange = gene_table$Gene
              # gene_table = gene_table %>% dplyr::mutate(
              #   Gene = factor(gene_table$Gene, levels = Gene_arrange))
              gene_table = gene_table %>% dplyr::mutate(
                Gene = factor(gene_table$Gene))

              # -Inf, Inf を置き換える新しい列を作る
              min_val <- min(gene_table$diff_LL[gene_table$diff_LL != 0], na.rm = TRUE)
              max_val <- max(gene_table$diff_UL[is.finite(gene_table$diff_UL)], na.rm = TRUE)
              gene_table$diff_LL_plot <- ifelse(gene_table$diff_LL != 0, gene_table$diff_LL, min_val)
              gene_table$diff_UL_plot <- ifelse(is.finite(gene_table$diff_UL), gene_table$diff_UL, max_val)
              gene_table$diff_median_plot <- ifelse(gene_table$diff_median >= gene_table$diff_LL_plot &
                                                      gene_table$diff_median <= gene_table$diff_UL_plot, gene_table$diff_median, NA)
              p_mid <-
                gene_table |>
                ggplot(aes(y = fct_rev(Gene))) +
                theme_classic() +  scale_x_log10() +
                geom_point(aes(x=diff_median_plot), shape=15, size=3) +
                geom_linerange(aes(xmin=diff_LL_plot, xmax=diff_UL_plot)) +
                geom_vline(xintercept = 1, linetype="dashed") +
                labs(x=paste0("Hazard ratio with ", paste(input$drug, collapse = ";")),
                     y="Mutation based cluster") +
                coord_cartesian(ylim=c(1,length(Gene_arrange) + 1),
                                xlim=c( min_val ^ 1.1,
                                        max_val ^ 1.1)) +
                annotate("text",
                         x = sqrt(max_val ^ 1.1),
                         y = length(Gene_arrange) + 1,
                         label = "cluster=short survival") +
                annotate("text",
                         x = sqrt(min_val ^ 1.1),
                         y = length(Gene_arrange) + 1,
                         label = "cluster=long survival") +
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
                    Gene = "Cluster",
                    estimate_lab_1 = "Survival, mut (+)",
                    estimate_lab_2 = "Survival, mut (-)",
                    estimate_lab_3 = "Hazard ratio",
                    patients = "Patients"
                  ), gene_table)
              Gene_arrange = gene_table$Gene
              gene_table = gene_table %>% dplyr::mutate(
                Gene = factor(gene_table$Gene, levels = Gene_arrange))

              p_left <- gene_table  %>% ggplot(aes(y = fct_rev(Gene)))
              p_left <- p_left + geom_text(aes(x = 0, label = Gene), hjust = 0,
                                           fontface = ifelse(gene_table$Gene == "Cluster", "bold", "bold.italic"))
              p_left <- p_left + geom_text(aes(x = 1.8, label = estimate_lab_1), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_1 == "Survival, mut (+)", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 3.3, label = estimate_lab_2), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_2 == "Survival, mut (-)", "bold", "plain"))
              p_left <- p_left + geom_text(aes(x = 4.8, label = estimate_lab_3), hjust = 0,
                                           fontface = ifelse(gene_table$estimate_lab_3 == "Hazard ratio", "bold", "plain"))
              p_left <- p_left + theme_void() + coord_cartesian(xlim = c(0, 7.5))
              p_right <- gene_table  |> ggplot() +
                geom_text(aes(x = 0, y = fct_rev(Gene), label = patients), hjust = 0,
                          fontface = ifelse(gene_table$patients == "Patients", "bold", "plain")) +
                theme_void()
              layout <- c(patchwork::area(t = 0, l = 0, b = 30, r = 7.5),
                          patchwork::area(t = 1, l = 7.5, b = 30, r = 12.5),
                          patchwork::area(t = 0, l = 13.5, b = 30, r = 14))
              OUTPUT_DATA$drug_analysis_h_clust = p_left + p_mid + p_right + plot_layout(design = layout)
              for(kk in k:32){
                hsc[[kk]] = ggsurvplot_empty()
              }
              OUTPUT_DATA$drug_analysis_gene_table_cluster = gene_table
            }
            ########################
            OUTPUT_DATA$drug_analysis_hs = hs
            OUTPUT_DATA$drug_analysis_hsc = hsc

            incProgress(1 / 13)
          }

          if(input$drug_volcano != "Yes"){

            Data_drug_RECIST = Data_drug_RECIST %>%
              dplyr::filter(Drug %in% input$drug) %>%
              dplyr::filter(治療ライン %in% input$target_line) %>%
              dplyr::arrange(ID, 治療ライン)


            Drug_summary_RECIST = Data_drug_RECIST %>%
              dplyr::select(
                ID,
                Drug,
                治療ライン,
                最良総合効果,
                終了理由,
                Adverse_effect,
                Overall_AE_grage,
                all_of(names(classification_rules)),
                Drug_length,
                Drug_length_censor,
                Early_AE
              )
            colnames(Drug_summary_RECIST)=
              c("ID", "Drugs", "CTx line", "RECIST", "Reason to finish treatment",
                "Any adverse effect (G3-G5)",
                "Max AE grade",
                names(classification_rules),
                "Treatment time", "Treatment finished or censored", "AE within 1-month, 3-month, or later")
            Drug_summary_RECIST$`Reason to finish treatment`[Drug_summary_RECIST$`Reason to finish treatment` == ""] = "入力なし"
            mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target_tmp$Hugo_Symbol)]
            if(length(mut_gene_)==0){
              mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
              ID_mutation_ = (Data_MAF_target_tmp %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
              Drug_summary_RECIST_line = Drug_summary_RECIST %>% dplyr::mutate(
                separation_value = case_when(
                  ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
                  TRUE ~ paste0(mut_gene_, " mut(-)")
                )
              )
            } else{
              if(is.null(input$gene_group_2)){
                ID_mutation_ = (Data_MAF_target_tmp %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
                Drug_summary_RECIST_line = Drug_summary_RECIST %>% dplyr::mutate(
                  separation_value = case_when(
                    ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
                    TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
                  )
                )
              } else{
                ID_mutation_1 = (Data_MAF_target_tmp %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
                ID_mutation_2 = (Data_MAF_target_tmp %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
                Drug_summary_RECIST_line = Drug_summary_RECIST %>% dplyr::mutate(
                  separation_value = case_when(
                    ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                    !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                  )
                )
              }
            }

            Drug_summary_RECIST_line$`Reason to finish treatment`[is.na(Drug_summary_RECIST_line$`Reason to finish treatment`)] = "Unknown"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "その他理由で中止"] = "Other reason"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "不明"] = "Unknown"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "副作用等で中止"] = "Side effect"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "本人希望により中止"] = "Patient will"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "無効中止"] = "Not effective"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "死亡中止"] = "Death"
            Drug_summary_RECIST_line$`Reason to finish treatment`[Drug_summary_RECIST_line$`Reason to finish treatment` == "計画通り終了"] = "As planned"
            Drug_summary_RECIST_line$`Treatment finished or censored` = as.character(Drug_summary_RECIST_line$`Treatment finished or censored`)
            Drug_summary_RECIST_line$`Treatment finished or censored`[is.na(Drug_summary_RECIST_line$`Treatment finished or censored`)] = "Unknown"
            Drug_summary_RECIST_line$`Treatment finished or censored`[Drug_summary_RECIST_line$`Treatment finished or censored` == "1"] = "Finished"
            Drug_summary_RECIST_line$`Treatment finished or censored`[Drug_summary_RECIST_line$`Treatment finished or censored` == "0"] = "Censored"

            ID_diagnosis_list = Data_case_target %>% dplyr::select(C.CAT調査結果.基本項目.ハッシュID,症例.基本情報.がん種.OncoTree.)
            Drug_summary_RECIST_line$diagnosis = unlist(lapply(list(Drug_summary_RECIST_line$ID), function(x) {
              as.vector(ID_diagnosis_list$症例.基本情報.がん種.OncoTree.[match(x, ID_diagnosis_list$C.CAT調査結果.基本項目.ハッシュID)])}))

            remove_cols <- names(Drug_summary_RECIST_line) %in% names(classification_rules) &
              sapply(Drug_summary_RECIST_line, is.logical)
            Drug_summary_RECIST_line <- Drug_summary_RECIST_line[, !remove_cols, with = FALSE]
            keep_cols_RECIST <- names(Drug_summary_RECIST_line)[names(Drug_summary_RECIST_line) %in% names(classification_rules)]
            Drug_summary_RECIST_line <- Drug_summary_RECIST_line %>%
              mutate(across(all_of(keep_cols_RECIST), as.factor))

            Drug_summary_for_table = Drug_summary_RECIST_line
            Drug_summary_for_table = Drug_summary_for_table %>%
              dplyr::left_join(
                Data_case_target %>% dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = T),
                by = c("ID" = "C.CAT調査結果.基本項目.ハッシュID")
              ) %>%
              dplyr::mutate(
                YoungOld = case_when(
                  YoungOld == "Younger" ~ paste0("Younger than ", input$mid_age + 1),
                  TRUE ~ paste0(input$mid_age + 1, " and older")
                ))

            translation_map <- list(
              "Family cancer history" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
              "Double cancer in a different organ" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown"),
              "Multiple tumor nodules in the organ" = c("あり" = "Yes", "なし" = "No", "不明" = "Unknown")
            )
            # 列名マッピング辞書
            column_map <- c(
              "CTx_lines_before_CGP" = "CTx lines before CGP",
              "症例.背景情報.初回治療前のステージ分類.名称." = "Stage at diagnosis",
              "症例.基本情報.がん種.OncoTree." = "Diagnosis",
              "症例.基本情報.がん種.OncoTree..名称." = "Diagnosis (OncoTree)",
              "症例.基本情報.がん種.OncoTree.LEVEL1." = "Diagnosis (OncoTree level 1)",
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
              "cluster" = "Mutation-based cluster",
              "EP_option" = "Treatment recommended by expert panel",
              "EP_treat" = "Indicated treatment administered",
              "Lymph_met" = "Lymphatic metastasis",
              "Brain_met" = "Brain metastasis",
              "Lung_met" = "Lung metastasis",
              "Bone_met" = "Bone metastasis",
              "Liver_met" = "Liver metastasis",
              "Other_met" = "Other organ metastasis"
            )
            # data.tableの場合
            setDT(Drug_summary_for_table)
            setnames(Drug_summary_for_table, names(column_map), column_map, skip_absent = TRUE)
            # YoungOld列の特別処理
            if ("YoungOld" %in% names(Drug_summary_for_table)) {
              setnames(Drug_summary_for_table, "YoungOld", paste0(input$mid_age + 1, " and older"))
            }
            # 対象列を特定
            target_cols <- intersect(names(translation_map), names(Drug_summary_for_table))
            # 列変換
            for (col in target_cols) {
              Drug_summary_for_table[is.na(get(col)), (col) := "Unknown"]
              Drug_summary_for_table[, (col) := translation_map[[col]][get(col)]]
              Drug_summary_for_table[is.na(get(col)), (col) := "Unknown"]
            }
            OUTPUT_DATA$drug_analysis_Drug_summary_for_table = Drug_summary_for_table

            Drug_summary_RECIST_line = Drug_summary_RECIST_line %>%
              dplyr::left_join(Data_case_age_sex, by = "ID")

            mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
            if(length(mut_gene_)==0){
              mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
              ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
              Drug_summary_RECIST_line = Drug_summary_RECIST_line %>% dplyr::mutate(
                mutation = case_when(
                  ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
                  TRUE ~ paste0(mut_gene_, " mut(-)")
                )
              )
            } else{
              if(is.null(input$gene_group_2)){
                ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
                Drug_summary_RECIST_line = Drug_summary_RECIST_line %>% dplyr::mutate(
                  mutation = case_when(
                    ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
                    TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
                  )
                )
              } else{
                mut_gene_ = input$gene
                ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
                ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
                Drug_summary_RECIST_line = Drug_summary_RECIST_line %>% dplyr::mutate(
                  mutation = case_when(
                    ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                    !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
                    !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
                  )
                )
              }
            }
            Drug_summary_whole = Drug_summary_RECIST_line
            Drug_summary_AE_line = Drug_summary_RECIST_line %>%
              dplyr::filter(!RECIST %in% c("NE", "Unknown")) %>%
              dplyr::filter(`Treatment time`>0 & !is.na(`Treatment time`) & is.finite(`Treatment time`))
            Drug_summary_RECIST_line = Drug_summary_RECIST_line %>%
              dplyr::filter(!RECIST %in% c("NE", "Unknown"))
            Data_forest = Data_drug_RECIST %>%
              dplyr::select(ID, Drug, Drug_length, Drug_length_censor, 最良総合効果, 治療ライン, Adverse_effect)
            colnames(Data_forest) = c("ID", "Drug", "Drug_length", "Drug_length_censor", "最良総合効果", "Lines", "Adverse_effect")
            col_added = 1
            if(!"LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
              Data_survival_tmp = Data_survival %>%
                dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                              YoungOld,
                              Lymph_met,
                              Brain_met,
                              Lung_met,
                              Bone_met,
                              Liver_met,
                              症例.基本情報.性別.名称.,
                              症例.基本情報.がん種.OncoTree.,
                              症例.基本情報.がん種.OncoTree.LEVEL1.,
                              HER2_IHC,
                              MSI_PCR,
                              MMR_IHC,
                              cluster#,
                              #症例.背景情報.喫煙歴有無.名称.,
                              #症例.背景情報.アルコール多飲有無.名称.
                )
            } else {
              Data_survival_tmp = Data_survival %>%
                dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                              YoungOld,
                              Lymph_met,
                              Brain_met,
                              Lung_met,
                              Bone_met,
                              Liver_met,
                              症例.基本情報.性別.名称.,
                              症例.基本情報.がん種.OncoTree.,
                              症例.基本情報.がん種.OncoTree.LEVEL1.,
                              HER2_IHC,
                              MSI_PCR,
                              MMR_IHC,
                              cluster,
                              PD_L1
                              #症例.背景情報.喫煙歴有無.名称.,
                              #症例.背景情報.アルコール多飲有無.名称.
                ) %>% dplyr::mutate(
                  PD_L1 = case_when(
                    PD_L1 == "Positive" ~ "Positive",
                    PD_L1 == "Negative" ~ "Negative",
                    TRUE ~ "Unknown"
                  )
                )
              col_added = col_added + 1
            }
            if(input$HER2 != "No"){
              col_added = col_added + 1
            } else {
              Data_survival_tmp = Data_survival_tmp %>%
                dplyr::select(-HER2_IHC)
            }
            if(input$MSI != "No"){
              col_added = col_added + 1
            } else {
              Data_survival_tmp = Data_survival_tmp %>%
                dplyr::select(-MSI_PCR)
            }
            if(input$MMR != "No"){
              col_added = col_added + 1
            } else {
              Data_survival_tmp = Data_survival_tmp %>%
                dplyr::select(-MMR_IHC)
            }

            Data_survival_tmp = Data_survival_tmp %>%
              dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = TRUE) %>%
              dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% Data_forest$ID)
            colnames_tmp = colnames(Data_survival_tmp)
            colnames_tmp[colnames_tmp == "C.CAT調査結果.基本項目.ハッシュID"] = "ID"
            colnames_tmp[colnames_tmp == "YoungOld"] = "Age"
            colnames_tmp[colnames_tmp == "症例.基本情報.性別.名称."] = "Sex"
            colnames_tmp[colnames_tmp == "症例.基本情報.がん種.OncoTree."] = "Histology"
            colnames_tmp[colnames_tmp == "cluster"] = "Cluster"
            colnames_tmp[colnames_tmp == "症例.基本情報.がん種.OncoTree.LEVEL1."] = "Histology_dummy"
            colnames(Data_survival_tmp) = colnames_tmp

            #Data_survival_tmp$Histology_dummy = paste0(Data_survival_tmp$Histology_dummy,"_NOS")
            Data_survival_tmp$Age[is.na(Data_survival_tmp$Age)] <- "不明"
            Data_forest = left_join(Data_forest, Data_survival_tmp, by="ID")
            Data_forest = Data_forest %>% dplyr::mutate(
              Lines = case_when(
                Lines == "1" ~ "1",
                Lines == "2" ~ "2",
                TRUE ~ "3~"
              ),
              ORR_data = case_when(
                最良総合効果 %in% c("CR", "PR") ~ 1,
                最良総合効果 %in% c("LongSD", "SD", "PD") ~ 0,
                最良総合効果 %in% c("NE", "Unknown") ~ -1,
                TRUE ~ -1
              ),
              DCR_data = case_when(
                最良総合効果 %in% c("CR", "PR", "LongSD", "SD") ~ 1,
                最良総合効果 %in% c("PD") ~ 0,
                最良総合効果 %in% c("NE", "Unknown") ~ -1,
                TRUE ~ -1
              )
            )

            Data_forest$Cluster = factor(Data_forest$Cluster)
            gene_names = unique(c(input$drug_multi_gene))
            for(i in seq_len(length(gene_names))){
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
              Data_forest = left_join(Data_forest, Data_forest_tmp, by="ID")
            }
            Data_forest_all_pts = Data_forest
            Data_forest = Data_forest %>%
              dplyr::filter(Sex != "Unknown") %>%
              dplyr::filter(Age != "不明") #%>%
            #dplyr::filter(Smoking_history != "Unknown" & Alcoholic_history != "Unknown")
            if(input$HER2 != "No"){
              Data_forest = Data_forest %>%
                dplyr::filter(HER2_IHC != "Unknown")
            }
            if(input$MSI != "No"){
              Data_forest = Data_forest %>%
                dplyr::filter(MSI_PCR != "Unknown")
            }
            if(input$MMR != "No"){
              Data_forest = Data_forest %>%
                dplyr::filter(MMR_IHC != "Unknown")
            }
            if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
              Data_forest = Data_forest %>%
                dplyr::filter(PD_L1 != "Unknown")
            }
            Data_forest_AE = Data_forest
            if(input$AE_analysis == "No"){
              Data_forest = Data_forest %>%
                dplyr::select(-Adverse_effect)
              col_added = col_added - 1
            }
            Factor_names = colnames(Data_forest)
            Factor_names = Factor_names[!Factor_names %in% c('ID', 'Drug',
                                                             "Drug_length_censor", "Drug_length",
                                                             "最良総合効果", "ORR_data", "DCR_data")]

            Data_forest_tmp = Data_forest %>% dplyr::filter(!is.na(Drug_length_censor))
            Disease_tmp = unique(sort(Data_forest_tmp$Histology))
            Disease_tmp2 = unique(sort(Data_forest_tmp$Histology_dummy))
            Disease_tmp = c(Disease_tmp[!Disease_tmp %in% Disease_tmp2], Disease_tmp2)
            for(i in length(gene_names):1){
              kk = 18+col_added+i
              if(sum(Data_forest_tmp[[kk]] == "mut(+)" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[-(11+col_added+i)]
              } else if(sum(Data_forest_tmp[[kk]] == "mut(-)" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[-(11+col_added+i)]
              }
            }
            for(i in 1:length(Disease_tmp)){
              Data_forest_tmp2 = Data_forest_tmp %>% dplyr::filter(Histology == Disease_tmp[i])
              if(sum(Data_forest_tmp2$Drug_length_censor == 1,na.rm = T) < 3){
                Data_forest_tmp = Data_forest_tmp %>%
                  dplyr::mutate(Histology = case_when(
                    Histology == Disease_tmp[i] ~ Histology_dummy,
                    TRUE ~ Histology
                  )
                  )
              } else {
                Data_forest_tmp2 = Data_forest_tmp %>% dplyr::filter(Histology != Disease_tmp[i])
                if(sum(Data_forest_tmp2$Drug_length_censor == 1,na.rm = T) < 3){
                  Data_forest_tmp = Data_forest_tmp %>%
                    dplyr::mutate(Histology = case_when(
                      Histology == Disease_tmp[i] ~ Histology_dummy,
                      TRUE ~ Histology
                    )
                    )
                }
              }
            }
            Data_forest_tmp$Histology = factor(Data_forest_tmp$Histology)
            if(length(unique(Data_forest_tmp$Histology))<2){
              Factor_names = Factor_names[!Factor_names %in% c('Histology')]
            } else {
              Data_forest_tmp$Histology = relevel(Data_forest_tmp$Histology,
                                                  ref=names(sort(table(Data_forest_tmp$Histology),decreasing = T))[[1]])
            }
            if(length(unique(Data_forest_tmp$Cluster))<2){
              Factor_names = Factor_names[!Factor_names %in% c('Cluster')]
            } else {
              Data_forest_tmp$Cluster = relevel(Data_forest_tmp$Cluster,
                                                ref=names(sort(table(Data_forest_tmp$Cluster),decreasing = T))[[1]])
            }
            if(sum(Data_forest_tmp$Age == "Younger" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Age')]
            } else if(sum(Data_forest_tmp$Age == "Older" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Age')]
            }
            if(sum(Data_forest_tmp$Sex == "Male" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Sex')]
            } else if(sum(Data_forest_tmp$Sex != "Male" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Sex')]
            }
            if(sum(Data_forest_tmp$Lines == "1" & Data_forest_tmp$Drug_length_censor == 1) < 3 &
               sum(Data_forest_tmp$Lines == "2" & Data_forest_tmp$Drug_length_censor == 1) < 3 &
               sum(Data_forest_tmp$Lines == "3~" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if((sum(Data_forest_tmp$Lines == "1" & Data_forest_tmp$Drug_length_censor == 1) >= 3 &
                       sum(Data_forest_tmp$Lines != "1" & Data_forest_tmp$Drug_length_censor == 1) < 3) |
                      (sum(Data_forest_tmp$Lines == "2" & Data_forest_tmp$Drug_length_censor == 1) >= 3 &
                       sum(Data_forest_tmp$Lines != "2" & Data_forest_tmp$Drug_length_censor == 1) < 3) |
                      (sum(Data_forest_tmp$Lines == "3~" & Data_forest_tmp$Drug_length_censor == 1) >= 3 &
                       sum(Data_forest_tmp$Lines != "3~" & Data_forest_tmp$Drug_length_censor == 1) < 3)){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if(sum(Data_forest_tmp$Lines == "1" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Data_forest_tmp$Lines[Data_forest_tmp$Lines == "1"] = "1~2"
              Data_forest_tmp$Lines[Data_forest_tmp$Lines == "2"] = "1~2"
            } else if(sum(Data_forest_tmp$Lines == "3~" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Data_forest_tmp$Lines[Data_forest_tmp$Lines == "3~"] = "2~"
              Data_forest_tmp$Lines[Data_forest_tmp$Lines == "2"] = "2~"
            }
            if(sum(Data_forest_tmp$Lymph_met == "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lymph_met')]
            } else if(sum(Data_forest_tmp$Lymph_met != "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lymph_met')]
            }
            if(sum(Data_forest_tmp$Lung_met == "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lung_met')]
            } else if(sum(Data_forest_tmp$Lung_met != "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lung_met')]
            }
            if(sum(Data_forest_tmp$Brain_met == "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Brain_met')]
            } else if(sum(Data_forest_tmp$Brain_met != "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Brain_met')]
            }
            if(sum(Data_forest_tmp$Bone_met == "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Bone_met')]
            } else if(sum(Data_forest_tmp$Bone_met != "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Bone_met')]
            }
            if(sum(Data_forest_tmp$Liver_met == "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Liver_met')]
            } else if(sum(Data_forest_tmp$Liver_met != "Yes" & Data_forest_tmp$Drug_length_censor == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Liver_met')]
            }
            if(input$HER2 != "No"){
              if(sum(Data_forest_tmp$HER2_IHC == "Positive" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('HER2_IHC')]
              } else if(sum(Data_forest_tmp$HER2_IHC != "Positive" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('HER2_IHC')]
              }
            }
            if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
              if(sum(Data_forest_tmp$PD_L1 == "Positive" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('PD_L1')]
              } else if(sum(Data_forest_tmp$PD_L1 != "Positive" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('PD_L1')]
              }
            }
            if(input$AE_analysis != "No"){
              if(sum(Data_forest_tmp$Adverse_effect == "AE (+)" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('Adverse_effect')]
              } else if(sum(Data_forest_tmp$Adverse_effect != "AE (+)" & Data_forest_tmp$Drug_length_censor == 1,na.rm = T) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('Adverse_effect')]
              }
            }
            if(input$MSI != "No"){
              if(sum(Data_forest_tmp$MSI_PCR == "Positive" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('MSI_PCR')]
              } else if(sum(Data_forest_tmp$MSI_PCR != "Positive" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('MSI_PCR')]
              }
            }
            if(input$MMR != "No"){
              if(sum(Data_forest_tmp$MMR_IHC == "dMMR" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('MMR_IHC')]
              } else if(sum(Data_forest_tmp$MMR_IHC != "dMMR" & Data_forest_tmp$Drug_length_censor == 1) < 3){
                Factor_names = Factor_names[!Factor_names %in% c('MMR_IHC')]
              }
            }
            if(length(Factor_names) > 1){
              for(i in length(Factor_names):1){
                kk = Factor_names[i]
                if(all(as.character(Data_forest_tmp[,..kk][[1]]) ==
                       as.character(Data_forest_tmp[1,..kk][[1]]))){
                  Factor_names = Factor_names[Factor_names != Factor_names[i]]
                }
              }
            }
            Factor_names = Factor_names[!Factor_names %in% c('Histology_dummy')]
            Factor_names_univariant = Factor_names
            if(length(Factor_names) > 1){
              Factor_names_tmp = Factor_names
              for(i in 1:(length(Factor_names_tmp) - 1)){
                for(j in (i+1):length(Factor_names_tmp)){
                  if(Factor_names_tmp[i] %in% colnames(Data_forest_tmp) &
                     Factor_names_tmp[j] %in% colnames(Data_forest_tmp)){
                    kk = Factor_names_tmp[i]
                    ll = Factor_names_tmp[j]
                    if(abs(cor(as.numeric(as.factor(Data_forest_tmp[,..kk][[1]])),
                               as.numeric(as.factor(Data_forest_tmp[,..ll][[1]])))) > 0.95){
                      Factor_names = Factor_names[Factor_names != ll]
                    }
                  }
                }
              }
            }
            Data_forest_tmp_x = Data_forest_tmp
            Factor_names_x = Factor_names
            Factor_names_x_univariant = Factor_names_univariant


            incProgress(1 / 13)

            Data_forest_tmp_1 = Data_forest  %>%
              dplyr::filter(ORR_data != -1) %>%
              dplyr::select(-ID, -Drug,
                            -Drug_length_censor, -Drug_length,
                            -最良総合効果, -DCR_data)
            Disease_tmp_1 = unique(sort(Data_forest_tmp_1$Histology))
            Disease_tmp_12 = unique(sort(Data_forest_tmp_1$Histology_dummy))
            Disease_tmp_1 = c(Disease_tmp_1[!Disease_tmp_1 %in% Disease_tmp_12], Disease_tmp_12)
            for(i in 1:length(Disease_tmp_1)){
              Data_forest_tmp_12 = Data_forest_tmp_1 %>% dplyr::filter(Histology == Disease_tmp_1[i])
              if(sum(Data_forest_tmp_12$ORR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::mutate(Histology = case_when(
                    Histology == Disease_tmp_1[i] ~ Histology_dummy,
                    TRUE ~ Histology
                  )
                  )
              } else {
                Data_forest_tmp_12 = Data_forest_tmp_1 %>% dplyr::filter(Histology == Disease_tmp_1[i])
                if(sum(Data_forest_tmp_12$ORR_data == 1,na.rm = T) < 3){
                  Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                    dplyr::mutate(Histology = case_when(
                      Histology == Disease_tmp_1[i] ~ Histology_dummy,
                      TRUE ~ Histology
                    )
                    )
                }
              }
            }
            for(i in length(gene_names):1){
              kk = 12+col_added+i
              if(sum(Data_forest_tmp_1[[kk]] == "mut(+)" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1[,-..kk]
              } else if(sum(Data_forest_tmp_1[[kk]] == "mut(-)" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1[,-..kk]
              }
            }
            if(sum(Data_forest_tmp_1$Age == "Younger" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Age)
            } else if(sum(Data_forest_tmp_1$Age == "Older" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Age)
            }
            if(sum(Data_forest_tmp_1$Sex == "Male" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Sex)
            } else if(sum(Data_forest_tmp_1$Sex != "Male" & Data_forest_tmp_1$ORR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Sex)
            }
            if(sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$ORR_data == 1) < 3 &
               sum(Data_forest_tmp_1$Lines == "2" & Data_forest_tmp_1$ORR_data == 1) < 3 &
               sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if((sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$ORR_data == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "1" & Data_forest_tmp_1$ORR_data == 1) < 3) |
                      (sum(Data_forest_tmp_1$Lines == "2" & Data_forest_tmp_1$ORR_data == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "2" & Data_forest_tmp_1$ORR_data == 1) < 3) |
                      (sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$ORR_data == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "3~" & Data_forest_tmp_1$ORR_data == 1) < 3)){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if(sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "1"] = "1~2"
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "2"] = "1~2"
            } else if(sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "3~"] = "2~"
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "2"] = "2~"
            }
            if(sum(Data_forest_tmp_1$Lymph_met == "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lymph_met)
            } else if(sum(Data_forest_tmp_1$Lymph_met != "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lymph_met)
            }
            if(sum(Data_forest_tmp_1$Lung_met == "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lung_met)
            } else if(sum(Data_forest_tmp_1$Lung_met != "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lung_met)
            }
            if(sum(Data_forest_tmp_1$Brain_met == "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Brain_met)
            } else if(sum(Data_forest_tmp_1$Brain_met != "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Brain_met)
            }
            if(sum(Data_forest_tmp_1$Bone_met == "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Bone_met)
            } else if(sum(Data_forest_tmp_1$Bone_met != "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Bone_met)
            }
            if(sum(Data_forest_tmp_1$Liver_met == "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Liver_met)
            } else if(sum(Data_forest_tmp_1$Liver_met != "Yes" & Data_forest_tmp_1$ORR_data == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Liver_met)
            }
            if(input$HER2 != "No"){
              if(sum(Data_forest_tmp_1$HER2_IHC == "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-HER2_IHC)
              } else if(sum(Data_forest_tmp_1$HER2_IHC != "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-HER2_IHC)
              }
            }
            if(input$MSI != "No"){
              if(sum(Data_forest_tmp_1$MSI_PCR == "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MSI_PCR)
              } else if(sum(Data_forest_tmp_1$MSI_PCR != "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MSI_PCR)
              }
            }
            if(input$MMR != "No"){
              if(sum(Data_forest_tmp_1$MMR_IHC == "dMMR" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MMR_IHC)
              } else if(sum(Data_forest_tmp_1$MMR_IHC != "dMMR" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MMR_IHC)
              }
            }
            if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
              if(sum(Data_forest_tmp_1$PD_L1 == "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-PD_L1)
              } else if(sum(Data_forest_tmp_1$PD_L1 != "Positive" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-PD_L1)
              }
            }
            if(input$AE_analysis != "No"){
              if(sum(Data_forest_tmp_1$Adverse_effect == "AE (+)" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-Adverse_effect)
              } else if(sum(Data_forest_tmp_1$Adverse_effect != "AE (+)" & Data_forest_tmp_1$ORR_data == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(Adverse_effect)
              }
            }
            Data_forest_tmp_1$Histology = factor(Data_forest_tmp_1$Histology)
            if(length(unique(Data_forest_tmp_1$Histology))<2){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Histology)
            } else {
              Data_forest_tmp_1$Histology = relevel(Data_forest_tmp_1$Histology,
                                                    ref=names(sort(table(Data_forest_tmp_1$Histology),decreasing = T))[[1]])
            }
            if(length(unique(Data_forest_tmp_1$Cluster))<2){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Cluster)
            } else {
              Data_forest_tmp_1$Cluster = relevel(Data_forest_tmp_1$Cluster,
                                                  ref=names(sort(table(Data_forest_tmp_1$Cluster),decreasing = T))[[1]])
            }
            Data_forest_tmp_1 = Data_forest_tmp_1 %>%
              dplyr::select(-Histology_dummy)
            if(length(colnames(Data_forest_tmp_1)) > 1){
              Factor_names_tmp_1 = colnames(Data_forest_tmp_1)
              for(i in length(Factor_names_tmp_1):1){
                kk = Factor_names_tmp_1[i]
                if(all(as.character(Data_forest_tmp_1[,..kk][[1]]) ==
                       as.character(Data_forest_tmp_1[1,..kk][[1]]))){
                  Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                    dplyr::select(-kk)
                }
              }
            }

            Data_forest_tmp_1_univariant = Data_forest_tmp_1
            if(length(colnames(Data_forest_tmp_1)) > 1){
              Factor_names_tmp_1 = colnames(Data_forest_tmp_1)
              for(i in 1:(length(Factor_names_tmp_1) - 1)){
                for(j in (i+1):length(Factor_names_tmp_1)){
                  if(Factor_names_tmp_1[i] %in% colnames(Data_forest_tmp_1) &
                     Factor_names_tmp_1[j] %in% colnames(Data_forest_tmp_1)){
                    kk = Factor_names_tmp_1[i]
                    ll = Factor_names_tmp_1[j]
                    if(abs(cor(as.numeric(as.factor(Data_forest_tmp_1[,..kk][[1]])),
                               as.numeric(as.factor(Data_forest_tmp_1[,..ll][[1]])))) > 0.95){
                      Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                        dplyr::select(-ll)
                    }
                  }
                }
              }
            }
            Data_forest_tmp_7 = Data_forest_tmp_1
            Data_forest_tmp_7_univariant = Data_forest_tmp_1_univariant

            incProgress(1 / 13)

            Data_forest_tmp_2 = Data_forest  %>%
              dplyr::filter(DCR_data != -1) %>%
              dplyr::select(-ID, -Drug,
                            -Drug_length_censor, -Drug_length,
                            -最良総合効果, -ORR_data, -Cluster)
            Disease_tmp_2 = unique(sort(Data_forest_tmp_2$Histology))
            Disease_tmp_22 = unique(sort(Data_forest_tmp_2$Histology_dummy))
            Disease_tmp_2 = c(Disease_tmp_2[!Disease_tmp_2 %in% Disease_tmp_22], Disease_tmp_22)
            for(i in 1:length(Disease_tmp_2)){
              Data_forest_tmp_22 = Data_forest_tmp_2 %>% dplyr::filter(Histology == Disease_tmp_2[i])
              if(sum(Data_forest_tmp_22$DCR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::mutate(Histology = case_when(
                    Histology == Disease_tmp_2[i] ~ Histology_dummy,
                    TRUE ~ Histology
                  )
                  )
              } else{
                Data_forest_tmp_22 = Data_forest_tmp_2 %>% dplyr::filter(Histology != Disease_tmp_2[i])
                if(sum(Data_forest_tmp_22$DCR_data == 1,na.rm = T) < 3){
                  Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                    dplyr::mutate(Histology = case_when(
                      Histology == Disease_tmp_2[i] ~ Histology_dummy,
                      TRUE ~ Histology
                    )
                    )
                }
              }
            }
            for(i in length(gene_names):1){
              kk = 11+col_added+i
              if(sum(Data_forest_tmp_2[[kk]] == "mut(+)" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2[,-..kk]
              } else if(sum(Data_forest_tmp_2[[kk]] == "mut(-)" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2[,-..kk]
              }
            }
            if(sum(Data_forest_tmp_2$Age == "Younger" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Age)
            } else if(sum(Data_forest_tmp_2$Age == "Older" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Age)
            }
            if(sum(Data_forest_tmp_2$Sex == "Male" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Sex)
            } else if(sum(Data_forest_tmp_2$Sex != "Male" & Data_forest_tmp_2$DCR_data == 1,na.rm = T) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Sex)
            }
            if(sum(Data_forest_tmp_2$Lines == "1" & Data_forest_tmp_2$DCR_data == 1) < 3 &
               sum(Data_forest_tmp_2$Lines == "2" & Data_forest_tmp_2$DCR_data == 1) < 3 &
               sum(Data_forest_tmp_2$Lines == "3~" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if((sum(Data_forest_tmp_2$Lines == "1" & Data_forest_tmp_2$DCR_data == 1) >= 3 &
                       sum(Data_forest_tmp_2$Lines != "1" & Data_forest_tmp_2$DCR_data == 1) < 3) |
                      (sum(Data_forest_tmp_2$Lines == "2" & Data_forest_tmp_2$DCR_data == 1) >= 3 &
                       sum(Data_forest_tmp_2$Lines != "2" & Data_forest_tmp_2$DCR_data == 1) < 3) |
                      (sum(Data_forest_tmp_2$Lines == "3~" & Data_forest_tmp_2$DCR_data == 1) >= 3 &
                       sum(Data_forest_tmp_2$Lines != "3~" & Data_forest_tmp_2$DCR_data == 1) < 3)){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if(sum(Data_forest_tmp_2$Lines == "1" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2$Lines[Data_forest_tmp_2$Lines == "1"] = "1~2"
              Data_forest_tmp_2$Lines[Data_forest_tmp_2$Lines == "2"] = "1~2"
            } else if(sum(Data_forest_tmp_2$Lines == "3~" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2$Lines[Data_forest_tmp_2$Lines == "3~"] = "2~"
              Data_forest_tmp_2$Lines[Data_forest_tmp_2$Lines == "2"] = "2~"
            }
            if(sum(Data_forest_tmp_2$Lymph_met == "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Lymph_met)
            } else if(sum(Data_forest_tmp_2$Lymph_met != "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Lymph_met)
            }
            if(sum(Data_forest_tmp_2$Lung_met == "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Lung_met)
            } else if(sum(Data_forest_tmp_2$Lung_met != "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Lung_met)
            }
            if(sum(Data_forest_tmp_2$Brain_met == "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Brain_met)
            } else if(sum(Data_forest_tmp_2$Brain_met != "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Brain_met)
            }
            if(sum(Data_forest_tmp_2$Bone_met == "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Bone_met)
            } else if(sum(Data_forest_tmp_2$Bone_met != "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Bone_met)
            }
            if(sum(Data_forest_tmp_2$Liver_met == "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Liver_met)
            } else if(sum(Data_forest_tmp_2$Liver_met != "Yes" & Data_forest_tmp_2$DCR_data == 1) < 3){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Liver_met)
            }
            if(input$HER2 != "No"){
              if(sum(Data_forest_tmp_2$HER2_IHC == "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-HER2_IHC)
              } else if(sum(Data_forest_tmp_2$HER2_IHC != "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-HER2_IHC)
              }
            }
            if(input$MSI != "No"){
              if(sum(Data_forest_tmp_2$MSI_PCR == "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-MSI_PCR)
              } else if(sum(Data_forest_tmp_2$MSI_PCR != "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-MSI_PCR)
              }
            }
            if(input$MMR != "No"){
              if(sum(Data_forest_tmp_2$MMR_IHC == "dMMR" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-MMR_IHC)
              } else if(sum(Data_forest_tmp_2$MMR_IHC != "dMMR" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-MMR_IHC)
              }
            }
            if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
              if(sum(Data_forest_tmp_2$PD_L1 == "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-PD_L1)
              } else if(sum(Data_forest_tmp_2$PD_L1 != "Positive" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-PD_L1)
              }
            }
            if(input$AE_analysis != "No"){
              if(sum(Data_forest_tmp_2$Adverse_effect == "AE (+)" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(-Adverse_effect)
              } else if(sum(Data_forest_tmp_2$Adverse_effect != "AE (+)" & Data_forest_tmp_2$DCR_data == 1) < 3){
                Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                  dplyr::select(Adverse_effect)
              }
            }
            Data_forest_tmp_2$Histology = factor(Data_forest_tmp_2$Histology)
            if(length(unique(Data_forest_tmp_2$Histology))<2){
              Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                dplyr::select(-Histology)
            }
            Data_forest_tmp_2 = Data_forest_tmp_2 %>%
              dplyr::select(-Histology_dummy)
            if(length(colnames(Data_forest_tmp_2)) > 1){
              Factor_names_tmp_2 = colnames(Data_forest_tmp_2)
              for(i in length(Factor_names_tmp_2):1){
                kk = Factor_names_tmp_2[i]
                if(all(as.character(Data_forest_tmp_2[,..kk][[1]]) ==
                       as.character(Data_forest_tmp_2[1,..kk][[1]]))){
                  Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                    dplyr::select(-kk)
                }
              }
            }
            Data_forest_tmp_2_univariant = Data_forest_tmp_2
            if(length(colnames(Data_forest_tmp_2)) > 1){
              Factor_names_tmp_2 = colnames(Data_forest_tmp_2)
              for(i in 1:(length(Factor_names_tmp_2) - 1)){
                for(j in (i+1):length(Factor_names_tmp_2)){
                  if(Factor_names_tmp_2[i] %in% colnames(Data_forest_tmp_2) &
                     Factor_names_tmp_2[j] %in% colnames(Data_forest_tmp_2)){
                    kk = Factor_names_tmp_2[i]
                    ll = Factor_names_tmp_2[j]
                    if(abs(cor(as.numeric(as.factor(Data_forest_tmp_2[,..kk][[1]])),
                               as.numeric(as.factor(Data_forest_tmp_2[,..ll][[1]])))) > 0.95){
                      Data_forest_tmp_2 = Data_forest_tmp_2 %>%
                        dplyr::select(-ll)
                    }
                  }
                }
              }
            }

            Data_forest_tmp_8 = Data_forest_tmp_2
            Data_forest_tmp_8_univariant = Data_forest_tmp_2_univariant

            incProgress(1 / 13)

            if(input$AE_analysis != "No"){
              col_added = col_added - 1
            }

            Data_forest_tmp_1 = Data_forest_AE %>%
              dplyr::filter(Drug_length>0 & !is.na(Drug_length) & is.finite(Drug_length)) %>%
              dplyr::filter(ORR_data != -1) %>%
              dplyr::mutate(Adverse_effect = case_when(
                Adverse_effect == "AE (+)" ~ 1,
                TRUE ~ 0
              )) %>%
              dplyr::select(-ID, -Drug,
                            -Drug_length_censor,
                            -最良総合効果, -DCR_data)
            Disease_tmp_1 = unique(sort(Data_forest_tmp_1$Histology))
            Disease_tmp_12 = unique(sort(Data_forest_tmp_1$Histology_dummy))
            Disease_tmp_1 = c(Disease_tmp_1[!Disease_tmp_1 %in% Disease_tmp_12], Disease_tmp_12)
            for(i in 1:length(Disease_tmp_1)){
              Data_forest_tmp_12 = Data_forest_tmp_1 %>% dplyr::filter(Histology == Disease_tmp_1[i])
              if(sum(Data_forest_tmp_12$ORR_data == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::mutate(Histology = case_when(
                    Histology == Disease_tmp_1[i] ~ Histology_dummy,
                    TRUE ~ Histology
                  )
                  )
              } else {
                Data_forest_tmp_12 = Data_forest_tmp_1 %>% dplyr::filter(Histology == Disease_tmp_1[i])
                if(sum(Data_forest_tmp_12$ORR_data == 1,na.rm = T) < 3){
                  Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                    dplyr::mutate(Histology = case_when(
                      Histology == Disease_tmp_1[i] ~ Histology_dummy,
                      TRUE ~ Histology
                    )
                    )
                }
              }
            }
            for(i in length(gene_names):1){
              kk = 14+col_added+i
              if(sum(Data_forest_tmp_1[[kk]] == "mut(+)" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1[,-..kk]
              } else if(sum(Data_forest_tmp_1[[kk]] == "mut(-)" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1[,-..kk]
              }
            }
            if(sum(Data_forest_tmp_1$Age == "Younger" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Age)
            } else if(sum(Data_forest_tmp_1$Age == "Older" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Age)
            }
            if(sum(Data_forest_tmp_1$Sex == "Male" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Sex)
            } else if(sum(Data_forest_tmp_1$Sex != "Male" & Data_forest_tmp_1$Adverse_effect == 1,na.rm = T) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Sex)
            }
            if(sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$Adverse_effect == 1) < 3 &
               sum(Data_forest_tmp_1$Lines == "2" & Data_forest_tmp_1$Adverse_effect == 1) < 3 &
               sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if((sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$Adverse_effect == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "1" & Data_forest_tmp_1$Adverse_effect == 1) < 3) |
                      (sum(Data_forest_tmp_1$Lines == "2" & Data_forest_tmp_1$Adverse_effect == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "2" & Data_forest_tmp_1$Adverse_effect == 1) < 3) |
                      (sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$Adverse_effect == 1) >= 3 &
                       sum(Data_forest_tmp_1$Lines != "3~" & Data_forest_tmp_1$Adverse_effect == 1) < 3)){
              Factor_names = Factor_names[!Factor_names %in% c('Lines')]
            } else if(sum(Data_forest_tmp_1$Lines == "1" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "1"] = "1~2"
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "2"] = "1~2"
            } else if(sum(Data_forest_tmp_1$Lines == "3~" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "3~"] = "2~"
              Data_forest_tmp_1$Lines[Data_forest_tmp_1$Lines == "2"] = "2~"
            }
            if(sum(Data_forest_tmp_1$ORR_data == 1 & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-ORR_data)
            } else if(sum(Data_forest_tmp_1$ORR_data != 1 & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-ORR_data)
            }
            if(sum(Data_forest_tmp_1$Lymph_met == "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lymph_met)
            } else if(sum(Data_forest_tmp_1$Lymph_met != "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lymph_met)
            }
            if(sum(Data_forest_tmp_1$Lung_met == "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lung_met)
            } else if(sum(Data_forest_tmp_1$Lung_met != "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Lung_met)
            }
            if(sum(Data_forest_tmp_1$Brain_met == "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Brain_met)
            } else if(sum(Data_forest_tmp_1$Brain_met != "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Brain_met)
            }
            if(sum(Data_forest_tmp_1$Bone_met == "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Bone_met)
            } else if(sum(Data_forest_tmp_1$Bone_met != "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Bone_met)
            }
            if(sum(Data_forest_tmp_1$Liver_met == "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Liver_met)
            } else if(sum(Data_forest_tmp_1$Liver_met != "Yes" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Liver_met)
            }
            if(input$HER2 != "No"){
              if(sum(Data_forest_tmp_1$HER2_IHC == "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-HER2_IHC)
              } else if(sum(Data_forest_tmp_1$HER2_IHC != "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-HER2_IHC)
              }
            }
            if(input$MSI != "No"){
              if(sum(Data_forest_tmp_1$MSI_PCR == "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MSI_PCR)
              } else if(sum(Data_forest_tmp_1$MSI_PCR != "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MSI_PCR)
              }
            }
            if(input$MMR != "No"){
              if(sum(Data_forest_tmp_1$MMR_IHC == "dMMR" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MMR_IHC)
              } else if(sum(Data_forest_tmp_1$MMR_IHC != "dMMR" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-MMR_IHC)
              }
            }
            if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
              if(sum(Data_forest_tmp_1$PD_L1 == "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-PD_L1)
              } else if(sum(Data_forest_tmp_1$PD_L1 != "Positive" & Data_forest_tmp_1$Adverse_effect == 1) < 3){
                Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                  dplyr::select(-PD_L1)
              }
            }
            Data_forest_tmp_1$Histology = factor(Data_forest_tmp_1$Histology)
            if(length(unique(Data_forest_tmp_1$Histology))<2){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Histology)
            } else {
              Data_forest_tmp_1$Histology = relevel(Data_forest_tmp_1$Histology,
                                                    ref=names(sort(table(Data_forest_tmp_1$Histology),decreasing = T))[[1]])
            }
            if(length(unique(Data_forest_tmp_1$Cluster))<2){
              Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                dplyr::select(-Cluster)
            } else {
              Data_forest_tmp_1$Cluster = relevel(Data_forest_tmp_1$Cluster,
                                                  ref=names(sort(table(Data_forest_tmp_1$Cluster),decreasing = T))[[1]])
            }
            Data_forest_tmp_1 = Data_forest_tmp_1 %>%
              dplyr::select(-Histology_dummy)
            if(length(colnames(Data_forest_tmp_1)) > 1){
              Factor_names_tmp_1 = colnames(Data_forest_tmp_1)
              for(i in length(Factor_names_tmp_1):1){
                kk = Factor_names_tmp_1[i]
                if(all(as.character(Data_forest_tmp_1[,..kk][[1]]) ==
                       as.character(Data_forest_tmp_1[1,..kk][[1]]))){
                  Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                    dplyr::select(-kk)
                }
              }
            }

            Data_forest_tmp_1_univariant = Data_forest_tmp_1
            if(length(colnames(Data_forest_tmp_1)) > 1){
              Factor_names_tmp_1 = colnames(Data_forest_tmp_1)
              for(i in 1:(length(Factor_names_tmp_1) - 1)){
                for(j in (i+1):length(Factor_names_tmp_1)){
                  if(Factor_names_tmp_1[i] %in% colnames(Data_forest_tmp_1) &
                     Factor_names_tmp_1[j] %in% colnames(Data_forest_tmp_1)){
                    kk = Factor_names_tmp_1[i]
                    ll = Factor_names_tmp_1[j]
                    if(abs(cor(as.numeric(as.factor(Data_forest_tmp_1[,..kk][[1]])),
                               as.numeric(as.factor(Data_forest_tmp_1[,..ll][[1]])))) > 0.95){
                      Data_forest_tmp_1 = Data_forest_tmp_1 %>%
                        dplyr::select(-ll)
                    }
                  }
                }
              }
            }
            Data_forest_tmp_9 = Data_forest_tmp_1
            Data_forest_tmp_9_univariant = Data_forest_tmp_1_univariant
          }

          incProgress(1 / 13)
         Data_drug_survival_1_save = Data_drug_original %>%
            dplyr::arrange(ID) %>%
            dplyr::filter(!is.na(TTD) & !is.na(TTE)) %>%
            dplyr::filter(TTD > 0 & TTE > 0 & TTD > TTE) %>%
            dplyr::mutate(Regimen_whole = case_when(
              治療ライン %in% input$target_line &
                Drug %in% input$drug ~ "Selected regimens",
              治療ライン %in% input$target_line ~ "Other regimens",
              TRUE ~ "Other"
            ))  %>%
           dplyr::filter(Regimen_whole != "Other")
          OUTPUT_DATA$Data_drug_Data_case_target = Data_case_target
          OUTPUT_DATA$Data_drug_Data_drug_survival_1_save = Data_drug_survival_1_save
          OUTPUT_DATA$Data_drug_Data_drug_original = Data_drug_original
          OUTPUT_DATA$Data_drug_Data_drug_TTF = Data_drug_TTF
          OUTPUT_DATA$Data_drug_Data_forest_all_pts = Data_forest_all_pts
          OUTPUT_DATA$Data_drug_Data_MAF_target = Data_MAF_target
          OUTPUT_DATA$Data_drug_Data_MAF_target_tmp = Data_MAF_target_tmp
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_7 = Data_forest_tmp_7
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_7_univariant = Data_forest_tmp_7_univariant
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_8 = Data_forest_tmp_8
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_8_univariant = Data_forest_tmp_8_univariant
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_9 = Data_forest_tmp_9
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_9_univariant = Data_forest_tmp_9_univariant
          OUTPUT_DATA$drug_analysis_gene_names = gene_names
          OUTPUT_DATA$drug_analysis_Factor_names_x = Factor_names_x
          OUTPUT_DATA$drug_analysis_Factor_names_x_univariant = Factor_names_x_univariant
          OUTPUT_DATA$drug_analysis_Data_forest_tmp_x = Data_forest_tmp_x
          OUTPUT_DATA$drug_analysis_Drug_summary_whole = Drug_summary_whole
          OUTPUT_DATA$drug_analysis_Drug_summary_line = Drug_summary_line
          OUTPUT_DATA$drug_analysis_Drug_summary_RECIST_line = Drug_summary_RECIST_line
          OUTPUT_DATA$drug_analysis_Drug_summary_AE_line = Drug_summary_AE_line
          OUTPUT_DATA$drug_analysis_keep_cols = keep_cols
          OUTPUT_DATA$drug_analysis_keep_cols_RECIST = keep_cols_RECIST

        }
      })

    }
  })
  rm(analysis_env)
  gc()
}

output$figure_drug_1 = renderPlot({
  req(input$ToT_var_1)  # NULL のときはスキップ
  # if(input$ToT_var_1 == "AE_matrix"){
  #   OUTPUT_DATA$drug_analysis_gb_AE_heatmap
  # } else
  if(input$ToT_var_1 == "previous"){
    OUTPUT_DATA$drug_analysis_gb_previous
  } else if(input$ToT_var_1 == "selected_line"){
    OUTPUT_DATA$drug_analysis_gb_all_selected_line
  } else if(input$ToT_var_1 == "mutation"){
    OUTPUT_DATA$drug_analysis_gb_mutation
  } else if(input$ToT_var_1 == "regimen"){
    OUTPUT_DATA$drug_analysis_gb_regimen
  } else if(input$ToT_var_1 == "regimen_detail"){
    OUTPUT_DATA$drug_analysis_gb_regimen_detail
  } else if(input$ToT_var_1 == "mutation_regimen"){
    OUTPUT_DATA$drug_analysis_gb_mutation_regimen
  } else if(input$ToT_var_1 == "mutation_detail"){
    OUTPUT_DATA$drug_analysis_gb_mutation_detail
  } else if(input$ToT_var_1 == "lymph"){
    OUTPUT_DATA$drug_analysis_gb_lymph
  } else if(input$ToT_var_1 == "gc_diagnosis"){
    OUTPUT_DATA$drug_analysis_gc_diagnosis
  } else if(input$ToT_var_1 == "gc_diagnosis_all"){
    OUTPUT_DATA$drug_analysis_gc_diagnosis_all
  } else if(input$ToT_var_1 == "gc_regimen"){
    OUTPUT_DATA$drug_analysis_gc_regimen
  } else if(input$ToT_var_1 == "AE"){
    OUTPUT_DATA$drug_analysis_gb_AE
  } else if(input$ToT_var_1 == "Early_AE"){
    OUTPUT_DATA$drug_analysis_gb_Early_AE
  } else if(input$ToT_var_1 == "gc_AE"){
    OUTPUT_DATA$drug_analysis_gc_AE
  } else if(input$ToT_var_1 == "gc_EAE"){
    OUTPUT_DATA$drug_analysis_gc_EAE
  }
})


output$select_table_var_drug_2 = renderUI({
  req(OUTPUT_DATA$Data_drug_Data_case_target,
      OUTPUT_DATA$drug_analysis_keep_cols,
      OUTPUT_DATA$drug_analysis_keep_cols_RECIST)
  choices_table = c("All" = "All",
                    "Mutation pattern" = "mutations",
                    "Diagnosis" = "diagnosis",
                    "Mutated gene" = "gene",
                    "Treatment line" = "line",
                    structure("age", names = paste0("Age <=", input$mid_age, " or older")),
                    "Sex" = "sex",
                    "RECIST" = "RECIST",
                    "Adverse effect" = "Adverse_effect",
                    "PD-L1 (lung only)" = "PD-L1",
                    "AE within 1-month, 3-month, or later" = "AE within 1-month, 3-month, or later",
                    structure(union(OUTPUT_DATA$drug_analysis_keep_cols, OUTPUT_DATA$drug_analysis_keep_cols_RECIST),
                              names = union(OUTPUT_DATA$drug_analysis_keep_cols, OUTPUT_DATA$drug_analysis_keep_cols_RECIST))
  )
  if(!"LUNG" %in% OUTPUT_DATA$Data_drug_Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. | input$PD_L1 == "No"){
    choices_table = choices_table[choices_table != "PD-L1"]
  }
  pickerInput(inputId = "table_var_drug_2",
              label = "Grouped by",
              choices = choices_table,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = FALSE)
})


# Helper function to get dataset configuration
table_drug_tot_get_dataset <- function(dataset_type) {
  switch(dataset_type,
         "all" = list(
           data = OUTPUT_DATA$drug_analysis_Drug_summary_whole,
           keep_cols = OUTPUT_DATA$drug_analysis_keep_cols,
           req_items = list(OUTPUT_DATA$drug_analysis_Drug_summary_whole, OUTPUT_DATA$drug_analysis_keep_cols)
         ),
         "time" = list(
           data = OUTPUT_DATA$drug_analysis_Drug_summary_line,
           keep_cols = OUTPUT_DATA$drug_analysis_keep_cols,
           req_items = list(OUTPUT_DATA$drug_analysis_Drug_summary_line, OUTPUT_DATA$drug_analysis_keep_cols)
         ),
         "RECIST" = list(
           data = OUTPUT_DATA$drug_analysis_Drug_summary_RECIST_line,
           keep_cols = OUTPUT_DATA$drug_analysis_keep_cols_RECIST,
           req_items = list(OUTPUT_DATA$drug_analysis_Drug_summary_RECIST_line, OUTPUT_DATA$drug_analysis_keep_cols_RECIST)
         ),
         "AE" = list(
           data = OUTPUT_DATA$drug_analysis_Drug_summary_AE_line,
           keep_cols = union(OUTPUT_DATA$drug_analysis_keep_cols, OUTPUT_DATA$drug_analysis_keep_cols_RECIST),
           req_items = list(OUTPUT_DATA$drug_analysis_Drug_summary_AE_line, OUTPUT_DATA$drug_analysis_keep_cols, OUTPUT_DATA$drug_analysis_keep_cols_RECIST)
         )
  )
}

# Helper function to get table configuration based on variable type
table_drug_tot_get_table_config <- function(table_var, mid_age = NULL) {
  switch(table_var,
         "All" = list(
           by_var = NULL,
           remove_cols = c("ID", "separation_value"),
           caption = "**Palliative CTx with treatment duration information, selected lines**"
         ),
         "AE within 1-month, 3-month, or later" = list(
           by_var = "AE within 1-month, 3-month, or later",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by adverse effect timing**"
         ),
         "mutations" = list(
           by_var = "separation_value",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by mutation pattern**"
         ),
         "RECIST" = list(
           by_var = "RECIST",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by treatment effect**"
         ),
         "age" = list(
           by_var = paste0("Age <=", mid_age, " or older"),
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by age**"
         ),
         "sex" = list(
           by_var = "Sex",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by sex**"
         ),
         "Adverse_effect" = list(
           by_var = "Any adverse effect (G3-G5)",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by adverse effect**"
         ),
         "PD-L1" = list(
           by_var = "PD-L1 (lung only)",
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by PD-L1 IHC**"
         ),
         "diagnosis" = list(
           by_var = "diagnosis",
           remove_cols = c("ID", "separation_value"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by diagnosis**"
         ),
         "gene" = list(
           by_var = "mutation",
           remove_cols = c("ID", "separation_value"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by mutation**"
         ),
         "line" = list(
           by_var = "CTx line",
           remove_cols = c("ID", "separation_value"),
           caption = "**Palliative CTx with treatment duration information, selected lines, by treatment line**"
         ),
         # Default for variables in keep_cols_select
         list(
           by_var = table_var,
           remove_cols = c("ID"),
           caption = "**Palliative CTx with treatment duration information, selected lines**"
         )
  )
}

# Helper function to create standardized tbl_summary
table_drug_tot_create_summary <- function(data, config) {
  # Prepare data
  processed_data <- data %>%
    dplyr::select(-all_of(config$remove_cols)) %>%
    select(where(~ any(!is.na(.))))

  # Create base summary
  if (is.null(config$by_var)) {
    summary_result <- processed_data %>%
      tbl_summary(
        statistic = list(all_continuous() ~ c("{N_nonmiss}",
                                              "{mean} ({sd})",
                                              "{median} ({p25}, {p75})",
                                              "{min}, {max}"),
                         all_categorical() ~ "{n} ({p}%)"),
        type = list(all_continuous() ~ "continuous2",
                    all_dichotomous() ~ "categorical"),
        digits = list(
          all_continuous() ~ list(N_nonmiss = 0, mean = 1, sd = 1, median = 1, p25 = 1, p75 = 1, min = 1, max = 1),
          all_categorical() ~ c(0, 1)
        ),
        missing = "ifany")
  } else {
    summary_result <- processed_data %>%
      tbl_summary(
        by = !!sym(config$by_var),
        statistic = list(all_continuous() ~ c("{N_nonmiss}",
                                              "{mean} ({sd})",
                                              "{median} ({p25}, {p75})",
                                              "{min}, {max}"),
                         all_categorical() ~ "{n} ({p}%)"),
        type = list(all_continuous() ~ "continuous2",
                    all_dichotomous() ~ "categorical"),
        digits = list(
          all_continuous() ~ list(N_nonmiss = 0, mean = 1, sd = 1, median = 1, p25 = 1, p75 = 1, min = 1, max = 1),
          all_categorical() ~ c(0, 1)
        ),
        missing = "ifany")
  }

  # Apply standard modifications
  summary_result %>%
    modify_table_body(
      ~ .x %>%
        mutate(across(starts_with("stat_"), ~ str_replace_all(.x, "\\(0\\.0%\\)", "(<0.1%)")))
    ) %>%
    modify_caption(config$caption) %>%
    add_n() %>%
    bold_labels() %>%
    as_gt()
}

# Main refactored render function
output$table_drug_ToT_1 <- render_gt({
  req(input$table_var_drug_2,
      input$table_var_drug_dataset,
      OUTPUT_DATA$drug_analysis_Drug_summary_whole,
      OUTPUT_DATA$drug_analysis_Drug_summary_line,
      OUTPUT_DATA$drug_analysis_Drug_summary_RECIST_line,
      OUTPUT_DATA$drug_analysis_Drug_summary_AE_line,
      OUTPUT_DATA$drug_analysis_keep_cols,
      OUTPUT_DATA$drug_analysis_keep_cols_RECIST
  )

  # Get dataset configuration
  dataset_config <- table_drug_tot_get_dataset(input$table_var_drug_dataset)

  # Apply additional req() for specific dataset
  do.call(req, dataset_config$req_items)

  # Get data and columns
  Drug_summary <- dataset_config$data
  keep_cols_select <- dataset_config$keep_cols

  # Determine table configuration
  if (input$table_var_drug_2 %in% keep_cols_select) {
    # Handle dynamic variables from keep_cols_select
    table_config <- list(
      by_var = input$table_var_drug_2,
      remove_cols = c("ID"),
      caption = "**Palliative CTx with treatment duration information, selected lines, by mutation pattern**"
    )
  } else {
    # Handle predefined variables
    table_config <- table_drug_tot_get_table_config(input$table_var_drug_2, input$mid_age)
  }

  # Create and return summary table
  table_drug_tot_create_summary(Drug_summary, table_config)
})


output$table_volcano_AE = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$Data_drug_df_volcano_AE)
  create_datatable_with_confirm(OUTPUT_DATA$Data_drug_df_volcano_AE, "FELIS downloaded raw data in Volcano plot for adverse effect tab")
})

output$table_outcome_3 = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$Data_drug_Data_case_target,
      OUTPUT_DATA$Data_drug_Data_MAF_target,
      OUTPUT_DATA$Data_drug_Data_forest_all_pts,
      OUTPUT_DATA$Data_drug_Data_MAF_target_tmp)
  req(input$drug_ORR_table_var)
  Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
  Data_MAF_target = OUTPUT_DATA$Data_drug_Data_MAF_target
  Data_MAF_target_tmp = OUTPUT_DATA$Data_drug_Data_MAF_target_tmp
  Data_forest_all_pts = OUTPUT_DATA$Data_drug_Data_forest_all_pts
  if(input$drug_ORR_table_var == "Mutation"){
    mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
    if(length(mut_gene_)==0){
      mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
      ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
      Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
        mutation = case_when(
          ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
          TRUE ~ paste0(mut_gene_, " mut(-)")
        )
      )
    } else{
      if(is.null(input$gene_group_2)){
        ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
        Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
          mutation = case_when(
            ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
            TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
          )
        )
      } else{
        ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
        ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
        Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
          mutation = case_when(
            ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
            ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
            !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
            !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
          )
        )
      }
    }
    column_name = input$drug_ORR_table_var
  } else if(input$drug_ORR_table_var == "Histology"){
    Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
      mutation = Histology
    )
    column_name = input$drug_ORR_table_var
  } else if(input$drug_ORR_table_var == "Cluster"){
    Data_forest_all_pts$mutation = ""
    for(i in sort(unique(Data_case_target$cluster))){
      ID_mutation = (Data_case_target %>% dplyr::filter(cluster == i))$C.CAT調査結果.基本項目.ハッシュID
      if(nrow(Data_forest_all_pts[Data_forest_all_pts$ID %in% ID_mutation,])>0){
        Data_forest_all_pts[Data_forest_all_pts$ID %in% ID_mutation,]$mutation = paste0("Cluster ", i)
      }
    }
    column_name = input$drug_ORR_table_var
  } else if(input$drug_ORR_table_var == "About the gene classified in detail in Settings"){
    if(!input$special_gene == ""){
      Data_MAF_target_tmp = Data_MAF_target_tmp %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
        dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)

      ID_special_gene_ = (Data_MAF_target_tmp %>%
                            dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_"))))$Tumor_Sample_Barcode
      ID_special_gene_mutation_1_ = (Data_MAF_target_tmp %>%
                                       dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                       amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
      ID_special_gene_mutation_2_ = (Data_MAF_target_tmp %>%
                                       dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                       amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
      Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
        separation_value = case_when(
          ID %in% ID_special_gene_mutation_1_ & ID %in% ID_special_gene_mutation_2_ ~
            paste0(input$special_gene_mutation_1_name, " & ", input$special_gene_mutation_2_name, " in ", input$special_gene),
          ID %in% ID_special_gene_mutation_1_ ~ paste0(input$special_gene_mutation_1_name, " in ", input$special_gene),
          ID %in% ID_special_gene_mutation_2_ ~ paste0(input$special_gene_mutation_2_name, " in ", input$special_gene),
          ID %in% ID_special_gene_ ~ paste0("Other mut in ", input$special_gene),
          TRUE ~ paste0("No mut in ", input$special_gene)
        )
      )
    } else{
      Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
        separation_value = "No mutation pattern")
    }
    Data_forest_all_pts = Data_forest_all_pts %>% dplyr::mutate(
      mutation = separation_value
    )
    column_name = input$special_gene
  }

  Table_outcome_mutation_total = NULL
  drug_list = list(input$drug, input$drug_group_1, input$drug_group_2, input$drug_group_3, input$drug_group_4)
  for(drug_count in 1:5){
    drug_name_2 = case_when(
      drug_count == 1 ~ paste(input$drug, collapse = ";"),
      drug_count == 2 ~ input$drug_group_1_name %||% NA_character_,
      drug_count == 3 ~ input$drug_group_2_name %||% NA_character_,
      drug_count == 4 ~ input$drug_group_3_name %||% NA_character_,
      drug_count == 5 ~ input$drug_group_4_name %||% NA_character_
    )
    Data_summary_mutation = Data_forest_all_pts %>% dplyr::filter(!最良総合効果 %in% c("NE", "Unknown")) %>% dplyr::filter(Drug %in% drug_list[[drug_count]])  %>% dplyr::arrange(-ORR_data, -DCR_data)
    if(nrow(Data_summary_mutation)>0){
      Table_outcome_mutation = data.frame(c("Total", as.character((Data_summary_mutation %>% group_by(mutation) %>% tally())$mutation)))
      colnames(Table_outcome_mutation) = column_name
      Table_outcome_mutation$Drug = drug_name_2
      Table_outcome_mutation$N = 0
      Table_outcome_mutation$OR = 0
      Table_outcome_mutation$ORR = 0
      Table_outcome_mutation$ORR_95CI_LL = 0
      Table_outcome_mutation$ORR_95CI_UL = 0
      Table_outcome_mutation$DC = 0
      Table_outcome_mutation$DCR = 0
      Table_outcome_mutation$DCR_95CI_LL = 0
      Table_outcome_mutation$DCR_95CI_UL = 0
      Table_outcome_mutation$OR_OddsRatio = 1
      Table_outcome_mutation$OR_Odds_pValue = 1
      Table_outcome_mutation$DC_OddsRatio = 1
      Table_outcome_mutation$DC_Odds_pValue = 1
      tmp2 = Data_summary_mutation %>% group_by(mutation) %>% tally()
      Table_outcome_mutation$N[1] = sum(tmp2$n)
      if(nrow(Table_outcome_mutation)>1){
        for(i in 1:(nrow(Table_outcome_mutation)-1)){
          Table_outcome_mutation$N[i+1] = tmp2$n[i]
        }

        tmp2 = Data_summary_mutation %>% group_by(mutation, ORR_data) %>% tally()
        Table_outcome_mutation$OR[1] = sum((tmp2 %>% dplyr::filter(ORR_data == 1))$n)
        for(i in 1:(nrow(Table_outcome_mutation)-1)){
          if(nrow(tmp2 %>% dplyr::filter(ORR_data == 1 & mutation == Table_outcome_mutation[[column_name]][i+1]))>0){
            Table_outcome_mutation$OR[i+1] = (tmp2 %>% dplyr::filter(ORR_data == 1 & mutation == Table_outcome_mutation[[column_name]][i+1]))$n
          }
        }
        tmp2 = Data_summary_mutation %>% group_by(mutation, DCR_data) %>% tally()
        Table_outcome_mutation$DC[1] = sum((tmp2 %>% dplyr::filter(DCR_data == 1))$n)
        for(i in 1:(nrow(Table_outcome_mutation)-1)){
          if(nrow(tmp2 %>% dplyr::filter(DCR_data == 1 & mutation == Table_outcome_mutation[[column_name]][i+1]))>0){
            Table_outcome_mutation$DC[i+1] = (tmp2 %>% dplyr::filter(DCR_data == 1 & mutation == Table_outcome_mutation[[column_name]][i+1]))$n
          }
        }

        for(i in 1:nrow(Table_outcome_mutation)){
          Table_outcome_mutation$ORR[i] = format_p(Table_outcome_mutation$OR[i] / Table_outcome_mutation$N[i], digits = 3)
          Table_outcome_mutation$DCR[i] = format_p(Table_outcome_mutation$DC[i] / Table_outcome_mutation$N[i], digits = 3)
          Table_outcome_mutation$ORR_95CI_LL[i] = format_p(exactci(Table_outcome_mutation$OR[i], Table_outcome_mutation$N[i], conf.level=0.95)$conf.int[1], digits = 3)
          Table_outcome_mutation$ORR_95CI_UL[i] = format_p(exactci(Table_outcome_mutation$OR[i], Table_outcome_mutation$N[i], conf.level=0.95)$conf.int[2], digits = 3)
          Table_outcome_mutation$DCR_95CI_LL[i] = format_p(exactci(Table_outcome_mutation$DC[i], Table_outcome_mutation$N[i], conf.level=0.95)$conf.int[1], digits = 3)
          Table_outcome_mutation$DCR_95CI_UL[i] = format_p(exactci(Table_outcome_mutation$DC[i], Table_outcome_mutation$N[i], conf.level=0.95)$conf.int[2], digits = 3)
          test_result <- fisher.test(matrix(c(Table_outcome_mutation$OR[i],
                                              Table_outcome_mutation$N[i] - Table_outcome_mutation$OR[i],
                                              Table_outcome_mutation$OR[1] - Table_outcome_mutation$OR[i],
                                              Table_outcome_mutation$N[1] - Table_outcome_mutation$N[i] -Table_outcome_mutation$OR[1] + Table_outcome_mutation$OR[i]), nrow = 2, byrow = TRUE))
          Table_outcome_mutation$OR_OddsRatio[i] = format_p(min(unname(test_result$estimate), 99999), digits = 3)
          Table_outcome_mutation$OR_Odds_pValue[i] = format_p(unname(test_result$p.value), digits = 3)
          test_result <- fisher.test(matrix(c(Table_outcome_mutation$DC[i],
                                              Table_outcome_mutation$N[i] - Table_outcome_mutation$DC[i],
                                              Table_outcome_mutation$DC[1] - Table_outcome_mutation$DC[i],
                                              Table_outcome_mutation$N[1] - Table_outcome_mutation$N[i] -Table_outcome_mutation$DC[1] + Table_outcome_mutation$DC[i]), nrow = 2, byrow = TRUE))
          Table_outcome_mutation$DC_OddsRatio[i] = format_p(min(unname(test_result$estimate), 99999), digits = 3)
          Table_outcome_mutation$DC_Odds_pValue[i] = format_p(unname(test_result$p.value), digits = 3)
        }
        Table_outcome_mutation_total = rbind(Table_outcome_mutation_total, Table_outcome_mutation)
      }
    }
  }
  create_datatable_with_confirm(Table_outcome_mutation_total, "FELIS downloaded raw data in Factors for response rate tab")
})


output$figure_drug_8 = render_gt({
  req(OUTPUT_DATA$drug_analysis_Data_forest_tmp_8,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_8_univariant
  )
  Data_forest_tmp_8 = OUTPUT_DATA$drug_analysis_Data_forest_tmp_8
  Data_forest_tmp_8_univariant = OUTPUT_DATA$drug_analysis_Data_forest_tmp_8_univariant

  if(length(colnames(Data_forest_tmp_8)) > 1){
    Data_forest_tmp_8 = data.frame( lapply(Data_forest_tmp_8, as.factor) )
    Data_forest_tmp_8_univariant = data.frame( lapply(Data_forest_tmp_8_univariant, as.factor) )
    if("Histology" %in% colnames(Data_forest_tmp_8)){
      Data_forest_tmp_8$Histology <- relevel(Data_forest_tmp_8$Histology, ref=names(sort(table(Data_forest_tmp_8$Histology),decreasing = T))[[1]])
    }
    if("Histology" %in% colnames(Data_forest_tmp_8_univariant)){
      Data_forest_tmp_8_univariant$Histology <- relevel(Data_forest_tmp_8_univariant$Histology, ref=names(sort(table(Data_forest_tmp_8_univariant$Histology),decreasing = T))[[1]])
    }
    colnames_tmp = colnames(Data_forest_tmp_8)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_8) = colnames_tmp
    colnames_tmp = colnames(Data_forest_tmp_8_univariant)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_8_univariant) = colnames_tmp
    univ_tab <- Data_forest_tmp_8_univariant %>%
      tbl_uvregression(                         ## 単変量解析の表を生成
        method = glm,                           ## 実行したい回帰（一般化線形モデル）を定義
        y = DCR_data,
        method.args = list(family = binomial),  ## 実行したい glm のタイプを定義（ここではロジスティック）
        exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
      ) |>
      add_global_p() |> # add global p-value
      add_n(location = "level") |>
      add_nevent(location = "level") |> # add number of events of the outcome
      add_q() |> # adjusts global p-values for multiple testing
      bold_p() |> # bold p-values under a given threshold (default 0.05)
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
      bold_labels()
    m1 <- glm(DCR_data ~., data = Data_forest_tmp_8, family = binomial(link = "logit"))
    final_mv_reg <- m1 %>%
      stats::step(direction = "backward", trace = FALSE)
    if(!is.null(final_mv_reg$xlevels)){
      mv_tab = final_mv_reg |>
        tbl_regression(                         ## 多変量解析の表を生成
          exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
        ) |>
        add_global_p() |> # add global p-value
        add_q() |> # adjusts global p-values for multiple testing
        bold_p() |> # bold p-values under a given threshold (default 0.05)
        apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
        bold_labels()
      tbl_merge(
        tbls = list(univ_tab, mv_tab),
        tab_spanner = c("**Univariate**", "**Multivariable**")) |>
        modify_caption(paste("Factors for disease control rate, ", paste(input$drug, collapse = ";"), "multivariable regression (Only patients with >2 events observed)")) |> as_gt()
    } else {
      univ_tab |>
        modify_caption(paste("Factors for disease control rate,", paste(input$drug, collapse = ";"), "(Only patients with >2 events observed), no significant factor in multivariable analysis")) |> as_gt()
    }
  }
})
# Helper function to create column name mapping
figure_drug_6_create_column_mapping <- function() {
  mapping <- list(
    "Lymph_met" = "Lymphatic metastasis",
    "Liver_met" = "Liver metastasis",
    "Lung_met" = "Lung metastasis",
    "Bone_met" = "Bone metastasis",
    "Brain_met" = "Brain metastasis",
    "EP_option" = "Treatment recommended",
    "PS" = "Performance status",
    "Lines" = "CTx line"
  )
  return(mapping)
}

# Helper function to rename columns and factors
figure_drug_6_rename_columns_and_factors <- function(data, factor_names_main, factor_names_univ) {
  mapping <- figure_drug_6_create_column_mapping()

  # Rename data columns
  colnames_tmp <- colnames(data)
  for (old_name in names(mapping)) {
    colnames_tmp[colnames_tmp == old_name] <- mapping[[old_name]]
  }
  colnames(data) <- colnames_tmp

  # Rename factor names
  rename_factor_vector <- function(factor_vec) {
    for (old_name in names(mapping)) {
      factor_vec[factor_vec == old_name] <- mapping[[old_name]]
    }
    return(factor_vec)
  }

  factor_names_main_renamed <- rename_factor_vector(factor_names_main)
  factor_names_univ_renamed <- rename_factor_vector(factor_names_univ)

  return(list(
    data = data,
    factor_names_main = factor_names_main_renamed,
    factor_names_univ = factor_names_univ_renamed
  ))
}

# Helper function to filter factor names
figure_drug_6_filter_factor_names <- function(factor_names_main, factor_names_univ,
                                              remove_cluster = TRUE, remove_genes = FALSE, gene_names = NULL) {

  # Remove Cluster if requested
  if (remove_cluster) {
    factor_names_main <- factor_names_main[factor_names_main != "Cluster"]
    factor_names_univ <- factor_names_univ[factor_names_univ != "Cluster"]
  }

  # Remove genes if requested
  if (remove_genes && !is.null(gene_names)) {
    factor_names_main <- factor_names_main[!factor_names_main %in% gene_names]
    factor_names_univ <- factor_names_univ[!factor_names_univ %in% gene_names]
  }

  return(list(
    main = factor_names_main,
    univ = factor_names_univ
  ))
}

# Helper function to perform Cox regression analysis
figure_drug_6_perform_cox_regression_analysis <- function(data, factor_names_main, factor_names_univ, drug_names) {

  # Create Cox regression formula
  formula_str <- paste0("Surv(Drug_length, Drug_length_censor) ~ ",
                        paste(paste0("`", factor_names_main, "`"), collapse = "+"))
  cox_formula <- as.formula(formula_str)

  # Fit Cox regression model
  linelistsurv_cox <- coxph(
    formula = cox_formula,
    data = data,
    control = coxph.control(iter.max = 50)
  )

  # Univariate analysis
  univ_tab <- data %>%
    tbl_uvregression(
      method = coxph,
      y = Surv(time = Drug_length, event = Drug_length_censor),
      include = factor_names_univ,
      exponentiate = TRUE
    ) %>%
    add_global_p() %>%
    add_n(location = "level") %>%
    add_nevent(location = "level") %>%
    add_q() %>%
    bold_p() %>%
    apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
    bold_labels()

  # Stepwise selection for multivariable model
  final_mv_reg <- linelistsurv_cox %>%
    stats::step(direction = "backward", trace = FALSE)

  # Generate final table
  caption_base <- paste0("Hazard ratio for Treatment discontinuation with ",
                         paste(drug_names, collapse = ";"))

  if (!is.null(final_mv_reg$xlevels)) {
    mv_tab <- final_mv_reg %>%
      tbl_regression(exponentiate = TRUE) %>%
      add_global_p() %>%
      add_q() %>%
      bold_p() %>%
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
      bold_labels()

    result <- tbl_merge(
      tbls = list(univ_tab, mv_tab),
      tab_spanner = c("**Univariate**", "**Multivariable**")
    ) %>%
      modify_caption(paste0(caption_base, " (Only factors with >2 events observed)")) %>%
      as_gt()
  } else {
    result <- univ_tab %>%
      modify_caption(paste0(caption_base, ", no significant factor in multivariable analysis (Only factors with >2 events observed)")) %>%
      as_gt()
  }

  return(result)
}

# Main function for treatment discontinuation analysis
figure_drug_6_generate_discontinuation_table <- function(data, factor_names_main, factor_names_univ, drug_names,
                                                         remove_cluster = TRUE, remove_genes = FALSE, gene_names = NULL) {

  # Filter factor names based on requirements
  filtered_factors <- figure_drug_6_filter_factor_names(
    factor_names_main = factor_names_main,
    factor_names_univ = factor_names_univ,
    remove_cluster = remove_cluster,
    remove_genes = remove_genes,
    gene_names = gene_names
  )

  # Check if we have factors to analyze
  if (length(filtered_factors$main) == 0) {
    return(NULL)
  }

  # Rename columns and factors
  renamed_data <- figure_drug_6_rename_columns_and_factors(
    data = data,
    factor_names_main = filtered_factors$main,
    factor_names_univ = filtered_factors$univ
  )

  # Perform Cox regression analysis
  result <- figure_drug_6_perform_cox_regression_analysis(
    data = renamed_data$data,
    factor_names_main = renamed_data$factor_names_main,
    factor_names_univ = renamed_data$factor_names_univ,
    drug_names = drug_names
  )

  return(result)
}

# Refactored output functions
output$figure_drug_6_1 <- render_gt({
  req(OUTPUT_DATA$drug_analysis_gene_names,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_x,
      OUTPUT_DATA$drug_analysis_Factor_names_x,
      OUTPUT_DATA$drug_analysis_Factor_names_x_univariant
  )
  figure_drug_6_generate_discontinuation_table(
    data = OUTPUT_DATA$drug_analysis_Data_forest_tmp_x,
    factor_names_main = OUTPUT_DATA$drug_analysis_Factor_names_x,
    factor_names_univ = OUTPUT_DATA$drug_analysis_Factor_names_x_univariant,
    drug_names = input$drug,
    remove_cluster = TRUE,
    remove_genes = FALSE
  )
})

output$figure_drug_6_2 <- render_gt({
  req(OUTPUT_DATA$drug_analysis_gene_names,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_x,
      OUTPUT_DATA$drug_analysis_Factor_names_x,
      OUTPUT_DATA$drug_analysis_Factor_names_x_univariant
  )
  figure_drug_6_generate_discontinuation_table(
    data = OUTPUT_DATA$drug_analysis_Data_forest_tmp_x,
    factor_names_main = OUTPUT_DATA$drug_analysis_Factor_names_x,
    factor_names_univ = OUTPUT_DATA$drug_analysis_Factor_names_x_univariant,
    drug_names = input$drug,
    remove_cluster = FALSE,
    remove_genes = TRUE,
    gene_names = OUTPUT_DATA$drug_analysis_gene_names
  )
})
# Helper function for data preprocessing
figure_drug_preprocess_analysis_data <- function(data_main, data_univ, remove_genes = FALSE, gene_names = NULL) {
  # Remove Cluster column if exists
  if("Cluster" %in% colnames(data_main) && !remove_genes) {
    data_main <- data_main %>% dplyr::select(-Cluster)
  }
  if("Cluster" %in% colnames(data_univ) && !remove_genes) {
    data_univ <- data_univ %>% dplyr::select(-Cluster)
  }

  # Remove gene columns if requested
  if(remove_genes && !is.null(gene_names)) {
    existing_genes <- intersect(colnames(data_main), gene_names)
    data_main[, (existing_genes) := NULL]

    existing_genes <- intersect(colnames(data_univ), gene_names)
    data_univ[, (existing_genes) := NULL]
  }

  return(list(main = data_main, univ = data_univ))
}

# Helper function for data type conversion and column renaming
figure_drug_convert_and_rename_data <- function(data_main, data_univ) {
  # Convert to factors
  data_main <- data.frame(lapply(data_main, as.factor))
  data_univ <- data.frame(lapply(data_univ, as.factor))

  # Convert Drug_length to months
  if("Drug_length" %in% colnames(data_main)) {
    data_main$Drug_length <- as.integer(data_main$Drug_length) / 365.25 * 12
  }
  if("Drug_length" %in% colnames(data_univ)) {
    data_univ$Drug_length <- as.integer(data_univ$Drug_length) / 365.25 * 12
  }

  # Relevel Histology factor
  if("Histology" %in% colnames(data_main)) {
    data_main$Histology <- relevel(data_main$Histology,
                                   ref = names(sort(table(data_main$Histology), decreasing = TRUE))[[1]])
  }
  if("Histology" %in% colnames(data_univ)) {
    data_univ$Histology <- relevel(data_univ$Histology,
                                   ref = names(sort(table(data_univ$Histology), decreasing = TRUE))[[1]])
  }

  # Relevel Cluster factor (for figure_drug_9_2)
  if("Cluster" %in% colnames(data_main)) {
    data_main$Cluster <- relevel(data_main$Cluster,
                                 ref = names(sort(table(data_main$Cluster), decreasing = TRUE))[[1]])
  }
  if("Cluster" %in% colnames(data_univ)) {
    data_univ$Cluster <- relevel(data_univ$Cluster,
                                 ref = names(sort(table(data_univ$Cluster), decreasing = TRUE))[[1]])
  }

  # Rename columns
  rename_columns <- function(data) {
    colnames_tmp <- colnames(data)
    colnames_tmp[colnames_tmp == "Drug_length"] <- "Treatment time (month)"
    colnames_tmp[colnames_tmp == "Lymph_met"] <- "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] <- "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] <- "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] <- "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] <- "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] <- "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] <- "CTx line"
    colnames_tmp[colnames_tmp == "ORR_data"] <- "Objective response achieved"
    colnames(data) <- colnames_tmp
    return(data)
  }

  data_main <- rename_columns(data_main)
  data_univ <- rename_columns(data_univ)

  return(list(main = data_main, univ = data_univ))
}

# Helper function for performing regression analysis
figure_drug_perform_adverse_effect_analysis <- function(data_main, data_univ, drug_names) {
  # Univariate analysis
  univ_tab <- data_univ %>%
    tbl_uvregression(
      method = glm,
      y = Adverse_effect,
      method.args = list(family = binomial),
      exponentiate = TRUE
    ) %>%
    add_global_p() %>%
    add_n(location = "level") %>%
    add_nevent(location = "level") %>%
    add_q() %>%
    bold_p() %>%
    apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
    bold_labels()

  # Multivariable analysis
  m1 <- glm(Adverse_effect ~ ., data = data_main, family = binomial(link = "logit"))
  final_mv_reg <- m1 %>% stats::step(direction = "backward", trace = FALSE)

  # Generate final table
  if(!is.null(final_mv_reg$xlevels)) {
    mv_tab <- final_mv_reg %>%
      tbl_regression(exponentiate = TRUE) %>%
      add_global_p() %>%
      add_q() %>%
      bold_p() %>%
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
      bold_labels()

    result <- tbl_merge(
      tbls = list(univ_tab, mv_tab),
      tab_spanner = c("**Univariate**", "**Multivariable**")
    ) %>%
      modify_caption(paste("Factors for adverse effect,",
                           paste(drug_names, collapse = ";"),
                           "multivariable regression (Only patients with >2 events observed)")) %>%
      as_gt()
  } else {
    result <- univ_tab %>%
      modify_caption(paste("Factors for adverse effect,",
                           paste(drug_names, collapse = ";"),
                           "(Only patients with >2 events observed), no significant factor in multivariable analysis")) %>%
      as_gt()
  }

  return(result)
}

# Main function for adverse effect analysis
figure_drug_generate_adverse_effect_table <- function(data_main, data_univ, drug_names, remove_genes = FALSE, gene_names = NULL) {
  if(length(colnames(data_main)) <= 1) {
    return(NULL)
  }

  # Preprocess data
  processed_data <- figure_drug_preprocess_analysis_data(data_main, data_univ, remove_genes, gene_names)

  # Convert data types and rename columns
  converted_data <- figure_drug_convert_and_rename_data(processed_data$main, processed_data$univ)

  # Perform analysis
  result <- figure_drug_perform_adverse_effect_analysis(converted_data$main, converted_data$univ, drug_names)

  return(result)
}

# Refactored output functions
output$figure_drug_9_1 <- render_gt({
  req(OUTPUT_DATA$drug_analysis_Data_forest_tmp_9,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_9_univariant,
      OUTPUT_DATA$drug_analysis_gene_names)
  figure_drug_generate_adverse_effect_table(
    data_main = OUTPUT_DATA$drug_analysis_Data_forest_tmp_9,
    data_univ = OUTPUT_DATA$drug_analysis_Data_forest_tmp_9_univariant,
    drug_names = input$drug,
    remove_genes = FALSE
  )
})

output$figure_drug_9_2 <- render_gt({
  req(OUTPUT_DATA$drug_analysis_Data_forest_tmp_9,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_9_univariant,
      OUTPUT_DATA$drug_analysis_gene_names)
  figure_drug_generate_adverse_effect_table(
    data_main = OUTPUT_DATA$drug_analysis_Data_forest_tmp_9,
    data_univ = OUTPUT_DATA$drug_analysis_Data_forest_tmp_9_univariant,
    drug_names = input$drug,
    remove_genes = TRUE,
    gene_names = unique(c(OUTPUT_DATA$drug_analysis_gene_names,
                          OUTPUT_DATA$drug_analysis_significant_genes_AE))
  )
})

output$select_table_var_volcano_2_AE = renderUI({
  req(OUTPUT_DATA$drug_analysis_regimen_choice_RECIST_AE)
  radioButtons(
    inputId = "table_var_volcano_2_AE",
    label = "Regimen name",
    choices = OUTPUT_DATA$drug_analysis_regimen_choice_RECIST_AE,
    selected = OUTPUT_DATA$drug_analysis_regimen_choice_RECIST_AE[1])
})

output$figure_survival_drug = renderPlot({
  req(figure_drug_trigger())
  figure_survival_drug_var = isolate(req(input$figure_survival_drug_var))
  drug_survival_gene = isolate(req(input$drug_survival_gene))
  Data_drug_survival_1_save = isolate(req(OUTPUT_DATA$Data_drug_Data_drug_survival_1_save))
  Data_MAF_target = isolate(req(OUTPUT_DATA$Data_drug_Data_MAF_target))
  Data_drug_survival_1_save$time_palliative_final = Data_drug_survival_1_save$TTD
  Data_drug_survival_1_save$time_palliative_enroll = Data_drug_survival_1_save$TTE
  Data_drug_survival_1_save$time_enroll_final = Data_drug_survival_1_save$TTD - Data_drug_survival_1_save$TTE
  withProgress(message = "Drawing...", {
    if(figure_survival_drug_var == "AE_OS"){
      Data_drug_survival_1 = Data_drug_survival_1_save %>%
        dplyr::mutate(mutation = case_when(
          Adverse_effect == "AE (+)" ~ "Adverse effect (+)",
          TRUE ~ "Adverse effect (-)"
        )) %>%
        dplyr::arrange(ID) %>%
        dplyr::filter(Regimen_whole == "Selected regimens")
      survival_drug_function(Data_drug_survival_1 = Data_drug_survival_1,
                             name_1 = "Adverse effect (+)",
                             name_2 = "Adverse effect (-)",
                             name_1_short = "AE(+)",
                             name_2_short = "AE(-)")
    } else if(figure_survival_drug_var == "regimen_OS"){
      Data_drug_survival_1 = Data_drug_survival_1_save %>% dplyr::mutate(
        mutation = case_when(
          Regimen_whole == "Selected regimens" ~ paste0(
              paste(head(input$drug, 2), collapse = ";"),
              if (length(input$drug) > 2) "..." else ""
            ),
          TRUE ~ "Without above, only other drugs"
        ))
      survival_drug_function(Data_drug_survival_1 = Data_drug_survival_1,
                             name_1 = paste0(
                               paste(head(input$drug, 2), collapse = ";"),
                               if (length(input$drug) > 2) "..." else ""
                             ),
                             name_2 = "Without above, only other drugs",
                             name_1_short = paste0(
                               paste(head(input$drug, 1), collapse = ";"),
                               if (length(input$drug) > 1) "..." else ""
                             ),
                             name_2_short = "Others"
      )
    } else if(!is.null(drug_survival_gene) && figure_survival_drug_var == "Mutation_OS") {
      Data_ID_target_mutation_original = unique((Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% drug_survival_gene))$Tumor_Sample_Barcode)
      Data_drug_survival_1 = Data_drug_survival_1_save %>%
        dplyr::mutate(mutation = case_when(
          ID %in% Data_ID_target_mutation_original ~ paste0(paste0(
            paste(head(drug_survival_gene, 2), collapse = ";"),
            if (length(drug_survival_gene) > 2) "..." else ""
          ), " mut(+)"),
          TRUE ~ paste0(paste0(
            paste(head(drug_survival_gene, 2), collapse = ";"),
            if (length(drug_survival_gene) > 2) "..." else ""
          ), " mut(-)")
        )) %>%
        dplyr::arrange(ID) %>%
        dplyr::filter(Regimen_whole == "Selected regimens")
      survival_drug_function(Data_drug_survival_1 = Data_drug_survival_1,
                             name_1 = paste0(paste0(
                               paste(head(drug_survival_gene, 2), collapse = ";"),
                               if (length(drug_survival_gene) > 2) "..." else ""
                             ), " mut(+)"),
                             name_2 = paste0(paste0(
                               paste(head(drug_survival_gene, 2), collapse = ";"),
                               if (length(drug_survival_gene) > 2) "..." else ""
                             ), " mut(-)"),
                             name_1_short = "mut(+)",
                             name_2_short = "mut(-)"
      )
    } else if(figure_survival_drug_var == "Drug_OS" &&
              !is.null(input$drug_group_1) &&
              !is.null(input$drug_group_2)) {
      Data_drug_survival_1 = Data_drug_survival_1_save %>% dplyr::mutate(
        mutation = case_when(
          Drug %in% input$drug_group_1 ~ input$drug_group_1_name,
          Drug %in% input$drug_group_2 ~ input$drug_group_2_name,
          TRUE ~ "Other"
        )) %>%
        dplyr::filter(mutation != "Other")
      survival_drug_function(Data_drug_survival_1 = Data_drug_survival_1,
                             name_1 = input$drug_group_1_name,
                             name_2 = input$drug_group_2_name,
                             name_1_short = "drug group 1",
                             name_2_short = "drug group 2"
      )
    } else {
      return(NULL)
    }
  })
})

survival_drug_function = function(Data_drug_survival_1, name_1, name_2, name_1_short, name_2_short){
  figure_survival_drug_var_2 = isolate(req(input$figure_survival_drug_var_2))
  if(is.null(figure_survival_drug_var_2)){
    return(NULL)
  } else if(figure_survival_drug_var_2 != "Bayes"){
    adjustment = figure_survival_drug_var_2 == "Number at risk"
    traditional_fit = survfit(Surv(event = censor,
                                   time = time_palliative_final) ~ mutation,
                              data = Data_drug_survival_1)
    if(length(unique(Data_drug_survival_1$mutation)) > 1){
      survival_compare_and_plot_CTx(
        data = Data_drug_survival_1,
        time_var1 = "time_palliative_enroll",
        time_var2 = "time_palliative_final",
        status_var = "censor",
        group_var = "mutation",
        plot_title = paste0("Survival from drug initiation, ", paste(traditional_fit[[1]],collapse = "/"), " patients"),
        adjustment = adjustment,
        color_var_surv_CTx_1 = input$color_var_surv_CTx_1,
        group_labels = c(paste0(name_1, ", Adjusted for bias"),
                         paste0(name_2, ", Adjusted for bias"))
      )
    }
  } else {
    if(sum(Data_drug_survival_1$censor) > 0){
      max_samples = isolate(input$figure_Drug_Bayes_max_samples)
      max_samples = ifelse(is.null(max_samples), 2000, max_samples)
      if (nrow(Data_drug_survival_1) > max_samples){
        sample_indices = sample(nrow(Data_drug_survival_1),max_samples)
        Data_drug_survival_1 = Data_drug_survival_1[sample_indices, ]
        ENTIRE = "subsampled"
      } else {
        Data_drug_survival_1 = Data_drug_survival_1
        ENTIRE = "whole"
      }

      traditional_fit = survfit(Surv(event = censor,
                                     time = time_palliative_final) ~ mutation,
                                data = Data_drug_survival_1)
      if(nrow(Data_drug_survival_1) != traditional_fit[[1]][[1]]){
        Data_drug_survival_1$gene_mut = Data_drug_survival_1$mutation == name_2
        Data_drug_survival_1_curve = Data_drug_survival_1
        Data_survival_tmp_pos = Data_drug_survival_1_curve %>%
          dplyr::filter(gene_mut == TRUE)
        if(sum(Data_survival_tmp_pos$censor) > 0){
          fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_pos)
          time_10 = min(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_enroll_final)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 0
          survival_rate_90 = 1
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_enroll_final)*0.9)
        Data_survival_tmp = Data_survival_tmp_pos %>% dplyr::filter(
          time_enroll_final >= time_10 &
            time_enroll_final <= time_90)
        Data_survival_tmp$time_enroll_final =
          Data_survival_tmp$time_enroll_final - (time_10 - 1)
        Data_survival_tmp2 = Data_survival_tmp_pos %>% dplyr::filter(
          time_enroll_final > time_90) %>%
          dplyr::arrange(time_enroll_final)

        Data_survival_tmp2$time_enroll_final = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp2)
        Data_survival_tmp = rbind(Data_survival_tmp, Data_survival_tmp2[idx[idx < max(idx)*(1-survival_rate_10)],])

        Data_survival_tmp_neg = Data_drug_survival_1_curve %>%
          dplyr::filter(gene_mut == FALSE)
        if(sum(Data_survival_tmp_neg$censor) > 0){
          fit = survfit(Surv(time_enroll_final, censor) ~ 1, data = Data_survival_tmp_neg)
          time_10 = min(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_enroll_final)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 0
          survival_rate_90 = 1
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_enroll_final)*0.9)
        Data_survival_tmp3 = Data_survival_tmp_neg %>% dplyr::filter(
          time_enroll_final >= time_10 &
            time_enroll_final <= time_90)
        Data_survival_tmp3$time_enroll_final =
          Data_survival_tmp3$time_enroll_final - (time_10 - 1)
        Data_survival_tmp4 = Data_survival_tmp_neg %>% dplyr::filter(
          time_enroll_final > time_90) %>%
          dplyr::arrange(time_enroll_final)

        Data_survival_tmp4$time_enroll_final = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp4)
        Data_survival_tmp3 = rbind(Data_survival_tmp3, Data_survival_tmp4[idx[idx < max(idx)*(1-survival_rate_10)],])
        Data_drug_survival_1 = rbind(Data_survival_tmp, Data_survival_tmp3)

        Data_survival_tmp_pos = Data_drug_survival_1_curve %>%
          dplyr::filter(gene_mut == TRUE)
        if(sum(Data_survival_tmp_pos$censor) > 0){
          fit = survfit(Surv(time_palliative_enroll, censor) ~ 1, data = Data_survival_tmp_pos)
          time_10 = min(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 0
          survival_rate_90 = 1
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_pos$time_palliative_enroll)*0.9)
        Data_survival_tmp5 = Data_survival_tmp_pos %>% dplyr::filter(
          time_palliative_enroll >= time_10 &
            time_palliative_enroll <= time_90)
        Data_survival_tmp5$time_palliative_enroll =
          Data_survival_tmp5$time_palliative_enroll - (time_10 - 1)
        Data_survival_tmp6 = Data_survival_tmp_pos %>% dplyr::filter(
          time_palliative_enroll > time_90) %>%
          dplyr::arrange(time_palliative_enroll)
        Data_survival_tmp6$time_palliative_enroll = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp6)
        Data_survival_tmp5 = rbind(Data_survival_tmp5, Data_survival_tmp6[idx[idx < max(idx)*(1-survival_rate_10)],])

        Data_survival_tmp_neg = Data_drug_survival_1_curve %>%
          dplyr::filter(gene_mut == FALSE)

        if(sum(Data_survival_tmp_neg$censor) > 0){
          fit = survfit(Surv(time_palliative_enroll, censor) ~ 1, data = Data_survival_tmp_neg)
          time_10 = min(90, quantile(fit, 0.10)$quantile[[1]],na.rm = T)
          survival_rate_10 = summary(fit,times = time_10)[[6]]
          survival_rate_90 = summary(fit,times = max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)[[6]]
        } else {
          time_10 = 90
          survival_rate_10 = 0
          survival_rate_90 = 1
        }
        time_90 = max(time_10 + 1, max(Data_survival_tmp_neg$time_palliative_enroll)*0.9)
        Data_survival_tmp7 = Data_survival_tmp_neg %>% dplyr::filter(
          time_palliative_enroll >= time_10 &
            time_palliative_enroll <= time_90)
        Data_survival_tmp7$time_palliative_enroll =
          Data_survival_tmp7$time_palliative_enroll - (time_10 - 1)
        Data_survival_tmp8 = Data_survival_tmp_neg %>% dplyr::filter(
          time_palliative_enroll > time_90) %>%
          dplyr::arrange(time_palliative_enroll)
        Data_survival_tmp8$time_palliative_enroll = time_90 - (time_10 - 1)
        idx <- 1:nrow(Data_survival_tmp8)
        Data_survival_tmp7 = rbind(Data_survival_tmp7, Data_survival_tmp8[idx[idx < max(idx)*(1-survival_rate_10)],])

        Median_entry_pos = median(Data_survival_tmp5$time_palliative_enroll)
        Median_entry_neg = median(Data_survival_tmp7$time_palliative_enroll)
        Data_drug_survival_1_2 = rbind(Data_survival_tmp5, Data_survival_tmp7)

        Data_drug_survival_1 = Data_drug_survival_1 %>%
          dplyr::mutate(delay_1 = case_when(
            time_palliative_enroll <= Median_entry_pos &
              gene_mut == TRUE ~ "SPT to entry early",
            time_palliative_enroll > Median_entry_pos &
              gene_mut == TRUE ~ "SPT to entry later",
            time_palliative_enroll <= Median_entry_neg &
              gene_mut == FALSE ~ "SPT to entry early",
            time_palliative_enroll > Median_entry_neg &
              gene_mut == FALSE ~ "SPT to entry later"))

        if(sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 == "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 == "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == TRUE & delay_1 != "SPT to entry later"))$censor) >= 1 &
           sum((Data_drug_survival_1 %>% dplyr::filter(gene_mut == FALSE & delay_1 != "SPT to entry later"))$censor) >= 1){

          stan_weibull_survival_model_data <-
            list(
              Median_pos = as.integer(Median_entry_pos),
              Median_neg = as.integer(Median_entry_neg),
              Nobs_early_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Nobs_late_pos = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Ncen_early_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Ncen_late_pos = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE),
              Nobs_early_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Nobs_late_neg = sum(Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Ncen_early_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Ncen_late_neg = sum(Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE),
              Nexp_pos = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE]),
              Nexp_neg = length(Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE]),
              ybef_exp_pos = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut != FALSE],
              ybef_exp_neg = Data_drug_survival_1_2$time_palliative_enroll[Data_drug_survival_1_2$gene_mut == FALSE],
              yobs_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              yobs_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              ycen_early_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              ycen_late_pos = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut != FALSE],
              yobs_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              yobs_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 1 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              ycen_early_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 != "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE],
              ycen_late_neg = Data_drug_survival_1$time_enroll_final[Data_drug_survival_1$censor == 0 & Data_drug_survival_1$delay_1 == "SPT to entry later" & Data_drug_survival_1$gene_mut == FALSE]
            )

          if(DOCKER){
            stan_model_compiled <- readRDS(stan_model_factor_save_path_cmdstan)
            stan_weibull_survival_model_fit <- stan_model_compiled$sample(
              data = stan_weibull_survival_model_data,
              seed = 1234,
              chains = 4,
              parallel_chains = PARALLEL,
              iter_sampling = ITER,  # (iter-warmup)
              iter_warmup = ITER / 2,adapt_delta=0.99, max_treedepth = 10,
              thin = 1,
              refresh = 500
            )
          } else {
            mod <- readRDS(stan_model_factor_save_path_rstan)
            stan_weibull_survival_model_fit <-
              rstan::sampling(object = mod,
                          data = stan_weibull_survival_model_data,
                          control=list(adapt_delta=0.99, max_treedepth = 10),
                          seed=1234, chains=4, iter=3 / 2 * ITER, warmup=ITER / 2, thin=1, refresh=500, verbose = FALSE)
          }

          stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_weibull_survival_model_fit)
          stan_weibull_survival_model_draws_total <-
            stan_weibull_survival_model_draws %>%
            dplyr::select(.chain, .iteration, .draw, starts_with("yhat_total")) %>%
            tidyr::gather(key = key, value = yhat_total, starts_with("yhat_total")) %>%
            dplyr::select(-key)
          stan_weibull_survival_model_draws_total_exp <-
            stan_weibull_survival_model_draws %>%
            dplyr::select(.chain, .iteration, .draw, starts_with("yhat_factor_total")) %>%
            tidyr::gather(key = key, value = yhat_total_exp, starts_with("yhat_factor_total")) %>%
            dplyr::select(-key)

          median_os_list_weibull_total_exp = matrix(stan_weibull_survival_model_draws_total_exp$yhat_total_exp, nrow=4 * ITER)
          median_os_list_weibull_total = matrix(stan_weibull_survival_model_draws_total$yhat_total, nrow=4 * ITER)
          median_os_list_weibull_total_exp = rowMedians(median_os_list_weibull_total_exp,na.rm = T)
          median_os_list_weibull_total = rowMedians(median_os_list_weibull_total,na.rm = T)
          median_os_list_weibull_bias = median_os_list_weibull_total - median_os_list_weibull_total_exp

          median_os_summary_weibull_total = quantile(median_os_list_weibull_total, probs = c(0.025, 0.5, 0.975))
          median_os_summary_weibull_total_exp = quantile(median_os_list_weibull_total_exp, probs = c(0.025, 0.5, 0.975))
          median_os_summary_weibull_bias = quantile(median_os_list_weibull_bias, probs = c(0.025, 0.5, 0.975))
          yhat_weibull_sort_total = sort(stan_weibull_survival_model_draws_total$yhat_total)
          yhat_weibull_sort_total_exp = sort(stan_weibull_survival_model_draws_total_exp$yhat_total_exp)
          yhat_weibull_sort_average_total = yhat_weibull_sort_total[1:nrow(Data_drug_survival_1_curve) * as.integer(length(yhat_weibull_sort_total)/nrow(Data_drug_survival_1_curve))]
          yhat_weibull_sort_average_total_exp = yhat_weibull_sort_total_exp[1:nrow(Data_drug_survival_1_curve) * as.integer(length(yhat_weibull_sort_total_exp)/nrow(Data_drug_survival_1_curve))]

          Data_drug_survival_1_curve$factor_neg = yhat_weibull_sort_average_total
          Data_drug_survival_1_curve$factor_pos = yhat_weibull_sort_average_total_exp
          Data_drug_survival_1_curve$factor_neg[Data_drug_survival_1_curve$factor_neg > 10000] = 10000
          Data_drug_survival_1_curve$factor_pos[Data_drug_survival_1_curve$factor_pos > 10000] = 10000

          traditional_fit = survfit(Surv(event = censor,
                                         time = time_palliative_final) ~ gene_mut,
                                    data = Data_drug_survival_1_curve)
          simulation_fit_neg = survfit(Surv(event = rep(1,nrow(Data_drug_survival_1_curve)),
                                            time = factor_neg) ~ 1,
                                       data = Data_drug_survival_1_curve)
          simulation_fit_pos = survfit(Surv(event = rep(1,nrow(Data_drug_survival_1_curve)),
                                            time = factor_pos) ~ 1,
                                       data = Data_drug_survival_1_curve)

          g_surv_drug_tmp = ggsurvplot(
            fit = list(traditional_fit = traditional_fit,
                       simulation_fit_neg = simulation_fit_neg,
                       simulation_fit_pos = simulation_fit_pos),
            combine = TRUE,
            data = Data_drug_survival_1_curve,
            xlab = "Time from CTx induction to final observation (months)",
            ylab = paste0("Survival Probability, ", ENTIRE, " cohort"),
            censor = TRUE,
            font.title = 8,
            font.subtitle = 8,
            font.main = 8,
            font.submain = 8,
            font.caption = 8,
            font.legend = 8,
            surv.median.line = "hv",
            conf.int = FALSE,
            pval = FALSE,
            risk.table = TRUE,
            risk.table.y.text = FALSE,
            tables.theme = clean_theme(),
            legend = c(0.8,0.8),
            xlim = c(0, max(Data_drug_survival_1_curve$time_palliative_final) * 1.05),
            cumevents = FALSE,
            cumcensor = FALSE,
            break.x.by = 365.25 * 2.5,
            xscale = "d_m",
            legend.labs = c(paste0(name_1, ", Unadjusted"),
                            paste0(name_2, ", Unadjusted"),
                            paste0(name_1, ", Adjusted"),
                            paste0(name_2, ", Adjusted"))
          ) +
            labs(title = paste(paste0("Treated regimens including ", paste(input$drug, collapse = ";")), "/n", sum(traditional_fit[[1]]), " patients ",
                               "Median OS ",
                               format_p(summary(traditional_fit)$table[[13]] / 365.25 * 12, digits = 1)," (",
                               format_p(summary(traditional_fit)$table[[15]] / 365.25 * 12, digits = 1),"-",
                               format_p(summary(traditional_fit)$table[[17]] / 365.25 * 12, digits = 1),") (",
                               name_1_short,", unadj.)/",
                               format_p(summary(traditional_fit)$table[[14]] / 365.25 * 12, digits = 1), " (",
                               format_p(summary(traditional_fit)$table[[16]] / 365.25 * 12, digits = 1),"-",
                               format_p(summary(traditional_fit)$table[[18]] / 365.25 * 12, digits = 1),") (",
                               name_2_short,", unadj.)/",
                               format_p(median_os_summary_weibull_total[[2]] / 365.25 * 12, digits = 1),
                               " (", format_p(median_os_summary_weibull_total[[1]] / 365.25 * 12, digits = 1), "-",
                               format_p(median_os_summary_weibull_total[[3]] / 365.25 * 12, digits = 1), ") (",
                               name_1_short,", adj.)/",
                               format_p(median_os_summary_weibull_total_exp[[2]] / 365.25 * 12, digits = 1), " (",
                               format_p(median_os_summary_weibull_total_exp[[1]] / 365.25 * 12, digits = 1), "-",
                               format_p(median_os_summary_weibull_total_exp[[3]] / 365.25 * 12, digits = 1),
                               ") (", name_2_short,", adj.) months",
                               sep=""),
                 subtitle = paste("Survival difference, ",name_1_short, "-", name_2_short, ": ",
                                  format_p(median_os_summary_weibull_bias[[2]] / 365.25 * 12, digits = 1), " (",
                                  format_p(median_os_summary_weibull_bias[[1]] / 365.25 * 12, digits = 1), "-",
                                  format_p(median_os_summary_weibull_bias[[3]] / 365.25 * 12, digits = 1),
                                  ") months, median (95% CI)",
                                  sep=""))
          g_surv_drug_tmp$table <- g_surv_drug_tmp$table + theme(plot.title = element_blank(),
                                                                           plot.subtitle = element_blank())
          return(g_surv_drug_tmp)
        }
      } else {return(NULL)}
    }
  }
}

output$figure_drug_7_1 = render_gt({
  req(OUTPUT_DATA$drug_analysis_Data_forest_tmp_7,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_7_univariant)
  Data_forest_tmp_7 = OUTPUT_DATA$drug_analysis_Data_forest_tmp_7
  Data_forest_tmp_7_univariant = OUTPUT_DATA$drug_analysis_Data_forest_tmp_7_univariant
  if("Cluster" %in% colnames(Data_forest_tmp_7)){
    Data_forest_tmp_7 = Data_forest_tmp_7 %>%
      dplyr::select(-Cluster)
  }
  if("Cluster" %in% colnames(Data_forest_tmp_7_univariant)){
    Data_forest_tmp_7_univariant = Data_forest_tmp_7_univariant %>%
      dplyr::select(-Cluster)
  }
  if(length(colnames(Data_forest_tmp_7)) > 1){
    Data_forest_tmp_7 = data.frame( lapply(Data_forest_tmp_7, as.factor) )
    Data_forest_tmp_7_univariant = data.frame( lapply(Data_forest_tmp_7_univariant, as.factor) )

    if("Histology" %in% colnames(Data_forest_tmp_7)){
      Data_forest_tmp_7$Histology <- relevel(Data_forest_tmp_7$Histology, ref=names(sort(table(Data_forest_tmp_7$Histology),decreasing = T))[[1]])
    }
    if("Histology" %in% colnames(Data_forest_tmp_7_univariant)){
      Data_forest_tmp_7_univariant$Histology <- relevel(Data_forest_tmp_7_univariant$Histology, ref=names(sort(table(Data_forest_tmp_7_univariant$Histology),decreasing = T))[[1]])
    }
    colnames_tmp = colnames(Data_forest_tmp_7)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_7) = colnames_tmp
    colnames_tmp = colnames(Data_forest_tmp_7_univariant)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_7_univariant) = colnames_tmp
    univ_tab <- Data_forest_tmp_7_univariant %>%
      tbl_uvregression(                         ## 単変量解析の表を生成
        method = glm,                           ## 実行したい回帰（一般化線形モデル）を定義
        y = ORR_data,
        method.args = list(family = binomial),  ## 実行したい glm のタイプを定義（ここではロジスティック）
        exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
      ) |>
      add_global_p() |> # add global p-value
      add_n(location = "level") |>
      add_nevent(location = "level") |> # add number of events of the outcome
      add_q() |> # adjusts global p-values for multiple testing
      bold_p() |> # bold p-values under a given threshold (default 0.05)
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
      bold_labels()
    m1 <- glm(ORR_data ~., data = Data_forest_tmp_7, family = binomial(link = "logit"))
    final_mv_reg <- m1 %>%
      stats::step(direction = "backward", trace = FALSE)
    if(!is.null(final_mv_reg$xlevels)){
      mv_tab = final_mv_reg |>
        tbl_regression(                         ## 多変量解析の表を生成
          exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
        ) |>
        add_global_p() |> # add global p-value
        add_q() |> # adjusts global p-values for multiple testing
        bold_p() |> # bold p-values under a given threshold (default 0.05)
        apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
        bold_labels()
      tbl_merge(
        tbls = list(univ_tab, mv_tab),
        tab_spanner = c("**Univariate**", "**Multivariable**")) |>
        modify_caption(paste("Factors for objective response rate, ", paste(input$drug, collapse = ";"), "multivariable regression (Only patients with >2 events observed)")) |> as_gt()
    } else {
      univ_tab |>
        modify_caption(paste("Factors for objective response rate,", paste(input$drug, collapse = ";"), "(Only patients with >2 events observed), no significant factor in multivariable analysis")) |> as_gt()
    }
  }
})

output$figure_drug_7_2 = render_gt({
  req(OUTPUT_DATA$drug_analysis_Data_forest_tmp_7,
      OUTPUT_DATA$drug_analysis_Data_forest_tmp_7_univariant,
      OUTPUT_DATA$drug_analysis_gene_names)
  Data_forest_tmp_7 = OUTPUT_DATA$drug_analysis_Data_forest_tmp_7
  Data_forest_tmp_7_univariant = OUTPUT_DATA$drug_analysis_Data_forest_tmp_7_univariant
  gene_names = unique(c(OUTPUT_DATA$drug_analysis_gene_names))
  existing_genes <- intersect(colnames(Data_forest_tmp_7), gene_names)
  Data_forest_tmp_7[, (existing_genes) := NULL]
  existing_genes <- intersect(colnames(Data_forest_tmp_7_univariant), gene_names)
  Data_forest_tmp_7_univariant[, (existing_genes) := NULL]
  if(length(colnames(Data_forest_tmp_7)) > 1){
    Data_forest_tmp_7 = data.frame( lapply(Data_forest_tmp_7, as.factor) )
    Data_forest_tmp_7_univariant = data.frame( lapply(Data_forest_tmp_7_univariant, as.factor) )

    if("Histology" %in% colnames(Data_forest_tmp_7)){
      Data_forest_tmp_7$Histology <- relevel(Data_forest_tmp_7$Histology, ref=names(sort(table(Data_forest_tmp_7$Histology),decreasing = T))[[1]])
    }
    if("Histology" %in% colnames(Data_forest_tmp_7_univariant)){
      Data_forest_tmp_7_univariant$Histology <- relevel(Data_forest_tmp_7_univariant$Histology, ref=names(sort(table(Data_forest_tmp_7_univariant$Histology),decreasing = T))[[1]])
    }
    colnames_tmp = colnames(Data_forest_tmp_7)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_7) = colnames_tmp
    colnames_tmp = colnames(Data_forest_tmp_7_univariant)
    colnames_tmp[colnames_tmp == "Lymph_met"] = "Lymphatic metastasis"
    colnames_tmp[colnames_tmp == "Liver_met"] = "Liver metastasis"
    colnames_tmp[colnames_tmp == "Lung_met"] = "Lung metastasis"
    colnames_tmp[colnames_tmp == "Bone_met"] = "Bone metastasis"
    colnames_tmp[colnames_tmp == "Brain_met"] = "Brain metastasis"
    colnames_tmp[colnames_tmp == "PS"] = "Performance status"
    colnames_tmp[colnames_tmp == "Lines"] = "CTx line"
    colnames(Data_forest_tmp_7_univariant) = colnames_tmp
    univ_tab <- Data_forest_tmp_7_univariant %>%
      tbl_uvregression(                         ## 単変量解析の表を生成
        method = glm,                           ## 実行したい回帰（一般化線形モデル）を定義
        y = ORR_data,
        method.args = list(family = binomial),  ## 実行したい glm のタイプを定義（ここではロジスティック）
        exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
      ) |>
      add_global_p() |> # add global p-value
      add_n(location = "level") |>
      add_nevent(location = "level") |> # add number of events of the outcome
      add_q() |> # adjusts global p-values for multiple testing
      bold_p() |> # bold p-values under a given threshold (default 0.05)
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
      bold_labels()
    m1 <- glm(ORR_data ~., data = Data_forest_tmp_7, family = binomial(link = "logit"))
    final_mv_reg <- m1 %>%
      stats::step(direction = "backward", trace = FALSE)
    if(!is.null(final_mv_reg$xlevels)){
      mv_tab = final_mv_reg |>
        tbl_regression(                         ## 多変量解析の表を生成
          exponentiate = TRUE                     ## 対数オッズ比ではなくオッズ比を得るために指数変換を指定
        ) |>
        add_global_p() |> # add global p-value
        add_q() |> # adjusts global p-values for multiple testing
        bold_p() |> # bold p-values under a given threshold (default 0.05)
        apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) %>%
        bold_labels()
      tbl_merge(
        tbls = list(univ_tab, mv_tab),
        tab_spanner = c("**Univariate**", "**Multivariable**")) |>
        modify_caption(paste("Factors for objective response rate, ", paste(input$drug, collapse = ";"), "multivariable regression (Only patients with >2 events observed)")) |> as_gt()
    } else {
      univ_tab |>
        modify_caption(paste("Factors for objective response rate,", paste(input$drug, collapse = ";"), "(Only patients with >2 events observed), no significant factor in multivariable analysis")) |> as_gt()
    }
  }
})

output$figure_volcano_1_AE = renderPlot({
  req(input$table_var_volcano_2_AE,
      OUTPUT_DATA$drug_analysis_g_volcano_AE)
  OUTPUT_DATA$drug_analysis_g_volcano_AE[[as.integer(input$table_var_volcano_2_AE)]]
})


output$select_table_var_volcano_2 = renderUI({
  req(OUTPUT_DATA$drug_analysis_regimen_choice_RECIST)
  radioButtons(
    inputId = "table_var_volcano_2",
    label = "Regimen name",
    choices = OUTPUT_DATA$drug_analysis_regimen_choice_RECIST,
    selected = OUTPUT_DATA$drug_analysis_regimen_choice_RECIST[1])
})

output$figure_volcano_1 = renderPlot({
  req(input$table_var_volcano_2,
      OUTPUT_DATA$drug_analysis_g_volcano)
  OUTPUT_DATA$drug_analysis_g_volcano[[as.integer(input$table_var_volcano_2)]]
})

output$table_volcano = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$Data_drug_df_volcano)
  create_datatable_with_confirm(OUTPUT_DATA$Data_drug_df_volcano, "FELIS downloaded raw data in Volcano plot for objective response rate tab")
})


output$select_table_var_volcano_1 = renderUI({
  req(OUTPUT_DATA$drug_analysis_regimen_choice_ToT)
  radioButtons(
    inputId = "table_var_volcano_1",
    label = "Regimen name",
    choices = OUTPUT_DATA$drug_analysis_regimen_choice_ToT,
    selected = OUTPUT_DATA$drug_analysis_regimen_choice_ToT[1])
})

output$figure_volcano_ToT_1 = renderPlot({
  req(input$table_var_volcano_1,
      OUTPUT_DATA$drug_analysis_g_volcano_ToT)
  OUTPUT_DATA$drug_analysis_g_volcano_ToT[[as.integer(input$table_var_volcano_1)]]
})

output$table_volcano_ToT = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$drug_analysis_df_volcano_ToT)
  DT::datatable(OUTPUT_DATA$drug_analysis_df_volcano_ToT,
                filter = 'top',
                extensions = c('Buttons'),
                options = list(pageLength = 100,
                               scrollX = TRUE,
                               scrollY = "1000px",
                               scrollCollapse = TRUE,
                               dom="Blfrtip",
                               buttons = c('csv', 'excel', 'copy')))
})

output$table_drug_all_summary = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$drug_analysis_Drug_summary_for_table)
  create_datatable_with_confirm(OUTPUT_DATA$drug_analysis_Drug_summary_for_table, "FELIS downloaded raw data in Drug usage data tab")
})
output$download_drug_data = downloadHandler(
  filename = "table_drug_data.csv",
  content = function(file) {
    req(OUTPUT_DATA$drug_analysis_Drug_summary_for_table)
    write_excel_csv(OUTPUT_DATA$drug_analysis_Drug_summary_for_table, file)
  }
)

output$figure_survival_drug_interactive_1 = renderPlot({
  req(OUTPUT_DATA$drug_analysis_Data_drug_TTF_comparison,
      OUTPUT_DATA$Data_drug_Data_MAF_target,
      OUTPUT_DATA$Data_drug_Data_case_target)
  Data_drug_TTF_comparison = OUTPUT_DATA$drug_analysis_Data_drug_TTF_comparison
  Data_MAF_target = OUTPUT_DATA$Data_drug_Data_MAF_target
  Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
  ID_1 = unique(Data_drug_TTF_comparison$ID)
  if(!all(is.null(input$gene_survival_drug_interactive_1_P_1))){
    ID_1 = intersect(ID_1, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_1_P_1))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_P_2))){
    ID_1 = intersect(ID_1, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_1_P_2))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_W))){
    ID_1 = setdiff(ID_1, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_1_W))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_W))){
    ID_1 = setdiff(ID_1, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_1_W))$Tumor_Sample_Barcode)
  }
  if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
    if(!all(is.null(input$gene_survival_drug_interactive_1_PD_L1))){
      ID_1 = intersect(ID_1, (Data_case_target %>% dplyr::filter(PD_L1 %in% input$gene_survival_drug_interactive_1_PD_L1))$C.CAT調査結果.基本項目.ハッシュID)
    }
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_A))){
    ID_1 = intersect(ID_1, (Data_case_target %>% dplyr::filter(YoungOld %in% input$gene_survival_drug_interactive_1_A))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_S))){
    ID_1 = intersect(ID_1, (Data_case_target %>% dplyr::filter(症例.基本情報.性別.名称. %in% input$gene_survival_drug_interactive_1_S))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_H))){
    ID_1 = intersect(ID_1, (Data_case_target %>% dplyr::filter(症例.基本情報.がん種.OncoTree. %in% input$gene_survival_drug_interactive_1_H))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_C))){
    ID_1 = intersect(ID_1, (Data_case_target %>% dplyr::filter(cluster %in% input$gene_survival_drug_interactive_1_C))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_M))){
    mets <- input$gene_survival_drug_interactive_1_M
    ID_M1 <- purrr::map(mets, function(met) {
      Data_case_target %>%
        dplyr::filter(.data[[met]] == "Yes") %>%
        dplyr::pull(C.CAT調査結果.基本項目.ハッシュID)
    }) %>% unlist() %>% unique()
    ID_1 = intersect(ID_1, ID_M1)
  }
  Data_survival_1 = Data_drug_TTF_comparison %>%
    dplyr::filter(ID %in% ID_1)
  Data_survival_1$Group = 1
  if(!all(is.null(input$gene_survival_drug_interactive_1_AE))){
    if("AE (+)" %in% input$gene_survival_drug_interactive_1_AE &
       !"AE (-)" %in% input$gene_survival_drug_interactive_1_AE){
      Data_survival_1 = Data_survival_1 %>% dplyr::filter(AE == 1)
    } else if(!"AE (+)" %in% input$gene_survival_drug_interactive_1_AE &
              "AE (-)" %in% input$gene_survival_drug_interactive_1_AE){
      Data_survival_1 = Data_survival_1 %>% dplyr::filter(AE == 0)
    }
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_EAE))){
    Data_survival_1 = Data_survival_1 %>% dplyr::filter(Early_AE %in% input$gene_survival_drug_interactive_1_EAE)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_D))){
    Data_survival_1 = Data_survival_1 %>% dplyr::filter(Drug %in% input$gene_survival_drug_interactive_1_D)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_R))){
    Data_survival_1 = Data_survival_1 %>% dplyr::filter(RECIST %in% input$gene_survival_drug_interactive_1_R)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_D)) & !all(is.null(input$gene_survival_drug_interactive_1_L))){
    Data_survival_1 = Data_survival_1 %>% dplyr::filter(Drug %in% input$gene_survival_drug_interactive_1_D & CTx_line %in% input$gene_survival_drug_interactive_1_L)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_AED))){
    AEs <- input$gene_survival_drug_interactive_1_AED
    Data_survival_1 <- Data_survival_1 %>%
      dplyr::filter(if_any(all_of(AEs), ~ . != "AE (-)"))
  }
  if(!all(is.null(input$gene_survival_drug_interactive_1_AEG))){
    Data_survival_1 = Data_survival_1 %>% dplyr::filter(Overall_AE_grage %in% input$gene_survival_drug_interactive_1_AEG)
  }

  ID_2 = unique(Data_drug_TTF_comparison$ID)
  if(!all(is.null(input$gene_survival_drug_interactive_2_P_1))){
    ID_2 = intersect(ID_2, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_2_P_1))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_P_2))){
    ID_2 = intersect(ID_2, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_2_P_2))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_W))){
    ID_2 = setdiff(ID_2, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_2_W))$Tumor_Sample_Barcode)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_W))){
    ID_2 = setdiff(ID_2, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_survival_drug_interactive_2_W))$Tumor_Sample_Barcode)
  }
  if("LUNG" %in% Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
    if(!all(is.null(input$gene_survival_drug_interactive_2_PD_L1))){
      ID_2 = intersect(ID_2, (Data_case_target %>% dplyr::filter(PD_L1 %in% input$gene_survival_drug_interactive_2_PD_L1))$C.CAT調査結果.基本項目.ハッシュID)
    }
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_A))){
    ID_2 = intersect(ID_2, (Data_case_target %>% dplyr::filter(YoungOld %in% input$gene_survival_drug_interactive_2_A))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_S))){
    ID_2 = intersect(ID_2, (Data_case_target %>% dplyr::filter(症例.基本情報.性別.名称. %in% input$gene_survival_drug_interactive_2_S))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_H))){
    ID_2 = intersect(ID_2, (Data_case_target %>% dplyr::filter(症例.基本情報.がん種.OncoTree. %in% input$gene_survival_drug_interactive_2_H))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_C))){
    ID_2 = intersect(ID_2, (Data_case_target %>% dplyr::filter(cluster %in% input$gene_survival_drug_interactive_2_C))$C.CAT調査結果.基本項目.ハッシュID)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_M))){
    mets <- input$gene_survival_drug_interactive_2_M
    ID_M2 <- purrr::map(mets, function(met) {
      Data_case_target %>%
        dplyr::filter(.data[[met]] == "Yes") %>%
        dplyr::pull(C.CAT調査結果.基本項目.ハッシュID)
    }) %>% unlist() %>% unique()
    ID_2 = intersect(ID_2, ID_M2)
  }
  Data_survival_2 = Data_drug_TTF_comparison %>%
    dplyr::filter(ID %in% ID_2)
  Data_survival_2$Group = 2
  if(!all(is.null(input$gene_survival_drug_interactive_2_AE))){
    if("AE (+)" %in% input$gene_survival_drug_interactive_2_AE &
       !"AE (-)" %in% input$gene_survival_drug_interactive_2_AE){
      Data_survival_2 = Data_survival_2 %>% dplyr::filter(AE == 1)
    } else if(!"AE (+)" %in% input$gene_survival_drug_interactive_2_AE &
              "AE (-)" %in% input$gene_survival_drug_interactive_2_AE){
      Data_survival_2 = Data_survival_2 %>% dplyr::filter(AE == 0)
    }
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_EAE))){
    Data_survival_2 = Data_survival_2 %>% dplyr::filter(Early_AE %in% input$gene_survival_drug_interactive_2_EAE)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_D))){
    Data_survival_2 = Data_survival_2 %>% dplyr::filter(Drug %in% input$gene_survival_drug_interactive_2_D)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_R))){
    Data_survival_2 = Data_survival_2 %>% dplyr::filter(RECIST %in% input$gene_survival_drug_interactive_2_R)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_D)) & !all(is.null(input$gene_survival_drug_interactive_2_L))){
    Data_survival_2 = Data_survival_2 %>% dplyr::filter(Drug %in% input$gene_survival_drug_interactive_2_D & CTx_line %in% input$gene_survival_drug_interactive_2_L)
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_AED))){
    AEs <- input$gene_survival_drug_interactive_2_AED
    Data_survival_2 <- Data_survival_2 %>%
      dplyr::filter(if_any(all_of(AEs), ~ . != "AE (-)"))
  }
  if(!all(is.null(input$gene_survival_drug_interactive_2_AEG))){
    Data_survival_2 = Data_survival_2 %>% dplyr::filter(Overall_AE_grage %in% input$gene_survival_drug_interactive_2_AEG)
  }
  Data_survival = rbind(Data_survival_1, Data_survival_2)
  if(nrow(Data_survival) > 0){
    survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~Group, data=Data_survival,type = "kaplan-meier", conf.type = "log-log")
    if(length(unique(Data_survival$Group)) > 1){
      diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~Group,
                        data=Data_survival, rho=0)
      diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~Group,
                        data=Data_survival, rho=1)
      diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~Group,
                     data=Data_survival)
      surv_curv_drug(survfit_t, Data_survival, paste0("Treatment time, Group 1"), diff_0, diff_1, diff_2)
    } else{
      surv_curv_drug(survfit_t, Data_survival, paste0("Treatment time, Group 2"), NULL, NULL, NULL)
    }
  }
})

output$figure_drug_3_forest = renderPlot({
  req(OUTPUT_DATA$drug_analysis_has)
  OUTPUT_DATA$drug_analysis_has
})

output$figure_drug_3_forest_table = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$drug_analysis_gene_table_3)
  create_datatable_with_confirm(OUTPUT_DATA$drug_analysis_gene_table_3, "FELIS downloaded raw data in Treatment time and clinical information tab")
})

output$figure_drug_cluster = renderPlot({
  req(OUTPUT_DATA$drug_analysis_hsc)
  arrange_ggsurvplots(OUTPUT_DATA$drug_analysis_hsc,print=TRUE,ncol=4,nrow=8,surv.plot.height = 0.8,risk.table.height = 0.2)
})
output$figure_drug_cluster_2 = renderPlot({
  req(OUTPUT_DATA$drug_analysis_h_clust)
  OUTPUT_DATA$drug_analysis_h_clust
})

output$figure_drug_cluster_table = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$drug_analysis_gene_table_cluster)
  create_datatable_with_confirm(OUTPUT_DATA$drug_analysis_gene_table_cluster, "FELIS downloaded raw data in Treatment time by gene mutation cluster tab")
})

output$figure_drug_4_table = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$drug_analysis_gene_table_gene)
  create_datatable_with_confirm(OUTPUT_DATA$drug_analysis_gene_table_gene, "FELIS downloaded raw data in Treatment time by mutated genes tab")
})

output$figure_drug_4 = renderPlot({
  OUTPUT_DATA$drug_analysis_h
})

output$figure_drug_5 = renderPlot({
  arrange_ggsurvplots(OUTPUT_DATA$drug_analysis_hs,print=TRUE,ncol=4,nrow=8,surv.plot.height = 0.8,risk.table.height = 0.2)
})

output$figure_drug_pattern = renderPlot({
  req(OUTPUT_DATA$Data_drug_Data_MAF_target,
      OUTPUT_DATA$Data_drug_Data_drug_TTF,
      OUTPUT_DATA$Data_drug_Data_case_target)
  Data_MAF_target = OUTPUT_DATA$Data_drug_Data_MAF_target
  Data_drug_TTF = OUTPUT_DATA$Data_drug_Data_drug_TTF
  Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
  hsa = list()
  mut_gene_ = input$gene[input$gene %in% unique(Data_MAF_target$Hugo_Symbol)]
  if(length(mut_gene_)==0){
    mut_gene_ = names(sort(table(Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
    ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
    Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
      mutation = case_when(
        ID %in% ID_mutation_ ~ paste0(mut_gene_, " mut(+)"),
        TRUE ~ paste0(mut_gene_, " mut(-)")
      )
    )
  } else{
    if(is.null(input$gene_group_2)){
      ID_mutation_ = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% mut_gene_))$Tumor_Sample_Barcode
      Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
        mutation = case_when(
          ID %in% ID_mutation_ ~ paste0(paste(mut_gene_, collapse = "/"), " mut(+)"),
          TRUE ~ paste0(paste(mut_gene_, collapse = "/"), " mut(-)")
        )
      )
    } else{
      ID_mutation_1 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_1))$Tumor_Sample_Barcode
      ID_mutation_2 = (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% input$gene_group_2))$Tumor_Sample_Barcode
      Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
        mutation = case_when(
          ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
          ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(+) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
          !ID %in% ID_mutation_1 & ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(+)"),
          !ID %in% ID_mutation_1 & !ID %in% ID_mutation_2 ~ paste0(paste(input$gene_group_1, collapse = "/"), " mut(-) & ", paste(input$gene_group_2, collapse = "/"), " mut(-)"),
        )
      )
    }
  }
  if(!input$special_gene == ""){
    ID_special_gene = (Data_MAF_target %>%
                         dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_"))))$Tumor_Sample_Barcode
    ID_special_gene_mutation_1 = (Data_MAF_target %>%
                                    dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                    amino.acid.change %in% input$special_gene_mutation_1))$Tumor_Sample_Barcode
    ID_special_gene_mutation_2 = (Data_MAF_target %>%
                                    dplyr::filter(Hugo_Symbol == input$special_gene | str_detect(Hugo_Symbol, paste0(input$special_gene, "_")) &
                                                    amino.acid.change %in% input$special_gene_mutation_2))$Tumor_Sample_Barcode
    Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
      separation_value = case_when(
        ID %in% ID_special_gene_mutation_1 ~ paste0(input$special_gene_mutation_1_name, " in ", input$special_gene),
        ID %in% ID_special_gene_mutation_2 ~ paste0(input$special_gene_mutation_2_name, " in ", input$special_gene),
        ID %in% ID_special_gene ~ paste0("Other mut in ", input$special_gene),
        TRUE ~ paste0("No mut in ", input$special_gene)
      )
    )
  } else{
    Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::mutate(
      separation_value = "No mutation pattern")
  }
  kkk = 1
  if(length(unique(Data_drug_Drug_length_tmp$mutation)) > 1){
    survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~mutation, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~mutation,
                      data=Data_drug_Drug_length_tmp, rho=0)
    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~mutation,
                      data=Data_drug_Drug_length_tmp, rho=1)
    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~mutation,
                   data=Data_drug_Drug_length_tmp)
    hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, mutation, all drugs"), diff_0, diff_1, diff_2)
    kkk = kkk + 1
  }
  if(length(unique(Data_drug_Drug_length_tmp$separation_value)) > 1){
    survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~separation_value, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~separation_value,
                      data=Data_drug_Drug_length_tmp, rho=0)
    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~separation_value,
                      data=Data_drug_Drug_length_tmp, rho=1)
    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~separation_value,
                   data=Data_drug_Drug_length_tmp)
    hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, mutation pattern, all drugs"), diff_0, diff_1, diff_2)
    kkk = kkk + 1
  }
  Data_drug_Drug_length_tmp = Data_drug_Drug_length_tmp %>% dplyr::filter(Drug %in% input$drug)
  if(length(unique(Data_drug_Drug_length_tmp$mutation)) > 1){
    survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~mutation, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~mutation,
                      data=Data_drug_Drug_length_tmp, rho=0)
    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~mutation,
                      data=Data_drug_Drug_length_tmp, rho=1)
    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~mutation,
                   data=Data_drug_Drug_length_tmp)
    hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, mutation, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
    kkk = kkk + 1
  }
  if(length(unique(Data_drug_Drug_length_tmp$separation_value)) > 1){
    survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~separation_value, data=Data_drug_Drug_length_tmp,type = "kaplan-meier", conf.type = "log-log")
    diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~separation_value,
                      data=Data_drug_Drug_length_tmp, rho=0)
    diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~separation_value,
                      data=Data_drug_Drug_length_tmp, rho=1)
    diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~separation_value,
                   data=Data_drug_Drug_length_tmp)
    hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp, paste0("Treatment time, mutation pattern, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
    kkk = kkk + 1
  }
  if(input$HER2 != "No"){
    ID_mutation_P = (Data_case_target %>% dplyr::filter(HER2_IHC == "Positive"))$C.CAT調査結果.基本項目.ハッシュID
    ID_mutation_N = (Data_case_target %>% dplyr::filter(HER2_IHC == "Negative"))$C.CAT調査結果.基本項目.ハッシュID
    Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
      HER2 = case_when(
        ID %in% ID_mutation_P ~ "HER2 IHC (+)",
        ID %in% ID_mutation_N ~ "HER2 IHC (-)",
        TRUE ~ "Unknown"
      )
    )
    Data_drug_Drug_length_tmp_ = Data_drug_Drug_length_tmp %>%
      dplyr::filter(HER2 != "Unknown")
    if(length(unique(Data_drug_Drug_length_tmp_$HER2)) > 1){
      survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~HER2, data=Data_drug_Drug_length_tmp_,type = "kaplan-meier", conf.type = "log-log")
      diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~HER2,
                        data=Data_drug_Drug_length_tmp_, rho=0)
      diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~HER2,
                        data=Data_drug_Drug_length_tmp_, rho=1)
      diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~HER2,
                     data=Data_drug_Drug_length_tmp_)
      hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp_, paste0("Treatment time, mutation pattern, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
      kkk = kkk + 1
    }
  }
  if(input$MSI != "No"){
    ID_mutation_P = (Data_case_target %>% dplyr::filter(MSI_PCR == "dMMR"))$C.CAT調査結果.基本項目.ハッシュID
    ID_mutation_N = (Data_case_target %>% dplyr::filter(MSI_PCR == "pMMR"))$C.CAT調査結果.基本項目.ハッシュID
    Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
      MSI = case_when(
        ID %in% ID_mutation_P ~ "MSI PCR dMMR",
        ID %in% ID_mutation_N ~ "MSI PCR pMMR",
        TRUE ~ "Unknown"
      )
    )
    Data_drug_Drug_length_tmp_ = Data_drug_Drug_length_tmp %>%
      dplyr::filter(MSI != "Unknown")
    if(length(unique(Data_drug_Drug_length_tmp_$MSI)) > 1){
      survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~HER2, data=Data_drug_Drug_length_tmp_,type = "kaplan-meier", conf.type = "log-log")
      diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~MSI,
                        data=Data_drug_Drug_length_tmp_, rho=0)
      diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~MSI,
                        data=Data_drug_Drug_length_tmp_, rho=1)
      diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~MSI,
                     data=Data_drug_Drug_length_tmp_)
      hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp_, paste0("Treatment time, mutation pattern, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
      kkk = kkk + 1
    }
  }
  if(input$MMR != "No"){
    ID_mutation_P = (Data_case_target %>% dplyr::filter(MMR_IHC == "Positive"))$C.CAT調査結果.基本項目.ハッシュID
    ID_mutation_N = (Data_case_target %>% dplyr::filter(MMR_IHC == "Negative"))$C.CAT調査結果.基本項目.ハッシュID
    Data_drug_Drug_length_tmp = Data_drug_TTF %>% dplyr::mutate(
      MMR = case_when(
        ID %in% ID_mutation_P ~ "MMR IHC (+)",
        ID %in% ID_mutation_N ~ "MMR IHC (-)",
        TRUE ~ "Unknown"
      )
    )
    Data_drug_Drug_length_tmp_ = Data_drug_Drug_length_tmp %>%
      dplyr::filter(MMR != "Unknown")
    if(length(unique(Data_drug_Drug_length_tmp_$MMR)) > 1){
      survfit_t <- survfit(Surv(Drug_length, Drug_length_censor)~HER2, data=Data_drug_Drug_length_tmp_,type = "kaplan-meier", conf.type = "log-log")
      diff_0 = survdiff(Surv(Drug_length, Drug_length_censor)~MMR,
                        data=Data_drug_Drug_length_tmp_, rho=0)
      diff_1 = survdiff(Surv(Drug_length, Drug_length_censor)~MMR,
                        data=Data_drug_Drug_length_tmp_, rho=1)
      diff_2 = coxph(Surv(Drug_length, Drug_length_censor)~MMR,
                     data=Data_drug_Drug_length_tmp_)
      hsa[[kkk]] = surv_curv_drug(survfit_t, Data_drug_Drug_length_tmp_, paste0("Treatment time, mutation pattern, treated with ", paste(input$drug, collapse = ";")), diff_0, diff_1, diff_2)
      kkk = kkk + 1
    }
  }
  if(kkk<=18){
    for(kkkk in kkk:18){
      hsa[[kkkk]] = ggsurvplot_empty()
    }
  }
  arrange_ggsurvplots(hsa,print=TRUE,ncol=min(3, length(hsa)), nrow=ceiling(length(hsa)/3),surv.plot.height = 0.8,risk.table.height = 0.2)
})

output$figure_CI_AE = renderPlot({
  req(OUTPUT_DATA$Data_drug_Data_drug_original,
      OUTPUT_DATA$Data_drug_Data_case_target,
      input$figure_CI_AE_var)
  Data_drug = OUTPUT_DATA$Data_drug_Data_drug_original
  Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
  Data_drug$ToT = Data_drug$ToT / 365.25 * 12
  Data_drug$TTAE = Data_drug$TTAE / 365.25 * 12
  Data_drug$TTD = Data_drug$TTD / 365.25 * 12

  max_samples = isolate(input$AE_max_samples)
  max_samples = ifelse(is.null(max_samples), 5000, max_samples)

  Data_drug_analysis = Data_drug %>%
    dplyr::filter(Drug %in% input$figure_CI_AE_var) %>%
    dplyr::mutate(
      ORR_data = case_when(
        最良総合効果 %in% c("CR", "PR") ~ 1,
        最良総合効果 %in% c("LongSD", "SD", "PD") ~ 0,
        最良総合効果 %in% c("NE", "Unknown") ~ -1,
        TRUE ~ NA
      ), TTAE_censor = case_when(
        is.na(TTAE) ~ 0,
        TRUE ~ 1
      )
    ) %>% dplyr::filter(TTD>0 & !is.na(TTD) & is.finite(TTD) &
                          ToT>0 & !is.na(ToT) & is.finite(ToT) &
                          ORR_data != -1)
  Data_case_age_sex_diagnosis = Data_case_target %>%
    dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                  YoungOld,
                  症例.基本情報.性別.名称.,
                  症例.背景情報.喫煙歴有無.名称.,
                  症例.基本情報.がん種.OncoTree..名称.
    ) %>% distinct()
  colnames(Data_case_age_sex_diagnosis) = c(
    "ID",
    "Age",
    "Sex",
    "Smoking",
    "Diagnosis"
  )
  Data_drug_analysis = Data_drug_analysis %>%
    dplyr::left_join(Data_case_age_sex_diagnosis, by = "ID")

  Data_drug_analysis <- Data_drug_analysis %>%
    mutate(
      # 各症例の観察終了時間を決定
      time_to_event = case_when(
        # 有害事象が発生：その時間まで観察
        TTAE_censor == 1 & !is.na(TTAE) ~ TTAE,
        TRUE ~ pmin(ToT + 90, TTD, na.rm = TRUE)
      ),
      # イベントタイプ（時間の整合性も確認）
      event_type = case_when(
        # 有害事象発生（観察期間内）
        TTAE_censor == 1 & !is.na(TTAE) & TTAE <= time_to_event ~ 1,
        # 投薬完了（有害事象なし）
        ToT_censor == 1 & (TTAE_censor == 0 | is.na(TTAE) | TTAE > time_to_event) ~ 2,
        # 打ち切り
        TRUE ~ 0
      ),
      # 投薬期間中か終了後かの詳細分類
      AE_timing_detailed = case_when(
        TTAE_censor == 1 & !is.na(TTAE) & !is.na(ToT) & TTAE <= ToT ~ "During_treatment",
        TTAE_censor == 1 & !is.na(TTAE) & !is.na(ToT) & TTAE > ToT ~ "Post_treatment",
        TTAE_censor == 1 & !is.na(TTAE) & is.na(ToT) ~ "Unknown_timing",
        TRUE ~ "No_AE"
      )
    ) %>%
    filter(time_to_event > 0 & !is.na(time_to_event))
  if (nrow(Data_drug_analysis) > max_samples){
    sample_indices = sample(nrow(Data_drug_analysis),max_samples)
    Data_drug_analysis = Data_drug_analysis[sample_indices, ]
    ENTIRE = "subsampled"
  } else {
    Data_drug_analysis = Data_drug_analysis
    ENTIRE = "whole"
  }


  # 多変量解析：共変量を含む
  covariates_matrix <- model.matrix(
    ~ ORR_data + Age + Sex + Smoking + 治療ライン + Diagnosis,
    data = Data_drug_analysis
  )[, -1]  # 切片を除去
  crr_multivariate <- crr(
    ftime = Data_drug_analysis$time_to_event,
    fstatus = Data_drug_analysis$event_type,
    cov1 = covariates_matrix,
    failcode = 1,
    cencode = 0
  )
  summary(crr_multivariate)
  # 結果の整理
  crr_results <- data.frame(
    HR = exp(crr_multivariate$coef),
    CI_lower = exp(crr_multivariate$coef - 1.96 * sqrt(diag(crr_multivariate$var))),
    CI_upper = exp(crr_multivariate$coef + 1.96 * sqrt(diag(crr_multivariate$var))),
    p_value = 2 * (1 - pnorm(abs(crr_multivariate$coef / sqrt(diag(crr_multivariate$var)))))
  )
  print(crr_results)
  # 累積発生率の計算
  cif_result <- cuminc(
    ftime = Data_drug_analysis$time_to_event,
    fstatus = Data_drug_analysis$event_type,
    group = Data_drug_analysis$ORR_data,
    cencode = 0
  )
  # Gray's test
  print(cif_result$Tests)
  # CIFデータをggplot用に変換
  cif_data <- list()
  group_names <- names(cif_result)[!names(cif_result) %in% "Tests"]
  for(i in seq_along(group_names)) {
    group_name <- group_names[i]
    parts <- str_split(group_name, " ")[[1]]
    effectiveness <- parts[1]
    event_type <- parts[2]
    cif_data[[i]] <- data.frame(
      time = cif_result[[group_name]]$time,
      est = cif_result[[group_name]]$est,
      var = cif_result[[group_name]]$var,
      effectiveness = factor(effectiveness,
                             levels = c("0", "1"),
                             labels = c("No Response", "Response")),
      event_type = factor(event_type,
                          levels = c("1", "2"),
                          labels = c("Adverse Event", "Treatment Complete"))
    )
  }
  cif_df <- bind_rows(cif_data) %>%
    mutate(
      ci_lower = pmax(0, est - 1.96 * sqrt(var)),
      ci_upper = pmin(1, est + 1.96 * sqrt(var))
    )
  # 有害事象のみに焦点を当てた可視化
  adverse_only <- cif_df %>%
    filter(event_type == "Adverse Event")
  p2 <- ggplot(adverse_only, aes(x = time, y = est, color = effectiveness)) +
    geom_step(size = 1.5) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = effectiveness),
                alpha = 0.3, linetype = 0) +
    scale_color_manual(values = c("No Response" = "#E31A1C", "Response" = "#1F78B4")) +
    scale_fill_manual(values = c("No Response" = "#E31A1C", "Response" = "#1F78B4")) +
    labs(
      title = paste0("Cumulative Incidence of Adverse Events, ", nrow(Data_drug_analysis), " ", ENTIRE, " cources"),
      subtitle = paste0("Accounting for Competing Risk of Treatment Completion\nDrugs: ",
                        paste0(
                          "Including ",
                          paste(head(input$figure_CI_AE_var, 2), collapse = ";"),
                          if (length(input$figure_CI_AE_var) > 2) "..." else ""
                        )),
      x = "Months from Treatment Start",
      y = "Cumulative Incidence of Adverse Events",
      color = "Drug Effectiveness",
      fill = "Drug Effectiveness"
    ) +
    theme_classic() +
    scale_y_continuous(labels = percent) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid.minor = element_blank()
    )
  # 統計結果を追加
  hr_text <- paste("Adjusted HR:",
                   round(exp(crr_multivariate$coef[1]), 2),
                   "\n95% CI:",
                   round(exp(crr_multivariate$coef[1] - 1.96 * sqrt(crr_multivariate$var[1,1])), 2), "-",
                   round(exp(crr_multivariate$coef[1] + 1.96 * sqrt(crr_multivariate$var[1,1])), 2),
                   "\np =", round(crr_results$p_value[1], 3))
  p2_final <- p2 +
    annotate("text",
             x = max(adverse_only$time) * 0.6,
             y = max(adverse_only$est) * 0.1,
             label = hr_text,
             size = 4,
             hjust = 0,
             vjust = 1,
             box.margin = unit(0.5, "lines"))
  print(p2_final)
})
