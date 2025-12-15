volcano_AE_logic <- function() {
  volcano_env <- new.env()
  with(volcano_env, {
    Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
    Data_drug_RECIST = OUTPUT_DATA$drug_analysis_Data_drug_RECIST
    Data_MAF_target_tmp = OUTPUT_DATA$Data_drug_Data_MAF_target_tmp
    significant_genes_AE = NULL

    Data_drug_volcano_all = Data_drug_RECIST %>%
      dplyr::filter(治療ライン %in% input$target_line) %>%
      dplyr::mutate(Adverse_effect = case_when(
        Adverse_effect == "AE (+)" ~ "1",
        TRUE ~ "0"
      )) %>%
      dplyr::arrange(ID, Adverse_effect) %>%
      dplyr::select(ID, Drug, Adverse_effect)
    Data_drug_volcano_all$Adverse_effect = as.integer(Data_drug_volcano_all$Adverse_effect)
    drugs_volcano = sort(table(Data_drug_volcano_all$Drug),decreasing = T)
    drugs_volcano = drugs_volcano[!names(drugs_volcano) %in% c("Others", "Unknown", "")]
    drugs_volcano = drugs_volcano[!is.na(names(drugs_volcano))]
    drugs_volcano = drugs_volcano[1:min(length(drugs_volcano), 199)]
    df_volcano_AE = data.frame(Regimen = as.character(),
                               Gene = as.character(),
                               Total_treated_patients = as.numeric(),
                               Mut_respond_patients = as.numeric(),
                               WT_respond_patients = as.numeric(),
                               Mut_nonrespond_patients = as.numeric(),
                               WT_nonrespond_patients = as.numeric(),
                               OddsRatio_considering_pathology = as.numeric(),
                               p_value_logistic_regression = as.numeric()
    )
    df_volcano_tmp = data.frame(Regimen = "",
                                Gene = "",
                                Total_treated_patients = 0,
                                Mut_respond_patients = 0,
                                WT_respond_patients = 0,
                                Mut_nonrespond_patients = 0,
                                WT_nonrespond_patients = 0,
                                OddsRatio_considering_pathology = 1,
                                p_value_logistic_regression = 1
    )
    figure_no = 1
    g_volcano_AE = list()
    drug_name_list_AE = NULL
    options(warn = -1)
    withProgress(message = "Volcano plot analysis for adverse effect", {
      for(drug_name in names(drugs_volcano)){
        Data_drug_volcano = Data_drug_volcano_all %>%
          dplyr::filter(Drug == drug_name) %>%
          dplyr::select(ID, Adverse_effect)
        if(nrow(Data_drug_volcano) >= 1){
          Data_MAF_target_volcano = Data_MAF_target_tmp %>%
            dplyr::filter(Tumor_Sample_Barcode %in%
                            Data_drug_volcano$ID) %>%
            dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
          Data_case_target_volcano = Data_case_target %>%
            dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                            Data_drug_volcano$ID) %>%
            dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.がん種.OncoTree..名称.)
          
          colnames(Data_drug_volcano) = c("ID", "Response")
          colnames(Data_MAF_target_volcano) = c("ID", "Gene")
          colnames(Data_case_target_volcano) = c("ID", "Diagnosis")
          Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
          Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_case_target_volcano, by="ID")
          
          df_volcano_tmp$Regimen = drug_name
          df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
          gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
          gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
          if(length(gene_table_volcano) > 0){
            for(gene_volcano in names(gene_table_volcano)){
              Data_drug_volcano_pos = Data_MAF_target_volcano %>%
                dplyr::filter(Gene == gene_volcano) %>%
                dplyr::select(-Gene)
              Data_drug_volcano_neg = Data_MAF_target_volcano %>%
                dplyr::select(-Gene) %>%
                dplyr::distinct()
              Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
              if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_pos$mutation = 1
              }
              if(nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg$mutation = 0
              }
              if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
              } else if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_neg_pos = Data_drug_volcano_pos
              } else{
                Data_drug_volcano_neg_pos = Data_drug_volcano_neg
              }
              df_volcano_tmp$Gene = gene_volcano
              df_volcano_tmp$Mut_respond_patients = sum(Data_drug_volcano_pos$Response)
              df_volcano_tmp$Mut_nonrespond_patients = nrow(Data_drug_volcano_pos) - df_volcano_tmp$Mut_respond_patients
              df_volcano_tmp$WT_respond_patients = sum(Data_drug_volcano_neg$Response)
              df_volcano_tmp$WT_nonrespond_patients = nrow(Data_drug_volcano_neg) - df_volcano_tmp$WT_respond_patients
              Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
              Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$Response)>1 &
                 sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else {
                major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
                Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                  dplyr::filter(Diagnosis %in% major_pathology)
                if(input$volcano_pathology == "Yes" &
                   min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                   max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                   sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$Response)>1 &
                   sum(Data_drug_volcano_neg_pos_major$Response)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                  Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos_major, family=binomial)
                  df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                  df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
                } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                          sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                          sum(Data_drug_volcano_neg_pos$Response)>1 &
                          sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                  Data_ORR = glm(Response ~ mutation, data=Data_drug_volcano_neg_pos, family=binomial)
                  df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                  df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
                } else {
                  df_volcano_tmp$OddsRatio_considering_pathology = 1
                  df_volcano_tmp$p_value_logistic_regression = 1
                }
                if((df_volcano_tmp$OddsRatio_considering_pathology > 16 | df_volcano_tmp$OddsRatio_considering_pathology < 1/16) &
                   df_volcano_tmp$Mut_respond_patients == 0 |
                   df_volcano_tmp$Mut_nonrespond_patients == 0 |
                   df_volcano_tmp$WT_respond_patients == 0 |
                   df_volcano_tmp$WT_nonrespond_patients == 0){
                  fisher_tmp = matrix(c(df_volcano_tmp$Mut_respond_patients,
                                        df_volcano_tmp$Mut_nonrespond_patients,
                                        df_volcano_tmp$WT_respond_patients,
                                        df_volcano_tmp$WT_nonrespond_patients),
                                      ncol = 2
                  )
                  df_volcano_tmp$p_value_logistic_regression = fisher.test(fisher_tmp)$p.value
                }
              }
              df_volcano_AE = rbind(df_volcano_AE, df_volcano_tmp)
            }
            
            if(figure_no<200){
              Figure_volcano_tmp = df_volcano_AE %>%
                dplyr::filter(Regimen == drug_name) %>%
                dplyr::mutate(Gene_effect = ifelse(p_value_logistic_regression >= 0.05 , "Not significant",
                                                   ifelse(OddsRatio_considering_pathology >= 2, "Relevant",
                                                          ifelse(OddsRatio_considering_pathology <= 0.5, "Irrelevant","Not significant"))))
              Figure_volcano_tmp$OddsRatio_considering_pathology[Figure_volcano_tmp$OddsRatio_considering_pathology>16]=16
              Figure_volcano_tmp$OddsRatio_considering_pathology[Figure_volcano_tmp$OddsRatio_considering_pathology<1/16]=1/16
              Figure_volcano_tmp$Gene_effect = factor(Figure_volcano_tmp$Gene_effect, levels = c("Irrelevant", "Not significant", "Relevant"))
              Figure_volcano_tmp$label = ifelse((Figure_volcano_tmp$OddsRatio_considering_pathology >=2 |
                                                   Figure_volcano_tmp$p_value_logistic_regression < 0.05 |
                                                   Figure_volcano_tmp$OddsRatio_considering_pathology <= 0.5), Figure_volcano_tmp$Gene, NA_character_)
              
              Figure_volcano_tmp$log2OR = log2(Figure_volcano_tmp$OddsRatio_considering_pathology)
              Figure_volcano_tmp$log10pval = -log10(Figure_volcano_tmp$p_value_logistic_regression)
              if(nrow(Figure_volcano_tmp)>0){
                g_volcano_AE[[figure_no]] =
                  ggplot(data = Figure_volcano_tmp, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
                  geom_point(size = 2) +
                  labs(x = expression("log"[2]*"odds ratio of adverse effect"), y = expression("-log"[10]*"p-value")) +
                  scale_color_manual(values = c("#bb0c00", "grey", "#00AFBB"), 
                                     labels = c("Irrelevant", "Not significant", "Relevant"),
                                     drop = FALSE) +
                  geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
                  geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  theme(legend.position = "none", panel.grid = element_blank()) +
                  geom_text_repel(max.overlaps = Inf) +
                  ggtitle(paste0(drug_name, " ",
                                 df_volcano_tmp$Total_treated_patients, " patients")) +
                  theme_classic()
              }
              figure_no = figure_no + 1
              drug_name_list_AE = c(drug_name_list_AE, drug_name)
            }
          }
        }
        incProgress(1 / length(names(drugs_volcano)))
      }
    })
    options(warn = 0)
    if(!is.null(input$drug)){
      Data_drug_volcano = Data_drug_volcano_all %>%
        dplyr::filter(Drug %in% input$drug) %>%
        dplyr::select(ID, Adverse_effect)
      if(nrow(Data_drug_volcano) >= 1){
        Data_MAF_target_volcano = Data_MAF_target_tmp %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_drug_volcano$ID) %>%
          dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
        Data_case_target_volcano = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                          Data_drug_volcano$ID) %>%
          dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.がん種.OncoTree..名称.)
        
        colnames(Data_drug_volcano) = c("ID", "Response")
        colnames(Data_MAF_target_volcano) = c("ID", "Gene")
        colnames(Data_case_target_volcano) = c("ID", "Diagnosis")
        Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
        Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_case_target_volcano, by="ID")
        
        df_volcano_tmp$Regimen = paste(input$drug, collapse = ";")
        df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
        gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
        gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
        if(length(gene_table_volcano) > 0){
          for(gene_volcano in names(gene_table_volcano)){
            Data_drug_volcano_pos = Data_MAF_target_volcano %>%
              dplyr::filter(Gene == gene_volcano) %>%
              dplyr::select(-Gene)
            Data_drug_volcano_neg = Data_MAF_target_volcano %>%
              dplyr::select(-Gene) %>%
              dplyr::distinct()
            Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
            if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_pos$mutation = 1
            }
            if(nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg$mutation = 0
            }
            if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
            } else if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_neg_pos = Data_drug_volcano_pos
            } else {
              Data_drug_volcano_neg_pos = Data_drug_volcano_neg
            }
            df_volcano_tmp$Gene = gene_volcano
            df_volcano_tmp$Mut_respond_patients = sum(Data_drug_volcano_pos$Response)
            df_volcano_tmp$Mut_nonrespond_patients = nrow(Data_drug_volcano_pos) - df_volcano_tmp$Mut_respond_patients
            df_volcano_tmp$WT_respond_patients = sum(Data_drug_volcano_neg$Response)
            df_volcano_tmp$WT_nonrespond_patients = nrow(Data_drug_volcano_neg) - df_volcano_tmp$WT_respond_patients
            Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
            Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
            if(input$volcano_pathology == "Yes" &
               min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
               max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$mutation)>1 &
               sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$Response)>1 &
               sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
              Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos, family=binomial)
              df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
              df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
            } else {
              major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
              Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                dplyr::filter(Diagnosis %in% major_pathology)
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$Response)>1 &
                 sum(Data_drug_volcano_neg_pos_major$Response)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos_major, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                        sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                        sum(Data_drug_volcano_neg_pos$Response)>1 &
                        sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = glm(Response ~ mutation, data=Data_drug_volcano_neg_pos, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else {
                df_volcano_tmp$OddsRatio_considering_pathology = 1
                df_volcano_tmp$p_value_logistic_regression = 1
              }
            }
            df_volcano_AE = rbind(df_volcano_AE, df_volcano_tmp)
          }
          Figure_volcano = df_volcano_AE %>%
            dplyr::filter(Regimen == paste(input$drug, collapse = ";")) %>%
            dplyr::mutate(Gene_effect = ifelse(p_value_logistic_regression >= 0.05 , "Not significant",
                                               ifelse(OddsRatio_considering_pathology >= 2, "Relevant",
                                                      ifelse(OddsRatio_considering_pathology <= 0.5, "Irrelevant","Not significant"))))
          Figure_volcano$label = ifelse((Figure_volcano$OddsRatio_considering_pathology >=2 |
                                           Figure_volcano$p_value_logistic_regression < 0.05 |
                                           Figure_volcano$OddsRatio_considering_pathology <= 0.5), Figure_volcano$Gene, NA_character_)
          Figure_volcano$Gene_effect = factor(Figure_volcano$Gene_effect, levels = c("Irrelevant", "Not significant", "Relevant"))
          Figure_volcano$OddsRatio_considering_pathology[Figure_volcano$OddsRatio_considering_pathology>16]=16
          Figure_volcano$OddsRatio_considering_pathology[Figure_volcano$OddsRatio_considering_pathology<1/16]=1/16
          
          Figure_volcano$log2OR = log2(Figure_volcano$OddsRatio_considering_pathology)
          Figure_volcano$log10pval = -log10(Figure_volcano$p_value_logistic_regression)
          if(nrow(Figure_volcano)>0){
            g_volcano_AE[[figure_no]] =
              ggplot(data = Figure_volcano, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
              geom_point(size = 2) +
              labs(x = expression("log"[2]*"odds ratio of adverse effect"), y = expression("-log"[10]*"p-value")) +
              scale_color_manual(values = c("#bb0c00", "grey", "#00AFBB"), 
                                 labels = c("Irrelevant", "Not significant", "Relevant"),
                                 drop = FALSE) +
              geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
              geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
              geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
              theme(legend.position = "none", panel.grid = element_blank()) +
              ggtitle(paste0(paste(input$drug, collapse = ";"), " ",
                             Figure_volcano$Total_treated_patients[[1]], " patients")) +
              geom_text_repel(max.overlaps = Inf) +
              theme_classic()
            figure_no = figure_no + 1
            drug_name_list_AE = c(drug_name_list_AE, paste(input$drug, collapse = ";"))
          }
          significant_genes_AE = c(significant_genes_AE, (Figure_volcano %>% dplyr::arrange(p_value_logistic_regression))$label)
          significant_genes_AE = significant_genes_AE[!is.na(significant_genes_AE)]
          
        }
      }
    }
    
    if(length(drug_name_list_AE) > 0){
      display_names_AE <- paste0(drug_name_list_AE, " (", seq_along(drug_name_list_AE), ")")
      regimen_choice_RECIST_AE = setNames(as.character(seq_along(drug_name_list_AE)), display_names_AE)
    } else {
      regimen_choice_RECIST_AE = NULL
    }
    OUTPUT_DATA$drug_analysis_regimen_choice_RECIST_AE = regimen_choice_RECIST_AE
    OUTPUT_DATA$drug_analysis_g_volcano_AE = g_volcano_AE
    OUTPUT_DATA$drug_analysis_df_volcano_AE = df_volcano_AE
    OUTPUT_DATA$drug_analysis_significant_genes_AE = significant_genes_AE
  })
  rm(volcano_env)
  gc()
}

volcano_ORR_logic <- function() {
  volcano_env <- new.env()
  with(volcano_env, {
    Data_case_target = OUTPUT_DATA$Data_drug_Data_case_target
    Data_drug_RECIST = OUTPUT_DATA$drug_analysis_Data_drug_RECIST
    Data_MAF_target_tmp = OUTPUT_DATA$Data_drug_Data_MAF_target_tmp
    significant_genes = NULL
    
    Data_drug_volcano_all = Data_drug_RECIST %>%
      dplyr::filter(治療ライン %in% input$target_line) %>%
      dplyr::filter(!最良総合効果 %in% c("NE", "Unknown")) %>%
      dplyr::mutate(最良総合効果 = case_when(
        最良総合効果 %in% c("CR", "PR") ~ "1",
        TRUE ~ "0"
      )) %>%
      dplyr::arrange(ID, 治療ライン) %>%
      dplyr::select(ID, Drug, 最良総合効果)
    Data_drug_volcano_all$最良総合効果 = as.integer(Data_drug_volcano_all$最良総合効果)
    drugs_volcano = sort(table(Data_drug_volcano_all$Drug),decreasing = T)
    drugs_volcano = drugs_volcano[!names(drugs_volcano) %in% c("Others", "Unknown", "")]
    drugs_volcano = drugs_volcano[!is.na(names(drugs_volcano))]
    drugs_volcano = drugs_volcano[1:min(length(drugs_volcano), 199)]
    df_volcano = data.frame(Regimen = as.character(),
                            Gene = as.character(),
                            Total_treated_patients = as.numeric(),
                            Mut_respond_patients = as.numeric(),
                            WT_respond_patients = as.numeric(),
                            Mut_nonrespond_patients = as.numeric(),
                            WT_nonrespond_patients = as.numeric(),
                            OddsRatio_considering_pathology = as.numeric(),
                            p_value_logistic_regression = as.numeric()
    )
    df_volcano_tmp = data.frame(Regimen = "",
                                Gene = "",
                                Total_treated_patients = 0,
                                Mut_respond_patients = 0,
                                WT_respond_patients = 0,
                                Mut_nonrespond_patients = 0,
                                WT_nonrespond_patients = 0,
                                OddsRatio_considering_pathology = 1,
                                p_value_logistic_regression = 1
    )
    figure_no = 1
    g_volcano = list()
    drug_name_list = NULL
    options(warn = -1)
    withProgress(message = "Volcano plot analysis for objective response", {
      for(drug_name in names(drugs_volcano)){
        Data_drug_volcano = Data_drug_volcano_all %>%
          dplyr::filter(Drug == drug_name) %>%
          dplyr::select(ID, 最良総合効果)
        if(nrow(Data_drug_volcano) >= 1){
          Data_MAF_target_volcano = Data_MAF_target_tmp %>%
            dplyr::filter(Tumor_Sample_Barcode %in%
                            Data_drug_volcano$ID) %>%
            dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
          Data_case_target_volcano = Data_case_target %>%
            dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                            Data_drug_volcano$ID) %>%
            dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.がん種.OncoTree..名称.)
          
          colnames(Data_drug_volcano) = c("ID", "Response")
          colnames(Data_MAF_target_volcano) = c("ID", "Gene")
          colnames(Data_case_target_volcano) = c("ID", "Diagnosis")
          Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
          Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_case_target_volcano, by="ID")
          
          df_volcano_tmp$Regimen = drug_name
          df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
          gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
          gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
          if(length(gene_table_volcano) > 0){
            for(gene_volcano in names(gene_table_volcano)){
              Data_drug_volcano_pos = Data_MAF_target_volcano %>%
                dplyr::filter(Gene == gene_volcano) %>%
                dplyr::select(-Gene)
              Data_drug_volcano_neg = Data_MAF_target_volcano %>%
                dplyr::select(-Gene) %>%
                dplyr::distinct()
              Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
              if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_pos$mutation = 1
              }
              if(nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg$mutation = 0
              }
              if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
              } else if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_neg_pos = Data_drug_volcano_pos
              } else{
                Data_drug_volcano_neg_pos = Data_drug_volcano_neg
              }
              df_volcano_tmp$Gene = gene_volcano
              df_volcano_tmp$Mut_respond_patients = sum(Data_drug_volcano_pos$Response)
              df_volcano_tmp$Mut_nonrespond_patients = nrow(Data_drug_volcano_pos) - df_volcano_tmp$Mut_respond_patients
              df_volcano_tmp$WT_respond_patients = sum(Data_drug_volcano_neg$Response)
              df_volcano_tmp$WT_nonrespond_patients = nrow(Data_drug_volcano_neg) - df_volcano_tmp$WT_respond_patients
              # print(gene_volcano)
              # print(drug_name)
              Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
              Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$Response)>1 &
                 sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else {
                major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
                Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                  dplyr::filter(Diagnosis %in% major_pathology)
                if(input$volcano_pathology == "Yes" &
                   min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                   max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                   sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$Response)>1 &
                   sum(Data_drug_volcano_neg_pos_major$Response)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                  Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos_major, family=binomial)
                  df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                  df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
                } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                          sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                          sum(Data_drug_volcano_neg_pos$Response)>1 &
                          sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                  Data_ORR = glm(Response ~ mutation, data=Data_drug_volcano_neg_pos, family=binomial)
                  df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                  df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
                } else {
                  df_volcano_tmp$OddsRatio_considering_pathology = 1
                  df_volcano_tmp$p_value_logistic_regression = 1
                }
                if((df_volcano_tmp$OddsRatio_considering_pathology > 16 | df_volcano_tmp$OddsRatio_considering_pathology < 1/16) &
                   df_volcano_tmp$Mut_respond_patients == 0 |
                   df_volcano_tmp$Mut_nonrespond_patients == 0 |
                   df_volcano_tmp$WT_respond_patients == 0 |
                   df_volcano_tmp$WT_nonrespond_patients == 0){
                  fisher_tmp = matrix(c(df_volcano_tmp$Mut_respond_patients,
                                        df_volcano_tmp$Mut_nonrespond_patients,
                                        df_volcano_tmp$WT_respond_patients,
                                        df_volcano_tmp$WT_nonrespond_patients),
                                      ncol = 2
                  )
                  df_volcano_tmp$p_value_logistic_regression = fisher.test(fisher_tmp)$p.value
                }
              }
              df_volcano = rbind(df_volcano, df_volcano_tmp)
            }
            
            if(figure_no<200){
              Figure_volcano_tmp = df_volcano %>%
                dplyr::filter(Regimen == drug_name) %>%
                dplyr::mutate(Gene_effect = ifelse(p_value_logistic_regression >= 0.05 , "Not significant",
                                                   ifelse(OddsRatio_considering_pathology >= 2, "Effective",
                                                          ifelse(OddsRatio_considering_pathology <= 0.5, "Not effective","Not significant"))))
              Figure_volcano_tmp$OddsRatio_considering_pathology[Figure_volcano_tmp$OddsRatio_considering_pathology>16]=16
              Figure_volcano_tmp$OddsRatio_considering_pathology[Figure_volcano_tmp$OddsRatio_considering_pathology<1/16]=1/16
              Figure_volcano_tmp$Gene_effect = factor(Figure_volcano_tmp$Gene_effect, levels = c("Not effective", "Not significant", "Effective"))
              # Figure_volcano_tmp$label = ifelse((Figure_volcano_tmp$OddsRatio_considering_pathology >=2 &
              #                                      Figure_volcano_tmp$p_value_logistic_regression < 0.05) |
              #                             (Figure_volcano_tmp$OddsRatio_considering_pathology <= 0.5 &
              #                                Figure_volcano_tmp$p_value_logistic_regression < 0.05), Figure_volcano_tmp$Gene, NA_character_)
              Figure_volcano_tmp$label = ifelse((Figure_volcano_tmp$OddsRatio_considering_pathology >=2 |
                                                   Figure_volcano_tmp$p_value_logistic_regression < 0.05 |
                                                   Figure_volcano_tmp$OddsRatio_considering_pathology <= 0.5), Figure_volcano_tmp$Gene, NA_character_)
              
              Figure_volcano_tmp$log2OR = log2(Figure_volcano_tmp$OddsRatio_considering_pathology)
              Figure_volcano_tmp$log10pval = -log10(Figure_volcano_tmp$p_value_logistic_regression)
              if(nrow(Figure_volcano_tmp)>0){
                g_volcano[[figure_no]] =
                  ggplot(data = Figure_volcano_tmp, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
                  geom_point(size = 2) +
                  labs(x = expression("log"[2]*"odds ratio of objective response"), y = expression("-log"[10]*"p-value")) +
                  scale_color_manual(values = c("#bb0c00", "grey", "#00AFBB"), 
                                     labels = c("Not effective", "Not significant", "Effective"),
                                     drop = FALSE) +
                  geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
                  geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  theme(legend.position = "none", panel.grid = element_blank()) +
                  geom_text_repel(max.overlaps = Inf) +
                  ggtitle(paste0(drug_name, " ",
                                 df_volcano_tmp$Total_treated_patients, " patients")) +
                  theme_classic()
              }
              figure_no = figure_no + 1
              drug_name_list = c(drug_name_list, drug_name)
            }
          }
        }
        incProgress(1 / length(names(drugs_volcano)))
      }
    })
    options(warn = 0)
    if(!is.null(input$drug)){
      Data_drug_volcano = Data_drug_volcano_all %>%
        dplyr::filter(Drug %in% input$drug) %>%
        dplyr::select(ID, 最良総合効果)
      if(nrow(Data_drug_volcano) >= 1){
        Data_MAF_target_volcano = Data_MAF_target_tmp %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_drug_volcano$ID) %>%
          dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
        Data_case_target_volcano = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                          Data_drug_volcano$ID) %>%
          dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.がん種.OncoTree..名称.)
        
        colnames(Data_drug_volcano) = c("ID", "Response")
        colnames(Data_MAF_target_volcano) = c("ID", "Gene")
        colnames(Data_case_target_volcano) = c("ID", "Diagnosis")
        Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
        Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_case_target_volcano, by="ID")
        
        df_volcano_tmp$Regimen = paste(input$drug, collapse = ";")
        df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
        gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
        gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
        if(length(gene_table_volcano) > 0){
          for(gene_volcano in names(gene_table_volcano)){
            Data_drug_volcano_pos = Data_MAF_target_volcano %>%
              dplyr::filter(Gene == gene_volcano) %>%
              dplyr::select(-Gene)
            Data_drug_volcano_neg = Data_MAF_target_volcano %>%
              dplyr::select(-Gene) %>%
              dplyr::distinct()
            Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
            if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_pos$mutation = 1
            }
            if(nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg$mutation = 0
            }
            if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
            } else if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_neg_pos = Data_drug_volcano_pos
            } else {
              Data_drug_volcano_neg_pos = Data_drug_volcano_neg
            }
            df_volcano_tmp$Gene = gene_volcano
            df_volcano_tmp$Mut_respond_patients = sum(Data_drug_volcano_pos$Response)
            df_volcano_tmp$Mut_nonrespond_patients = nrow(Data_drug_volcano_pos) - df_volcano_tmp$Mut_respond_patients
            df_volcano_tmp$WT_respond_patients = sum(Data_drug_volcano_neg$Response)
            df_volcano_tmp$WT_nonrespond_patients = nrow(Data_drug_volcano_neg) - df_volcano_tmp$WT_respond_patients
            # print(gene_volcano)
            # print(drug_name)
            Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
            Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
            if(input$volcano_pathology == "Yes" &
               min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
               max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$mutation)>1 &
               sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$Response)>1 &
               sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
              Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos, family=binomial)
              df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
              df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
            } else {
              major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
              Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                dplyr::filter(Diagnosis %in% major_pathology)
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$Response)>1 &
                 sum(Data_drug_volcano_neg_pos_major$Response)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                Data_ORR = glm(Response ~ mutation + Diagnosis, data=Data_drug_volcano_neg_pos_major, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                        sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                        sum(Data_drug_volcano_neg_pos$Response)>1 &
                        sum(Data_drug_volcano_neg_pos$Response)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = glm(Response ~ mutation, data=Data_drug_volcano_neg_pos, family=binomial)
                df_volcano_tmp$OddsRatio_considering_pathology = exp(coef(Data_ORR))["mutation"]
                df_volcano_tmp$p_value_logistic_regression = summary(Data_ORR)$coefficients["mutation",4]
              } else {
                df_volcano_tmp$OddsRatio_considering_pathology = 1
                df_volcano_tmp$p_value_logistic_regression = 1
              }
            }
            df_volcano = rbind(df_volcano, df_volcano_tmp)
          }
          Figure_volcano = df_volcano %>%
            dplyr::filter(Regimen == paste(input$drug, collapse = ";")) %>%
            dplyr::mutate(Gene_effect = ifelse(p_value_logistic_regression >= 0.05 , "Not significant",
                                               ifelse(OddsRatio_considering_pathology >= 2, "Effective",
                                                      ifelse(OddsRatio_considering_pathology <= 0.5, "Not effective","Not significant"))))
          Figure_volcano$label = ifelse((Figure_volcano$OddsRatio_considering_pathology >=2 |
                                           Figure_volcano$p_value_logistic_regression < 0.05 |
                                           Figure_volcano$OddsRatio_considering_pathology <= 0.5), Figure_volcano$Gene, NA_character_)
          Figure_volcano$Gene_effect = factor(Figure_volcano$Gene_effect, levels = c("Not effective", "Not significant", "Effective"))
          Figure_volcano$OddsRatio_considering_pathology[Figure_volcano$OddsRatio_considering_pathology>16]=16
          Figure_volcano$OddsRatio_considering_pathology[Figure_volcano$OddsRatio_considering_pathology<1/16]=1/16
          
          Figure_volcano$log2OR = log2(Figure_volcano$OddsRatio_considering_pathology)
          Figure_volcano$log10pval = -log10(Figure_volcano$p_value_logistic_regression)
          if(nrow(Figure_volcano)>0){
            g_volcano[[figure_no]] =
              ggplot(data = Figure_volcano, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
              geom_point(size = 2) +
              labs(x = expression("log"[2]*"odds ratio of objective response"), y = expression("-log"[10]*"p-value")) +
              scale_color_manual(values = c("#bb0c00", "grey", "#00AFBB"), 
                                 labels = c("Not effective", "Not significant", "Effective"),
                                 drop = FALSE) +
              geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
              geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
              geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
              theme(legend.position = "none", panel.grid = element_blank()) +
              ggtitle(paste0(paste(input$drug, collapse = ";"), " ",
                             Figure_volcano$Total_treated_patients[[1]], " patients")) +
              geom_text_repel(max.overlaps = Inf) +
              theme_classic()
            figure_no = figure_no + 1
            drug_name_list = c(drug_name_list, paste(input$drug, collapse = ";"))
          }
          significant_genes = c(significant_genes, (Figure_volcano %>% dplyr::arrange(p_value_logistic_regression))$label)
          significant_genes = significant_genes[!is.na(significant_genes)]
        }
      }
    }
    
    if(length(drug_name_list) > 0){
      display_names <- paste0(drug_name_list, " (", seq_along(drug_name_list), ")")
      regimen_choice_RECIST = setNames(as.character(seq_along(drug_name_list)), display_names)
    } else {
      regimen_choice_RECIST = NULL
    }
    
    OUTPUT_DATA$drug_analysis_regimen_choice_RECIST = regimen_choice_RECIST
    OUTPUT_DATA$drug_analysis_g_volcano = g_volcano
    OUTPUT_DATA$drug_analysis_df_volcano = df_volcano
    OUTPUT_DATA$drug_analysis_significant_genes = significant_genes
  })
  rm(volcano_env)
  gc()
}


volcano_ToT_logic <- function() {
  volcano_env <- new.env()
  with(volcano_env, {
    Data_drug_TTF = OUTPUT_DATA$Data_drug_Data_drug_TTF
    Data_MAF_target_tmp = OUTPUT_DATA$Data_drug_Data_MAF_target_tmp
    significant_genes = NULL
    
    Data_drug_volcano_all_ToT = Data_drug_TTF
    drugs_volcano = sort(table(Data_drug_volcano_all_ToT$Drug),decreasing = T)
    drugs_volcano = drugs_volcano[!names(drugs_volcano) %in% c("Others", "Unknown", "")]
    drugs_volcano = drugs_volcano[!is.na(names(drugs_volcano))]
    drugs_volcano = drugs_volcano[1:min(200, length(drugs_volcano))]
    df_volcano_ToT = data.frame(Regimen = as.character(),
                                Gene = as.character(),
                                Total_treated_patients = as.numeric(),
                                Mut_patients = as.numeric(),
                                HazardRatio_considering_pathology = as.numeric(),
                                p_value = as.numeric()
    )
    df_volcano_tmp = data.frame(Regimen = "",
                                Gene = "",
                                Total_treated_patients = 0,
                                Mut_patients = 0,
                                HazardRatio_considering_pathology = 1,
                                p_value = 1
    )
    figure_no = 1
    g_volcano_ToT = list()
    options(warn = -1)
    drug_name_list = NULL
    withProgress(message = "volcane plot for treatment time", {
      for(drug_name in names(drugs_volcano)){
        Data_drug_volcano = Data_drug_volcano_all_ToT %>%
          dplyr::filter(Drug == drug_name) %>%
          dplyr::select(ID, Drug_length, Drug_length_censor, diagnosis)
        if(nrow(Data_drug_volcano) >= 1){
          Data_MAF_target_volcano = Data_MAF_target_tmp %>%
            dplyr::filter(Tumor_Sample_Barcode %in%
                            Data_drug_volcano$ID) %>%
            dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
          
          colnames(Data_drug_volcano) = c("ID", "Drug_length", "Drug_length_censor", "Diagnosis")
          colnames(Data_MAF_target_volcano) = c("ID", "Gene")
          Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
          
          df_volcano_tmp$Regimen = drug_name
          df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
          gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
          gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
          if(length(gene_table_volcano) > 0){
            for(gene_volcano in names(gene_table_volcano)){
              Data_drug_volcano_pos = Data_MAF_target_volcano %>%
                dplyr::filter(Gene == gene_volcano) %>%
                dplyr::select(-Gene)
              Data_drug_volcano_neg = Data_MAF_target_volcano %>%
                dplyr::select(-Gene) %>%
                dplyr::distinct()
              Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
              if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_pos$mutation = 1
              }
              if(nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg$mutation = 0
              }
              if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
                Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
              } else if(nrow(Data_drug_volcano_pos) > 0){
                Data_drug_volcano_neg_pos = Data_drug_volcano_pos
              } else{
                Data_drug_volcano_neg_pos = Data_drug_volcano_neg
              }
              df_volcano_tmp$Gene = gene_volcano
              df_volcano_tmp$Mut_patients = nrow(Data_drug_volcano_pos)
              Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
              Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                 sum(Data_drug_volcano_neg_pos$Drug_length_censor)>1 &
                 sum(Data_drug_volcano_neg_pos$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation+Diagnosis,
                                 data=Data_drug_volcano_neg_pos)
                data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
                df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
                df_volcano_tmp$p_value = data_tmp$p.value
              } else {
                major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
                Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                  dplyr::filter(Diagnosis %in% major_pathology)
                Data_drug_volcano_neg_pos_major$Diagnosis = factor(Data_drug_volcano_neg_pos_major$Diagnosis)
                Data_drug_volcano_neg_pos_major$Diagnosis = relevel(Data_drug_volcano_neg_pos_major$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
                if(input$volcano_pathology == "Yes" &
                   min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                   max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                   sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                   sum(Data_drug_volcano_neg_pos_major$Drug_length_censor)>1 &
                   sum(Data_drug_volcano_neg_pos_major$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                  Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation+Diagnosis,
                                   data=Data_drug_volcano_neg_pos)
                  data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
                  df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
                  df_volcano_tmp$p_value = data_tmp$p.value
                } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                          sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                          sum(Data_drug_volcano_neg_pos$Drug_length_censor)>1 &
                          sum(Data_drug_volcano_neg_pos$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos)-1)){
                  Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation,
                                   data=Data_drug_volcano_neg_pos)
                  data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
                  df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
                  df_volcano_tmp$p_value = data_tmp$p.value
                } else {
                  df_volcano_tmp$HazardRatio_considering_pathology = 1
                  df_volcano_tmp$p_value = 1
                }
              }
              df_volcano_ToT = rbind(df_volcano_ToT, df_volcano_tmp)
            }
            
            if(figure_no<200){
              Figure_volcano_tmp = df_volcano_ToT %>%
                dplyr::filter(Regimen == drug_name) %>%
                dplyr::mutate(Gene_effect = ifelse(p_value >= 0.05 , "Not significant",
                                                   ifelse(HazardRatio_considering_pathology > 1, "Not effective",
                                                          ifelse(HazardRatio_considering_pathology < 1, "Effective","Not significant"))))
              Figure_volcano_tmp$HazardRatio_considering_pathology[Figure_volcano_tmp$HazardRatio_considering_pathology>16]=16
              Figure_volcano_tmp$HazardRatio_considering_pathology[Figure_volcano_tmp$HazardRatio_considering_pathology<1/16]=1/16
              Figure_volcano_tmp$Gene_effect = factor(Figure_volcano_tmp$Gene_effect, levels = c("Effective", "Not significant", "Not effective"))
              Figure_volcano_tmp$label = ifelse((Figure_volcano_tmp$p_value < 0.05), Figure_volcano_tmp$Gene, NA_character_)
              
              Figure_volcano_tmp$log2OR = log2(Figure_volcano_tmp$HazardRatio_considering_pathology)
              Figure_volcano_tmp$log10pval = -log10(Figure_volcano_tmp$p_value)
              if(nrow(Figure_volcano_tmp)>0){
                g_volcano_ToT[[figure_no]] =
                  ggplot(data = Figure_volcano_tmp, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
                  geom_point(size = 2) +
                  labs(x = expression("log"[2]*"hazard ratio of treatment time"), y = expression("-log"[10]*"p-value")) +
                  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                                     labels = c("Effective", "Not significant", "Not effective"),
                                     drop = FALSE) +
                  geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
                  geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
                  theme(legend.position = "none", panel.grid = element_blank()) +
                  geom_text_repel(max.overlaps = Inf) +
                  ggtitle(paste0(drug_name, " ",
                                 df_volcano_tmp$Total_treated_patients, " patients")) +
                  theme_classic()
              }
              figure_no = figure_no + 1
              drug_name_list = c(drug_name_list, drug_name)
            }
          }
        }
        incProgress(1 / length(names(drugs_volcano)))
      }
    })
    options(warn = 0)
    if(!is.null(input$drug)){
      Data_drug_volcano = Data_drug_volcano_all_ToT %>%
        dplyr::filter(Drug %in% input$drug) %>%
        dplyr::select(ID, Drug_length, Drug_length_censor, diagnosis) %>%
        dplyr::distinct(ID, .keep_all = T)
      if(nrow(Data_drug_volcano) >= 1){
        Data_MAF_target_volcano = Data_MAF_target_tmp %>%
          dplyr::filter(Tumor_Sample_Barcode %in%
                          Data_drug_volcano$ID) %>%
          dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol)
        colnames(Data_drug_volcano) = c("ID", "Drug_length", "Drug_length_censor", "Diagnosis")
        colnames(Data_MAF_target_volcano) = c("ID", "Gene")
        Data_MAF_target_volcano = dplyr::left_join(Data_MAF_target_volcano, Data_drug_volcano, by="ID")
        df_volcano_tmp$Regimen = paste(input$drug, collapse = ";")
        df_volcano_tmp$Total_treated_patients = nrow(Data_drug_volcano)
        gene_table_volcano = sort(table(Data_MAF_target_volcano$Gene),decreasing = T)
        gene_table_volcano = gene_table_volcano[gene_table_volcano>=8]
        if(length(gene_table_volcano) > 0){
          for(gene_volcano in names(gene_table_volcano)){
            Data_drug_volcano_pos = Data_MAF_target_volcano %>%
              dplyr::filter(Gene == gene_volcano) %>%
              dplyr::select(-Gene)
            Data_drug_volcano_neg = Data_MAF_target_volcano %>%
              dplyr::select(-Gene) %>%
              dplyr::distinct()
            Data_drug_volcano_neg = setdiff(Data_drug_volcano_neg, Data_drug_volcano_pos)
            if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_pos$mutation = 1
            }
            if(nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg$mutation = 0
            }
            if(nrow(Data_drug_volcano_pos) > 0 & nrow(Data_drug_volcano_neg) > 0){
              Data_drug_volcano_neg_pos = rbind(Data_drug_volcano_pos, Data_drug_volcano_neg)
            } else if(nrow(Data_drug_volcano_pos) > 0){
              Data_drug_volcano_neg_pos = Data_drug_volcano_pos
            } else {
              Data_drug_volcano_neg_pos = Data_drug_volcano_neg
            }
            
            df_volcano_tmp$Gene = gene_volcano
            df_volcano_tmp$Mut_patients = nrow(Data_drug_volcano_pos)
            Data_drug_volcano_neg_pos$Diagnosis = factor(Data_drug_volcano_neg_pos$Diagnosis)
            Data_drug_volcano_neg_pos$Diagnosis = relevel(Data_drug_volcano_neg_pos$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
            if(input$volcano_pathology == "Yes" &
               min(table(Data_drug_volcano_neg_pos$Diagnosis))>1 &
               max(table(Data_drug_volcano_neg_pos$Diagnosis))<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$mutation)>1 &
               sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
               sum(Data_drug_volcano_neg_pos$Drug_length_censor)>1 &
               sum(Data_drug_volcano_neg_pos$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos)-1)){
              Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation+Diagnosis,
                               data=Data_drug_volcano_neg_pos)
              data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
              df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
              df_volcano_tmp$p_value = data_tmp$p.value
            } else {
              major_pathology = names(table(Data_drug_volcano_neg_pos$Diagnosis)[table(Data_drug_volcano_neg_pos$Diagnosis) > 1])
              Data_drug_volcano_neg_pos_major = Data_drug_volcano_neg_pos %>%
                dplyr::filter(Diagnosis %in% major_pathology)
              Data_drug_volcano_neg_pos_major$Diagnosis = factor(Data_drug_volcano_neg_pos_major$Diagnosis)
              Data_drug_volcano_neg_pos_major$Diagnosis = relevel(Data_drug_volcano_neg_pos_major$Diagnosis, ref=names(sort(table(Data_drug_volcano_neg_pos$Diagnosis),decreasing = T))[[1]])
              if(input$volcano_pathology == "Yes" &
                 min(table(Data_drug_volcano_neg_pos_major$Diagnosis))>1 &
                 max(table(Data_drug_volcano_neg_pos_major$Diagnosis))<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$mutation)>1 &
                 sum(Data_drug_volcano_neg_pos_major$mutation)<(nrow(Data_drug_volcano_neg_pos_major)-1) &
                 sum(Data_drug_volcano_neg_pos_major$Drug_length_censor)>1 &
                 sum(Data_drug_volcano_neg_pos_major$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos_major)-1)){
                Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation+Diagnosis,
                                 data=Data_drug_volcano_neg_pos)
                data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
                df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
                df_volcano_tmp$p_value = data_tmp$p.value
              } else if(sum(Data_drug_volcano_neg_pos$mutation)>1 &
                        sum(Data_drug_volcano_neg_pos$mutation)<(nrow(Data_drug_volcano_neg_pos)-1) &
                        sum(Data_drug_volcano_neg_pos$Drug_length_censor)>1 &
                        sum(Data_drug_volcano_neg_pos$Drug_length_censor)<(nrow(Data_drug_volcano_neg_pos)-1)){
                Data_ORR = coxph(Surv(Drug_length, Drug_length_censor)~mutation,
                                 data=Data_drug_volcano_neg_pos)
                data_tmp = tidy(Data_ORR, exponentiate=TRUE, conf.int=TRUE)[1,]
                df_volcano_tmp$HazardRatio_considering_pathology = data_tmp$estimate
                df_volcano_tmp$p_value = data_tmp$p.value
              } else {
                df_volcano_tmp$HazardRatio_considering_pathology = 1
                df_volcano_tmp$p_value = 1
              }
            }
            df_volcano_ToT = rbind(df_volcano_ToT, df_volcano_tmp)
          }
          Figure_volcano = df_volcano_ToT  %>%
            dplyr::filter(Regimen == paste(input$drug, collapse = ";")) %>%
            dplyr::mutate(Gene_effect = ifelse(p_value >= 0.05 , "Not significant",
                                               ifelse(HazardRatio_considering_pathology > 1, "Not effective",
                                                      ifelse(HazardRatio_considering_pathology < 1, "Effective","Not significant"))))
          Figure_volcano$HazardRatio_considering_pathology[Figure_volcano$HazardRatio_considering_pathology>16]=16
          Figure_volcano$HazardRatio_considering_pathology[Figure_volcano$HazardRatio_considering_pathology<1/16]=1/16
          Figure_volcano$Gene_effect = factor(Figure_volcano$Gene_effect, levels = c("Effective", "Not significant", "Not effective"))
          Figure_volcano$label = ifelse((Figure_volcano$p_value < 0.05), Figure_volcano$Gene, NA_character_)
          Figure_volcano$log2OR = log2(Figure_volcano$HazardRatio_considering_pathology)
          Figure_volcano$log10pval = -log10(Figure_volcano$p_value)
          if(nrow(Figure_volcano)>0){
            g_volcano_ToT[[figure_no]] =
              ggplot(data = Figure_volcano, aes(x = log2OR, y = log10pval, col = Gene_effect, label = label)) +
              geom_point(size = 2) +
              labs(x = expression("log"[2]*"hazard ratio of treatment time"), y = expression("-log"[10]*"p-value")) +
              scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                                 labels = c("Effective", "Not significant", "Not effective"),
                                 drop = FALSE) +
              geom_hline(yintercept = -log10(0.05), size = 0.2, color = "dark green") + 
              geom_vline(xintercept = -1, size = 0.2, color = "dark green", linetype = 'dashed') +
              geom_vline(xintercept = 1, size = 0.2, color = "dark green", linetype = 'dashed') +
              theme(legend.position = "none", panel.grid = element_blank()) +
              ggtitle(paste0(paste(input$drug, collapse = ";"), " ",
                             Figure_volcano$Total_treated_patients[[1]], " patients")) +
              geom_text_repel(max.overlaps = Inf) +
              theme_classic()
            figure_no = figure_no + 1
            drug_name_list = c(drug_name_list, paste(input$drug, collapse = ";"))
          }
          significant_genes = (Figure_volcano %>% dplyr::arrange(p_value))$label
          significant_genes = significant_genes[!is.na(significant_genes)]
        }
      }
    }
    if(length(drug_name_list) > 0){
      display_names <- paste0(drug_name_list, " (", seq_along(drug_name_list), ")")
      regimen_choice_ToT = setNames(as.character(seq_along(drug_name_list)), display_names)
    } else {
      regimen_choice_ToT = NULL
    }
    
    OUTPUT_DATA$drug_analysis_regimen_choice_ToT = regimen_choice_ToT
    OUTPUT_DATA$drug_analysis_g_volcano_ToT = g_volcano_ToT
    OUTPUT_DATA$drug_analysis_df_volcano_ToT = df_volcano_ToT
    OUTPUT_DATA$drug_analysis_significant_genes_ToT = significant_genes
  })
  rm(volcano_env)
  gc()
}