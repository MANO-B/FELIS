# Mutually exclusive or co-occurrence
mutually_exclusive_logic <- function() {
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
                      症例.基本情報.性別.名称.,
                      症例.基本情報.がん種.OncoTree.,
                      症例.基本情報.がん種.OncoTree..名称.,
                      症例.基本情報.がん種.OncoTree.LEVEL1.,
                      症例.検体情報.パネル.名称.,
                      症例.背景情報.ECOG.PS.名称.,
                      症例.背景情報.喫煙歴有無.名称.,
                      症例.背景情報.アルコール多飲有無.名称.,
                      症例.背景情報.重複がん有無.異なる臓器..名称.,
                      症例.背景情報.多発がん有無.同一臓器..名称.
        )
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
      if(nrow(Data_MAF[Data_MAF$TMB > 30,]) > 0){
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
      Data_MAF_target_dt <- as.data.table(Data_MAF_target)
      Data_MAF_target_dt <- unique(Data_MAF_target_dt, by = c("Tumor_Sample_Barcode", "Hugo_Symbol"))
      Gene_names <- unique(c(input$gene, names(sort(table(Data_MAF_target_dt$Hugo_Symbol), decreasing = TRUE))))
      Gene_names <- Gene_names[Gene_names %in% unique(Data_MAF_target_dt$Hugo_Symbol)]
      Gene_names <- Gene_names[1:min(length(unique(Data_MAF_target_dt$Hugo_Symbol)), input$gene_no)]
      Patient_names <- unique(Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
      Mutation_matrix_dt <- dcast(Data_MAF_target_dt[Hugo_Symbol %in% Gene_names],
                                  Hugo_Symbol ~ Tumor_Sample_Barcode,
                                  fun.aggregate = length,
                                  value.var = "Tumor_Sample_Barcode")
      setDF(Mutation_matrix_dt)
      rownames(Mutation_matrix_dt) <- Mutation_matrix_dt$Hugo_Symbol
      Mutation_matrix_dt$Hugo_Symbol <- NULL
      Mutation_matrix <- as.matrix(Mutation_matrix_dt)
      #Mutation_matrix <- Mutation_matrix[Gene_names, , drop = FALSE]
      # NAを除外
      valid_genes <- Gene_names[!is.na(Gene_names)]
      
      # 行数・列数が0であっても安全にサブセット
      if (length(valid_genes) > 0 && nrow(Mutation_matrix) > 0 && ncol(Mutation_matrix) > 0) {
        Mutation_matrix <- Mutation_matrix[valid_genes, , drop = FALSE]
      } else {
        # 空の行列を明示的に定義（必要に応じて調整）
        Mutation_matrix <- matrix(nrow = 0, ncol = ncol(Mutation_matrix),
                                  dimnames = list(NULL, colnames(Mutation_matrix)))
      }
      all_patients <- unique(Data_case_target$C.CAT調査結果.基本項目.ハッシュID)
      missing_patients <- setdiff(all_patients, colnames(Mutation_matrix))
      if(length(missing_patients) > 0){
        additional_cols <- matrix(0, nrow = nrow(Mutation_matrix), ncol = length(missing_patients))
        colnames(additional_cols) <- missing_patients
        rownames(additional_cols) <- rownames(Mutation_matrix)
        Mutation_matrix <- cbind(Mutation_matrix, additional_cols)
      }
      Mutation_matrix <- Mutation_matrix[, all_patients, drop = FALSE]
      if(dim(Mutation_matrix)[1]>0){
        PMA <- getPM(Mutation_matrix)
        mutually_exclusive <- getMutex(A=Mutation_matrix,
                                       method = "ShiftedBinomial",
                                       #method = "Exact",
                                       #method = "Binomial",
                                       #method = "RefinedNormal",
                                       PM=PMA,
                                       lower.tail=TRUE
        )
        colnames(mutually_exclusive) = Gene_names
        rownames(mutually_exclusive) = Gene_names
        mutually_exclusive[lower.tri(mutually_exclusive,diag=TRUE)] <- NA
        Data_mut_ex = expand.grid(Gene_names, Gene_names)
        colnames(Data_mut_ex) = c("Gene_1", "Gene_2")
        Data_mut_ex$P_value = mutually_exclusive@x
        Data_mut_ex = Data_mut_ex %>%
          dplyr::mutate(P_value = case_when(
            P_value > 1-10^-10 ~ 1-10^-10,
            P_value < 10^-10 ~ 10^-10,
            TRUE ~ P_value
          ))
        Data_mut_ex$Gene_1 <- factor(Data_mut_ex$Gene_1, levels = Gene_names)
        Data_mut_ex$Gene_2 <- factor(Data_mut_ex$Gene_2, levels = Gene_names)
        Data_mut_ex_table = Data_mut_ex %>%
          dplyr::select(Gene_1, Gene_2, P_value)
        Data_mut_ex_table$Gene1mut_Gene2mut = 0
        Data_mut_ex_table$Gene1mut_Gene2wt = 0
        Data_mut_ex_table$Gene1wt_Gene2mut = 0
        Data_mut_ex_table$Gene1wt_Gene2wt = 0
        Data_mut_ex_table$OddsRatio = 0
        Data_mut_ex_table$OddsRatio_999permil_CI_LL = 0
        Data_mut_ex_table$OddsRatio_999permil_CI_UL = 0
        Gene_num = length(Gene_names)
        Single_mut = rowSums(Mutation_matrix)
        Total_patients = ncol(Mutation_matrix)
        for(i in 1:Gene_num){
          for(j in 1:Gene_num){
            tmp = sum(Mutation_matrix[i,] & Mutation_matrix[j,])
            Data_mut_ex_table$Gene1mut_Gene2mut[(i-1)*Gene_num+j] = tmp
            Data_mut_ex_table$Gene1mut_Gene2wt[(i-1)*Gene_num+j] = Single_mut[j] - tmp
            Data_mut_ex_table$Gene1wt_Gene2mut[(i-1)*Gene_num+j] = Single_mut[i] - tmp
            Data_mut_ex_table$Gene1wt_Gene2wt[(i-1)*Gene_num+j] = Total_patients - Single_mut[i] - Single_mut[j] + tmp
            tmp = fisher.test(matrix(c(tmp, Single_mut[i] - tmp, Single_mut[j] - tmp, Total_patients - Single_mut[i] - Single_mut[j] + tmp), nrow = 2), conf.level = 0.999)
            Data_mut_ex_table$OddsRatio[(i-1)*Gene_num+j] = tmp$estimate[[1]]
            Data_mut_ex_table$OddsRatio_999permil_CI_LL[(i-1)*Gene_num+j] = tmp$conf.int[[1]]
            Data_mut_ex_table$OddsRatio_999permil_CI_UL[(i-1)*Gene_num+j] = tmp$conf.int[[2]]
          }
        }
        Data_mut_ex_table = Data_mut_ex_table %>%
          dplyr::filter(!is.na(P_value))
        Data_mut_ex$OddsRatio = NA
        no_tmp = ((length(Data_mut_ex[,1]))^0.5)
        for(i in 1:length(Data_mut_ex[,1])){
          if(i %% no_tmp != 0 & i %% no_tmp <= floor(i/no_tmp)){
            Data_mut_ex[i,"OddsRatio"] = (Data_mut_ex_table %>% dplyr::filter(Gene_1 == Data_mut_ex[i,"Gene_1"] & Gene_2 == Data_mut_ex[i,"Gene_2"]))$OddsRatio
          }
        }
        OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex_table = format_numeric_columns(Data_mut_ex_table)
        
        Data_mut_ex$log10_OddsRatio = -log10(Data_mut_ex$OddsRatio)
        Data_mut_ex = Data_mut_ex %>%
          dplyr::mutate(OddsRatio_figure = case_when(
            log10_OddsRatio > 3 ~ 3,
            log10_OddsRatio < -3 ~ -3,
            TRUE ~ log10_OddsRatio
          ),
          comut_exclusive_OR = case_when(
            log10_OddsRatio <= 0 ~ "Co-occurrence based on odds ratio",
            TRUE ~ "Mutually exclusive based on odds ratio"
          ),
          comut_exclusive_OR_NB = case_when(
            P_value >= 0.5 ~ "Co-occurrence based on negative binomial distribution",
            TRUE ~ "Mutually exclusive based on odds ratio based on negative binomial distribution"
          ),
          P_value_OR = case_when(
            log10_OddsRatio <= 0 ~ 1 - P_value,
            TRUE ~ P_value
          ),
          P_value_NB = case_when(
            P_value >= 0.5 ~ 1 - P_value,
            TRUE ~ P_value
          ),
          P_co_or_ex_OR = case_when(
            log10_OddsRatio >= 0 ~ log(P_value_OR,base=10),
            TRUE ~ -log(P_value_OR,base=10)
          ),
          P_figure_OR = case_when(
            P_co_or_ex_OR > 3 ~ 3,
            P_co_or_ex_OR < -3 ~ -3,
            TRUE ~ P_co_or_ex_OR
          ),
          P_signif_OR = case_when(
            P_co_or_ex_OR > 3 ~ "#",
            P_co_or_ex_OR < -3 ~ "*",
            TRUE ~ ""
          ),
          P_co_or_ex_NB = case_when(
            P_value >= 0.5 ~ log(1-P_value,base=10),
            TRUE ~ -log(P_value,base=10)
          ),
          P_figure_NB = case_when(
            P_co_or_ex_NB > 3 ~ 3,
            P_co_or_ex_NB < -3 ~ -3,
            TRUE ~ P_co_or_ex_NB
          ),
          P_signif_NB = case_when(
            P_co_or_ex_NB <= -3 ~ "#",
            P_co_or_ex_NB >= 3 ~ "*",
            TRUE ~ ""
          )
        )
        OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex = Data_mut_ex
      }
      incProgress(1 / 1)
    })
  })
  rm(analysis_env)
  gc()
} 

output$table_mutually_exclusive = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex_table)
  create_datatable_with_confirm(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex_table, "FELIS downloaded forest plot data in Mutual exclusivity tab")
})


output$figure_mutually_exclusive_1 = renderGirafe({
  req(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex)
  ggiraph::girafe(ggobj = OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex %>%
                    ggplot(aes(x = Gene_1, y = Gene_2, fill = P_figure_NB)) +
                    geom_tile_interactive(color = "white") +
                    scale_fill_gradient2_interactive(low = "red", high = "blue", mid = "white",
                                                     midpoint = 0, limit = c(-3, 3),  breaks = -3:3,
                                                     labels = c("<-3", "-2", "-1", "0", "1", "2", ">3"),
                                                     name = "-log10(P-value)", space = "Lab",
                                                     na.value = "white") +
                    geom_text_interactive(aes(label = P_signif_NB),
                                          color = "black", size = 3)+
                    scale_y_discrete(limits = rev(levels(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex$Gene_2))) +
                    theme_classic() +
                    ggtitle("P value: Mutually exclusive (blue) or co-occurrence (red)\n*: p<0.001 for mutually exclusive, #: p<0.001 for co-occurrence") +
                    theme(legend.key.width = unit(0.4, "cm"),  # 凡例の色バーの幅を狭く
                          plot.title = element_text(size = 8, face = "bold"),         # タイトルのフォントサイズ
                          axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size=10),
                          axis.text.y = element_text(face = "italic", size=10),
                          axis.title.x = element_text(size=10),
                          axis.title.y = element_text(size=10),
                          legend.title = element_text(size = 8),                      # 凡例タイトルのフォントサイズ
                          legend.text = element_text(size = 8)),
                  width_svg = 10, height_svg = 10, options = list(opts_sizing(rescale = FALSE)))})
output$figure_mutually_exclusive_2 = renderGirafe({
  req(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex)
  ggiraph::girafe(ggobj = OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex %>% 
                    ggplot(aes(x = Gene_1, y = Gene_2, fill = OddsRatio_figure)) + 
                    geom_tile_interactive(color = "white") + 
                    scale_fill_gradient2_interactive(low = "red", high = "blue", mid = "white", 
                                                     midpoint = 0, limit = c(-3, 3), breaks = -3:3,
                                                     labels = c("<-3", "-2", "-1", "0", "1", "2", ">3"),
                                                     name = "-log10(OddsRatio)", space = "Lab",
                                                     na.value = "white")+
                    geom_text_interactive(aes(label = P_signif_OR),
                                          color = "black", size = 3)+
                    scale_y_discrete(limits = rev(levels(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex$Gene_2)))+
                    theme_classic() +
                    ggtitle("Odds ratio: Mutually exclusive (blue) or co-occurrence (red)\n*: p<0.001 for mutually exclusive, #: p<0.001 for co-occurrence") +
                    theme(legend.key.width = unit(0.4, "cm"),  # 凡例の色バーの幅を狭く
                          plot.title = element_text(size = 8, face = "bold"),         # タイトルのフォントサイズ
                          axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size=10),
                          axis.text.y = element_text(face = "italic", size=10),
                          axis.title.x = element_text(size=10),
                          axis.title.y = element_text(size=10),
                          legend.title = element_text(size = 8),                      # 凡例タイトルのフォントサイズ
                          legend.text = element_text(size = 8)),
                  width_svg = 10, height_svg = 10, options = list(opts_sizing(rescale = FALSE)))})
