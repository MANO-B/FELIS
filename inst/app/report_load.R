if (CCAT_FLAG & file.exists(file.path(app_dir, "source", "variant_data_whole.qs"))) {
  initial_data_report <- QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(app_dir, "source", "variant_data_whole.qs")) %>%
    dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
  Data_report_raw_original <- reactive({ initial_data_report })
} else {
  Data_report_raw_original =  reactive({
    withProgress(message = "Mutation data loading.", {
      if(!is.null(input$new_analysis) && input$new_analysis =="No, use the previous dataset" &&
         file.exists(file.path(tempdir(), "variant_data.qs"))){
        Data_MAF = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=file.path(tempdir(), "variant_data.qs")) %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
      } else {
        req(input$report_files)
        Encode = guess_encoding(input$report_files[[1, 'datapath']])[[1]][[1]]
        if(Encode == "Shift_JIS"){
          Encode = "CP932"
        }
        clin_tmp <- read_csv(col_names = TRUE,
                             file = input$report_files$datapath,
                             locale = locale(encoding=Encode),
                             num_threads=max(1, parallel::detectCores() - 1, na.rm = TRUE),
                             progress=FALSE,
                             show_col_types=FALSE,
                             name_repair=make.names,
                             col_select = c(
                               C.CAT調査結果.基本項目.ハッシュID,
                               C.CAT調査結果.変異情報.変異種類,
                               C.CAT調査結果.変異情報.変異由来,
                               C.CAT調査結果.変異情報.マーカー,
                               C.CAT調査結果.変異情報.変異内容,
                               C.CAT調査結果.変異情報.マーカー詳細,
                               C.CAT調査結果.変異情報.変異アレル頻度,
                               C.CAT調査結果.変異情報.コピー数変化,
                               C.CAT調査結果.変異情報.クロモソーム番号,
                               C.CAT調査結果.変異情報.物理位置,
                               C.CAT調査結果.変異情報.リファレンス塩基,
                               C.CAT調査結果.変異情報.塩基変化,
                               C.CAT調査結果.変異情報.変異種類,
                               C.CAT調査結果.変異情報.変異アレル頻度,
                               C.CAT調査結果.変異情報.変異内容,
                               C.CAT調査結果.変異情報TMB.TMB,
                               C.CAT調査結果.エビデンス.エビデンスレベル,
                               C.CAT調査結果.エビデンス.薬剤
                             ),
                             col_types = cols(
                               C.CAT調査結果.基本項目.ハッシュID = col_character(),
                               C.CAT調査結果.変異情報.変異種類 = col_character(),
                               C.CAT調査結果.変異情報.変異由来 = col_character(),
                               C.CAT調査結果.変異情報.マーカー = col_character(),
                               C.CAT調査結果.変異情報.変異内容 = col_character(),
                               C.CAT調査結果.変異情報.マーカー詳細 = col_character(),
                               C.CAT調査結果.変異情報.変異アレル頻度 = col_character(),
                               C.CAT調査結果.変異情報.コピー数変化 = col_character(),
                               C.CAT調査結果.変異情報.クロモソーム番号 = col_character(),
                               C.CAT調査結果.変異情報.物理位置 = col_character(),
                               C.CAT調査結果.変異情報.リファレンス塩基 = col_character(),
                               C.CAT調査結果.変異情報.塩基変化 = col_character(),
                               C.CAT調査結果.変異情報.変異種類 = col_character(),
                               C.CAT調査結果.変異情報.変異アレル頻度 = col_character(),
                               C.CAT調査結果.変異情報.変異内容 = col_character(),
                               C.CAT調査結果.変異情報TMB.TMB = col_character(),
                               C.CAT調査結果.エビデンス.エビデンスレベル = col_character(),
                               C.CAT調査結果.エビデンス.薬剤 = col_character(),
                               .default = col_guess()
                             )) %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
        incProgress(1 / 4)
        #clin_tmp = as.data.frame(clin_tmp)
        clin_tmp = data.table(clin_tmp) %>%
          dplyr::filter(
            !str_detect(C.CAT調査結果.変異情報.マーカー, ",") &
              C.CAT調査結果.変異情報.マーカー != "" &
              C.CAT調査結果.変異情報.変異種類 != "expression" &
              !(C.CAT調査結果.エビデンス.エビデンスレベル == "C" &
                  C.CAT調査結果.変異情報.変異内容 == "high not detected" &
                  C.CAT調査結果.変異情報.マーカー == "MSI")
          ) %>% dplyr::mutate(
            C.CAT調査結果.変異情報.変異種類 = case_when(
              C.CAT調査結果.変異情報.変異種類 == "rearrangement" &
                C.CAT調査結果.変異情報.変異内容 == "rearrangements" ~ "rearrangement",
              str_detect(C.CAT調査結果.変異情報.変異種類, "exon skipping") ~ "exon_skipping",
              str_detect(C.CAT調査結果.変異情報.変異内容, "exon skipping") ~ "exon_skipping",
              C.CAT調査結果.変異情報.変異種類 == "rearrangement" &
                C.CAT調査結果.変異情報.変異内容 == "long deletion" ~ "deletion",
              C.CAT調査結果.変異情報.変異種類 == "rearrangement" ~ C.CAT調査結果.変異情報.変異内容,
              C.CAT調査結果.変異情報.変異種類 == "copy_number_alteration" ~ C.CAT調査結果.変異情報.変異内容,
              C.CAT調査結果.変異情報.変異種類 == "small_scale_variant" &
                str_detect(C.CAT調査結果.変異情報.変異内容, "\\*") ~ "nonsense_frameshift",
              C.CAT調査結果.変異情報.変異種類 == "small_scale_variant"  ~ "small_scale_variant",
              C.CAT調査結果.変異情報.変異種類 == "other_biomarker" ~ "TMB_MSI_high",
              TRUE ~ C.CAT調査結果.変異情報.変異種類
            ))
        QS_SAVE(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), clin_tmp, file=file.path(tempdir(), "variant_data.qs"))
      }
    })
  })
}
Data_report_raw =  reactive({
  req(Data_report_raw_original())
  req(input$fusion)
  req(input$T_N)
  withProgress(message = "Mutation data loading.", {
    clin_tmp = Data_report_raw_original()
    if(!is.null(input$reaanotation)){
      reaanotation_list = data.frame(NULL)
      for(i in 1:length(input$reaanotation[,1])){
        reaanotation_list <- rbind(reaanotation_list,
                                   read.csv(header = TRUE,
                                            file(input$reaanotation[[i, 'datapath']],
                                                 encoding='UTF-8-BOM')))
      }
      clin_tmp$Gene_alt = paste0(clin_tmp$C.CAT調査結果.変異情報.マーカー, "_", clin_tmp$C.CAT調査結果.変異情報.変異内容)
      clin_tmp$C.CAT調査結果.エビデンス.エビデンスレベル_tmp = unlist(lapply(list(clin_tmp$Gene_alt), function(x) {
        as.vector(reaanotation_list$Level[match(x, reaanotation_list$Gene_alt)])}))
      clin_tmp = clin_tmp %>% dplyr::mutate(
        C.CAT調査結果.エビデンス.エビデンスレベル = case_when(
          !is.na(C.CAT調査結果.エビデンス.エビデンスレベル_tmp) ~ C.CAT調査結果.エビデンス.エビデンスレベル_tmp,
          TRUE ~ C.CAT調査結果.エビデンス.エビデンスレベル
        ))
      clin_tmp = clin_tmp %>% dplyr::select(-Gene_alt, -C.CAT調査結果.エビデンス.エビデンスレベル_tmp)
    }
    clin_tmp = clin_tmp %>%
      dplyr::mutate(
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "NKX2\\-1", "NKX2XXXX1"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "NKX3\\-1", "NKX3XXXX1"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H1\\-2", "H1XXXX2"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3\\-3A", "H3XXXX3A"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3\\-3B", "H3XXXX3B"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3\\-4", "H3XXXX4"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3\\-5", "H3XXXX5"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "HLA\\-A", "HLAXXXXA"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "HLA\\-DRB1", "HLAXXXXDRB1")
      ) %>%
      dplyr::mutate(
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "-", "_")
      )

    if(input$fusion == "Split into two and treat as mutation in each gene"){
      if(any(str_detect(clin_tmp$C.CAT調査結果.変異情報.マーカー, "_"))){
        clin_tmp_base = clin_tmp %>%
          dplyr::filter(!str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion = clin_tmp %>%
          dplyr::filter(str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion_1 = clin_tmp_fusion %>%
          dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー,
                        C.CAT調査結果.変異情報.マーカー = str_split(C.CAT調査結果.変異情報.マーカー, "_",simplify = T)[,1])
        clin_tmp_fusion_2 = clin_tmp_fusion %>%
          dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー,
                        C.CAT調査結果.変異情報.マーカー = str_split(C.CAT調査結果.変異情報.マーカー, "_",simplify = T)[,2])
        clin_tmp = rbind(clin_tmp_base, clin_tmp_fusion_1, clin_tmp_fusion_2)
      }
    } else if(input$fusion == "Split into two and treat like EML-fusion, ALK-fusion"){
      if(any(str_detect(clin_tmp$C.CAT調査結果.変異情報.マーカー, "_"))){
        clin_tmp_base = clin_tmp %>%
          dplyr::filter(!str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion = clin_tmp %>%
          dplyr::filter(str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion_1 = clin_tmp_fusion %>%
          dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー,
                        C.CAT調査結果.変異情報.マーカー = paste0(str_split(C.CAT調査結果.変異情報.マーカー, "_",simplify = T)[,1], "_fusion"))
        clin_tmp_fusion_2 = clin_tmp_fusion %>%
          dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー,
                        C.CAT調査結果.変異情報.マーカー = paste0(str_split(C.CAT調査結果.変異情報.マーカー, "_",simplify = T)[,2], "_fusion"))
        clin_tmp = rbind(clin_tmp_base, clin_tmp_fusion_1, clin_tmp_fusion_2)
      }
    } else if(input$fusion == "All fusion genes are named as 'fusion_genes'"){
      if(any(str_detect(clin_tmp$C.CAT調査結果.変異情報.マーカー, "_"))){
        clin_tmp_base = clin_tmp %>%
          dplyr::filter(!str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion = clin_tmp %>%
          dplyr::filter(str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_tmp_fusion = clin_tmp_fusion %>%
          dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー,
                        C.CAT調査結果.変異情報.マーカー = "fusion_genes")
        clin_tmp = rbind(clin_tmp_base, clin_tmp_fusion)
      }
    }
    if(input$amplification == "Treat as indipendent gene"){
      clin_tmp_base = clin_tmp %>%
        dplyr::filter(C.CAT調査結果.変異情報.変異内容 != "amplification" &
                        C.CAT調査結果.変異情報.変異内容 != "loss")
      clin_tmp_amplification = clin_tmp %>%
        dplyr::filter(C.CAT調査結果.変異情報.変異内容 == "amplification") %>%
        dplyr::mutate(C.CAT調査結果.変異情報.マーカー = paste0(C.CAT調査結果.変異情報.マーカー, "_amp"))
      clin_tmp_loss = clin_tmp %>%
        dplyr::filter(C.CAT調査結果.変異情報.変異内容 == "loss") %>%
        dplyr::mutate(C.CAT調査結果.変異情報.マーカー = paste0(C.CAT調査結果.変異情報.マーカー, "_loss"))
      clin_tmp = rbind(clin_tmp_base, clin_tmp_amplification, clin_tmp_loss)
    }
    clin_tmp = clin_tmp %>%
      dplyr::mutate(
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "NKX2XXXX1", "NKX2_1"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "NKX3XXXX1", "NKX3_1"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H1XXXX2", "H1_2"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3XXXX4", "H3_4"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3XXXX5", "H3_5"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3XXXX3B", "H3_3B"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "HLAXXXXA", "HLA_A"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "HLAXXXXDRB1", "HLA_DRB1"),
        C.CAT調査結果.変異情報.マーカー = str_replace(C.CAT調査結果.変異情報.マーカー, "H3XXXX3A", "H3_3A")
      )
    incProgress(1 / 4)
    if(!is.null(input$mutation_rename)){
      mutation_rename_list = data.frame(NULL)
      for(i in 1:length(input$mutation_rename[,1])){
        mutation_rename_list <- rbind(mutation_rename_list,
                                      read.csv(header = TRUE,
                                               file(input$mutation_rename[[i, 'datapath']],
                                                    encoding='UTF-8-BOM')))
      }
      Gene_list = sort(unique(mutation_rename_list$Gene))
      clin = clin_tmp %>% dplyr::filter(!C.CAT調査結果.変異情報.マーカー %in% Gene_list)
      for(i in 1:length(Gene_list)){
        mutation_rename_list_tmp = mutation_rename_list[mutation_rename_list$Gene == Gene_list[i],]
        clin_tmp2 = clin_tmp %>% dplyr::filter(C.CAT調査結果.変異情報.マーカー == Gene_list[i])
        clin_tmp2$C.CAT調査結果.変異情報.変異内容_tmp = unlist(lapply(list(clin_tmp2$C.CAT調査結果.変異情報.変異内容), function(x) {
          as.vector(mutation_rename_list_tmp$Rename[match(x, mutation_rename_list_tmp$Mutation)])}))
        clin_tmp2 = clin_tmp2 %>% dplyr::mutate(
          C.CAT調査結果.変異情報.変異内容 = case_when(
            !is.na(C.CAT調査結果.変異情報.変異内容_tmp) ~ C.CAT調査結果.変異情報.変異内容_tmp,
            TRUE ~ "Other"
          )
        )
        clin_tmp2 = clin_tmp2 %>% dplyr::select(-C.CAT調査結果.変異情報.変異内容_tmp)
        clin = rbind(clin, clin_tmp2)
      }
      clin_tmp = clin
    }
    Data_MAF = clin_tmp
    Data_MAF$C.CAT調査結果.変異情報.物理位置2 =
      Data_MAF$C.CAT調査結果.変異情報.物理位置
    Data_MAF = Data_MAF %>%
      dplyr::select(C.CAT調査結果.変異情報.マーカー,
                    C.CAT調査結果.変異情報.クロモソーム番号,
                    C.CAT調査結果.変異情報.物理位置,
                    C.CAT調査結果.変異情報.物理位置2,
                    C.CAT調査結果.変異情報.リファレンス塩基,
                    C.CAT調査結果.変異情報.塩基変化,
                    C.CAT調査結果.変異情報.変異種類,
                    C.CAT調査結果.エビデンス.エビデンスレベル,
                    C.CAT調査結果.エビデンス.薬剤,
                    C.CAT調査結果.基本項目.ハッシュID,
                    C.CAT調査結果.変異情報.変異アレル頻度,
                    C.CAT調査結果.変異情報.変異内容,
                    C.CAT調査結果.変異情報TMB.TMB,
                    C.CAT調査結果.変異情報.変異由来)
    colnames(Data_MAF) = c("Hugo_Symbol",
                           "Chromosome",
                           "Start_Position",
                           "End_Position",
                           "Reference_Allele",
                           "Tumor_Seq_Allele2",
                           "Variant_Classification",
                           "Evidence_level",
                           "Drug",
                           "Tumor_Sample_Barcode",
                           "VAF",
                           "amino.acid.change",
                           "TMB",
                           "Somatic_germline")
    Data_MAF = Data_MAF %>% dplyr::mutate(
      Variant_Type = case_when(
        Hugo_Symbol == "MSI" ~ "MSI",
        Hugo_Symbol == "TMB" ~ "TMB",
        nchar(Reference_Allele) == 1 &
          nchar(Tumor_Seq_Allele2) == 1 ~ "SNP",
        nchar(Reference_Allele) == 2 &
          nchar(Tumor_Seq_Allele2) == 2 ~ "DNP",
        nchar(Reference_Allele) == 3 &
          nchar(Tumor_Seq_Allele2) == 3 ~ "TNP",
        nchar(Reference_Allele) == 1 &
          Reference_Allele != "-" &
          nchar(Tumor_Seq_Allele2) > 1 ~ "INS",
        nchar(Reference_Allele) > 1 &
          Tumor_Seq_Allele2 != "-" ~ Variant_Classification,
        TRUE ~ "Other"
      ))
    incProgress(1 / 4)
    TMB_tmp = str_split(Data_MAF$TMB, "Muts", simplify = TRUE)[,1]
    TMB_tmp = str_split(TMB_tmp, "mutations", simplify = TRUE)[,1]
    Data_MAF$TMB = as.numeric(TMB_tmp)
    if(nrow(Data_MAF[is.na(Data_MAF$TMB),]) > 0){
      Data_MAF[is.na(Data_MAF$TMB),]$TMB = 0
    }
    Data_report_TMB = Data_MAF %>%
      dplyr::filter(Hugo_Symbol == "TMB") %>%
      dplyr::distinct(Tumor_Sample_Barcode, TMB, .keep_all = T) %>%
      dplyr::arrange(Tumor_Sample_Barcode) %>%
      dplyr::mutate(
        Evidence_level = case_when(
          TMB < input$TMB_threshold ~ "",
          TMB >= input$TMB_threshold ~ "F",
          TRUE ~ ""
        ),
        amino.acid.change = case_when(
          TMB < input$TMB_threshold ~ "",
          TMB >= input$TMB_threshold ~ "high",
          TRUE ~ ""
        )
      )
    Data_MAF = Data_MAF %>%
      dplyr::arrange(Tumor_Sample_Barcode) %>%
      dplyr::filter(Hugo_Symbol != "TMB") %>%
      dplyr::mutate(Evidence_level = case_when(
        Hugo_Symbol == "MSI" & amino.acid.change == "high" & Evidence_level == "" ~ "F",
        TRUE ~ Evidence_level
      ))
    Data_MAF = rbind(Data_MAF, Data_report_TMB)
    if(nrow(Data_MAF[is.na(Data_MAF$Evidence_level),]) > 0){
      Data_MAF[is.na(Data_MAF$Evidence_level),]$Evidence_level = ""
    }
    if(!is.null(input$T_N)){
      if(input$T_N == "Only somatic for T/N panel"){
        Data_MAF = Data_MAF %>%
          dplyr::filter(Somatic_germline != "Germline")
      } else if(input$T_N == "Only germline"){
        Data_MAF = Data_MAF %>%
          dplyr::filter(Somatic_germline == "Germline")
      }
    }
    incProgress(1 / 4)
  })
  return(Data_MAF)
})

Data_report = reactive({
  req(input$special_gene_annotation)
  Data_MAF = Data_report_raw()
  if(!is.null(input$special_gene_annotation)){
    if(input$special_gene_annotation == "This gene is also analyzed for VUS and Benign" &
       !is.null(input$special_gene)){
      Data_MAF = Data_MAF %>% dplyr::mutate(
        Evidence_level = case_when(
          Hugo_Symbol == input$special_gene &
            Evidence_level == "" ~ "F",
          TRUE ~ Evidence_level
        )
      )
    }
  }
  if(!is.null(input$special_gene_independent)){
    if(input$special_gene_independent == "Treat as independent genes in analyses" &
       !is.null(input$special_gene) &
       !is.null(input$special_gene_mutation_1)){
      Data_MAF = Data_MAF %>% dplyr::mutate(
        Hugo_Symbol = case_when(
          Hugo_Symbol == input$special_gene &
            amino.acid.change %in% input$special_gene_mutation_1 ~ paste0(Hugo_Symbol, "_", input$special_gene_mutation_1_name),
          Hugo_Symbol == input$special_gene &
            amino.acid.change %in% input$special_gene_mutation_2 ~ paste0(Hugo_Symbol, "_", input$special_gene_mutation_2_name),
          Hugo_Symbol == input$special_gene ~ paste0(Hugo_Symbol, "_NOS"),
          TRUE ~ Hugo_Symbol
        )
      )
    }
  }
  if(!is.null(input$gene_group_analysis)){
    if(input$gene_group_analysis == "Only cases with mutations in the gene set are analyzed" &
       !is.null(input$gene)){
      if(input$patho == "Only pathogenic muts"){
        ID_geneset = unique((Data_MAF %>% dplyr::filter(
          Hugo_Symbol %in% input$gene & Evidence_level == "F"
        ))$Tumor_Sample_Barcode)
      } else if(input$patho == "All muts"){
        ID_geneset = unique((Data_MAF %>% dplyr::filter(
          Hugo_Symbol %in% input$gene
        ))$Tumor_Sample_Barcode)
      }
      Data_MAF = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        ID_geneset)
    }
    if(input$gene_group_analysis == "Only cases without mutations in the gene set are analyzed" &
       !is.null(input$gene)){
      if(input$patho == "Only pathogenic muts"){
        ID_geneset = unique((Data_MAF %>% dplyr::filter(
          Hugo_Symbol %in% input$gene & Evidence_level == "F"
        ))$Tumor_Sample_Barcode)
      } else if(input$patho == "All muts"){
        ID_geneset = unique((Data_MAF %>% dplyr::filter(
          Hugo_Symbol %in% input$gene
        ))$Tumor_Sample_Barcode)
      }
      Data_MAF = Data_MAF %>%
        dplyr::filter(!Tumor_Sample_Barcode %in% ID_geneset)
      ALL_ID =  (Data_case_raw() %>% dplyr::filter(
        !C.CAT調査結果.基本項目.ハッシュID %in% ID_geneset
      ))$C.CAT調査結果.基本項目.ハッシュID
      DEFF_ID = ALL_ID[!ALL_ID %in% Data_MAF$Tumor_Sample_Barcode]
      if(length(DEFF_ID)>0){
        TEST_ALL = NULL
        TEST = (Data_MAF %>% dplyr::filter(
          Hugo_Symbol == "TMB"
        ))[1,]
        TEST$TMB=0
        TEST$Evidence_level = ""
        TEST$Drug = ""
        for(i in 1:length(DEFF_ID)){
          TEST$Tumor_Sample_Barcode = DEFF_ID[i]
          TEST_ALL = rbind(TEST_ALL, TEST)
        }
        Data_MAF = rbind(Data_MAF, TEST_ALL)
      }
    }
  }
  Data_MAF = Data_MAF %>%
    dplyr::distinct() %>%
    dplyr::arrange(Tumor_Sample_Barcode) %>%
    dplyr::select(-Somatic_germline)
  return(Data_MAF)
})
