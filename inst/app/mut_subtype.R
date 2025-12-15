mut_subtype_logic <- function() {
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
      if(length(Data_MAF[Data_MAF$TMB > 30,]$TMB) > 0){
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
      Data_MAF_target = Data_MAF_target  %>%
        dplyr::distinct(Tumor_Sample_Barcode,
                        Hugo_Symbol,
                        .keep_all = TRUE)
      gene_to_analyze = unique(c(input$gene,
                                 names(sort(table(Data_MAF_target$Hugo_Symbol),
                                            decreasing = T))))
      gene_to_analyze = gene_to_analyze[gene_to_analyze %in% unique(Data_MAF_target$Hugo_Symbol)]
      gene_to_analyze = gene_to_analyze[1:min(length(unique(Data_MAF_target$Hugo_Symbol)),input$gene_no)]
      Diseases = sort(unique(Data_case_target$症例.基本情報.がん種.OncoTree.))
      Summary_Gene_alteration = data.frame(matrix(0,
                                                  nrow=length(Diseases),
                                                  ncol=length(gene_to_analyze)))
      colnames(Summary_Gene_alteration) = gene_to_analyze
      rownames(Summary_Gene_alteration) = Diseases
      
      for(i in Diseases){
        Data_tmp = Data_case_target %>%
          dplyr::filter(症例.基本情報.がん種.OncoTree. == i)
        Data_tmp = unique(Data_tmp$C.CAT調査結果.基本項目.ハッシュID)
        patient_no = length(Data_tmp)
        for(j in gene_to_analyze){
          Summary_Gene_alteration[i,j] =
            length(unique((Data_MAF_target %>% dplyr::filter(
              Hugo_Symbol == j &
                Tumor_Sample_Barcode %in% Data_tmp))$Tumor_Sample_Barcode)) / patient_no * 100
        }
      }
      
      Data_tmp = Summary_Gene_alteration
      Data_tmp$Disease = rownames(Data_tmp)
      if(dim(Data_tmp)[1]>0){
        Data_tmp = Data_tmp %>% gather(value = "freq", key = Gene, -Disease)
        Data_tmp <- transform(Data_tmp, Disease= factor(Disease, levels = sort(unique(Disease), decreasing = TRUE)))
        row_height <- 20
        column_width <- 20
        n_diseases <- length(unique(Data_tmp$Disease))
        n_genes <- length(unique(Data_tmp$Gene))
        OUTPUT_DATA$figure_mut_subtype_plotHeight = max(min(5000, (row_height * (n_diseases + 4))), 100)
        OUTPUT_DATA$figure_mut_subtype_plotWidth = max(min(5000, (column_width * (n_genes + 15))), 100)
        OUTPUT_DATA$figure_mut_subtype_plot_Data_tmp = Data_tmp
      }
    })
  })
  rm(analysis_env)
  gc()
}

output$figure_mut_subtype_plot <- renderGirafe({
  req(OUTPUT_DATA$figure_mut_subtype_plot_Data_tmp)
  ggiraph::girafe(ggobj = ggplot(OUTPUT_DATA$figure_mut_subtype_plot_Data_tmp,
                                 aes(as.factor(Gene), as.factor(Disease))) +
                    geom_tile_interactive(aes(fill = freq), color = "black",
                                          linetype = 1) + 
                    geom_text_interactive(aes(label = as.integer(freq)), size=3.2) +
                    scale_fill_gradient_interactive(low = "white", high = "red", limit=c(0,100), name="Frequency (%)") +
                    theme_classic() +
                    labs(title = "Frequently mutated genes and diagnoses") +
                    theme(legend.key.width = unit(0.4, "cm"),  # 凡例の色バーの幅を狭く
                          plot.title = element_text(size = 8, face = "bold"),         # タイトルのフォントサイズ
                          legend.title = element_text(size = 8),                      # 凡例タイトルのフォントサイズ
                          legend.text = element_text(size = 8),                       
                          axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size=10),
                          axis.text.y = element_text(size=10),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.line.x = element_blank(),
                          axis.line.y = element_blank()), options = list(opts_sizing(rescale = FALSE)),
                  width_svg = OUTPUT_DATA$figure_mut_subtype_plotWidth / 70,
                  height_svg = OUTPUT_DATA$figure_mut_subtype_plotHeight / 70)
})
