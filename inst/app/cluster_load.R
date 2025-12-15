Data_cluster_ID = reactive({
  if(is.null(input$patho)){
    return(NULL)
  } else {
    withProgress(message = "Clustering based on gene mutation pattern.", {
      if(input$clustering == "Fixed to the whole C-CAT cohort" &
         file.exists("source/Data_mutation_cord_whole.qs")){
        Data_mutation_cord = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE),file="source/Data_mutation_cord_whole.qs") %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
      } else if(input$clustering == "Fixed to the 1st analysis" &
                file.exists(file.path(tempdir(), "Data_mutation_cord.qs"))){
        Data_mutation_cord = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE),file=file.path(tempdir(), "Data_mutation_cord.qs")) %>%
          dplyr::filter(!C.CAT調査結果.基本項目.ハッシュID %in% ID_exclude)
      } else {
        Data_case_target = Data_case() %>%
          dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = T)
        if(input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed", "Only cases without mutations in the gene set are analyzed")){
          Data_case_target = Data_case_target %>%
            dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                            Data_report()$Tumor_Sample_Barcode)
        }
        Data_MAF_target = Data_report() %>%
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
        Data_MAF_target = Data_MAF_target %>%
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
        Gene_list = unique(names(sort(table(Data_MAF_target$Hugo_Symbol),
                                      decreasing = T)))
        Data_mutation = data.frame(matrix(0,
                                          nrow=nrow(Data_case_target),
                                          ncol=length(Gene_list)))
        colnames(Data_mutation) = Gene_list
        rownames(Data_mutation) = Data_case_target$C.CAT調査結果.基本項目.ハッシュID
        for(i in 1:nrow(Data_MAF_target)){
          Data_mutation[Data_MAF_target$Tumor_Sample_Barcode[i],
                        Data_MAF_target$Hugo_Symbol[i]] = 1
        }
        if(length(Data_mutation[is.na(Data_mutation)])> 0){
          Data_mutation[is.na(Data_mutation)] <- 0
        }
        n_neighbors_ = min(15, dim(Data_mutation)[1])
        set.seed(1)
        Data_mutation_umap = umap::umap(Data_mutation, random_state=1, n_neighbors = n_neighbors_)
        Data_mutation_cord = Data_mutation_umap$layout %>% as.data.frame()
        EPS = input$eps
        MINPTS = 3
        set.seed(1)
        dbscan_res_changed <- dbscan::dbscan(Data_mutation_cord, eps = EPS, minPts = MINPTS)
        if(min(dbscan_res_changed$cluster) == 0){
          dbscan_res_changed$cluster = (dbscan_res_changed$cluster + 1)
        }
        Data_mutation_cord$C.CAT調査結果.基本項目.ハッシュID = rownames(Data_mutation)
        Data_mutation_cord$driver_mutations = rowSums(Data_mutation)
        Data_mutation_cord$cluster = dbscan_res_changed$cluster
        if(input$clustering == "Fixed to the 1st analysis"){
          QS_SAVE(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), Data_mutation_cord, file=file.path(tempdir(), "Data_mutation_cord.qs"))
        }
      }
      return(Data_mutation_cord)
    })
  }
})
