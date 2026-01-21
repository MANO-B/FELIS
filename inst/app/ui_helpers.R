output$select_evidence_level <- renderUI({
  checkboxGroupInput(
    inputId = "evidence_level",
    label   = "Select evidence levels for CGP benefit analysis",
    choices = c("A", "B", "C", "D", "E"),
    selected = c("A", "B", "C")
  )
})
output$select_minimum_courses = renderUI({
  numericInput("minimum_courses", "Minimum number of courses required for an independent regimen (numbers below this are classified as Others)",
               value = 10,
               min = 0,
               max = 1000,
               step = 10)
})
output$select_drug_table_layout = renderUI({
  radioButtons(
    "drug_table_layout", "Display percentage",
    choices = c("Yes", "No"),
    selected = "No")
})

output$select_figure_CI_AE_var = renderUI({
  req(OUTPUT_DATA$drug_table_Drug_summary)
  pickerInput("figure_CI_AE_var", "Choose drugs for adverse effect analysis",
              choices = sort(unique(OUTPUT_DATA$drug_table_Drug_summary$Drug)),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = names(sort(table(OUTPUT_DATA$drug_table_Drug_summary$Drug), decreasing = TRUE))[[1]], multiple = TRUE)
})

output$select_drug = renderUI({
  req(OUTPUT_DATA$drug_table_Drug_summary)
  pickerInput("drug", "Choose drugs for treatment effect analysis",
              choices = sort(unique(OUTPUT_DATA$drug_table_Drug_summary$Drug)),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = names(sort(table(OUTPUT_DATA$drug_table_Drug_summary$Drug), decreasing = TRUE))[[1]], multiple = TRUE)
})
output$select_drug_group_1 = renderUI({
  req(input$drug)
  pickerInput("drug_group_1", "Drug set 1 to be combined (if any)",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$drug, multiple = TRUE)
})
output$select_drug_group_1_name = renderUI({
  req(input$drug_group_1)
  textInput(inputId = "drug_group_1_name", label = "Name for Drug set 1",
            value = input$drug_group_1[[1]])
})
output$select_drug_group_2 = renderUI({
  req(input$drug_group_1)
  pickerInput("drug_group_2", "Drug set 2 to be combined (if any)",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$drug[-which(input$drug %in% input$drug_group_1)], multiple = TRUE)
})
output$select_drug_group_2_name = renderUI({
  req(input$drug_group_2)
  textInput(inputId = "drug_group_2_name", label = "Name for Drug set 2",
            value = input$drug_group_2[[1]])
})
output$select_drug_group_3 = renderUI({
  req(input$drug_group_2)
  pickerInput("drug_group_3", "Drug set 3 to be combined (if any)",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$drug[-which(input$drug %in%
                                  c(input$drug_group_1, input$drug_group_2))], multiple = TRUE)
})
output$select_drug_group_3_name = renderUI({
  req(input$drug_group_3)
  textInput(inputId = "drug_group_3_name", label = "Name for Drug set 3",
            value = input$drug_group_3[[1]])
})
output$select_drug_group_4 = renderUI({
  req(input$drug_group_3)
  pickerInput("drug_group_4", "Drug set 4 to be combined (if any)",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$drug[-which(input$drug %in%
                                  c(input$drug_group_1, input$drug_group_2, input$drug_group_3))], multiple = TRUE)
})
output$select_drug_group_4_name = renderUI({
  req(input$drug_group_4)
  textInput(inputId = "drug_group_4_name", label = "Name for Drug set 4",
            value = input$drug_group_4[[1]])
})

output$select_organ = renderUI({
  pickerInput("organ", "Filter by organ",
              choices = sort(unique(Data_case_raw()$症例.基本情報.がん種.OncoTree.LEVEL1.)),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = sort(unique(Data_case_raw()$症例.基本情報.がん種.OncoTree.LEVEL1.)), multiple = TRUE)
})
output$select_histology = renderUI({
  pickerInput("histology", "Filter by histology",
              choices = sort(unique((Data_case_raw() %>% filter(症例.基本情報.がん種.OncoTree.LEVEL1. %in% input$organ))$症例.基本情報.がん種.OncoTree..名称.)),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = sort(unique((Data_case_raw() %>% filter(症例.基本情報.がん種.OncoTree.LEVEL1. %in% input$organ))$症例.基本情報.がん種.OncoTree..名称.)), multiple = TRUE)
})
output$select_histology_group_1 = renderUI({
  req(input$histology)
  pickerInput("histology_group_1", "Histology type 1 to be combined (if none, not selected)",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$histology, multiple = TRUE)
})

output$select_histology_group_1_name = renderUI({
  req(input$histology_group_1)
  pickerInput("histology_group_1_name", "Name for histology type 1",
              options = list(`actions-box` = TRUE),
              input$histology_group_1, multiple = FALSE)
})
output$select_histology_group_2 = renderUI({
  req(input$histology_group_1)
  pickerInput("histology_group_2", "Histology type 2 to be combined",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$histology[-which(input$histology %in%
                                       input$histology_group_1)], multiple = TRUE)
})
output$select_histology_group_2_name = renderUI({
  req(input$histology_group_2)
  pickerInput("histology_group_2_name", "Name for histology type 2",
              input$histology_group_2, multiple = FALSE)
})

output$select_histology_group_3 = renderUI({
  req(input$histology_group_2)
  pickerInput("histology_group_3", "Histology type 3 to be combined",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$histology[-which(input$histology %in%
                                       c(input$histology_group_1,
                                         input$histology_group_2))], multiple = TRUE)
})

output$select_histology_group_3_name = renderUI({
  req(input$histology_group_3)
  pickerInput("histology_group_3_name", "Name for histology type 3",
              options = list(`actions-box` = TRUE),
              input$histology_group_3, multiple = FALSE)
})
output$select_histology_group_4 = renderUI({
  req(input$histology_group_3)
  pickerInput("histology_group_4", "Histology type 4 to be combined",
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              input$histology[-which(input$histology %in%
                                       c(input$histology_group_1,
                                         input$histology_group_2,
                                         input$histology_group_3))], multiple = TRUE)
})
output$select_histology_group_4_name = renderUI({
  req(input$histology_group_4)
  pickerInput("histology_group_4_name", "Name for histology type 4",
              options = list(`actions-box` = TRUE),
              input$histology_group_4, multiple = FALSE)
})
output$select_sex = renderUI({
  pickerInput("sex", "Filter by sex",
              choices = sort(unique(Data_case_raw()$症例.基本情報.性別.名称.)),
              selected = sort(unique(Data_case_raw()$症例.基本情報.性別.名称.)), multiple = TRUE)
})
output$select_minimum_pts = renderUI({
  numericInput("minimum_pts", "Minimum patients for each histology",
               value = 1,
               min = 1,
               max = nrow((Data_case_raw() %>% dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID))),
               step = 1)
})

output$select_eps = renderUI({
  numericInput("eps", "Distance value for DBSCAN clustering",
               value = 1.0,
               min = 0,
               max = 10,
               step = 0.1)
})
output$select_panel = renderUI({
  pickerInput("panel", "Filter by panel",
              choices = sort(unique(Data_case_raw()$症例.検体情報.パネル.名称.)),
              selected = sort(unique(Data_case_raw()$症例.検体情報.パネル.名称.)), multiple = TRUE)
})
output$select_year = renderUI({
  pickerInput("year", "Filter by test-year",
              choices = sort(unique(as.POSIXlt(as.Date(Data_case_raw()$症例.管理情報.登録日))$year + 1900)),
              selected = sort(unique(as.POSIXlt(as.Date(Data_case_raw()$症例.管理情報.登録日))$year + 1900)), multiple = TRUE)
})
output$select_age = renderUI({
  sliderInput(inputId = "age", label = "Age for analysis",
              value = c(
                min(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE),
                max(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE)),
              min = min(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE),
              max = max(Data_case_raw()$症例.基本情報.年齢, na.rm = TRUE),
              step = 1)
})
output$select_mid_age = renderUI({
  tmp = Data_case_raw() %>% dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.年齢)
  sliderInput(inputId = "mid_age", label = "Threshold age for analysis",
              value = round(median(tmp$症例.基本情報.年齢, na.rm = TRUE) / 5, digits = 0) * 5,
              min = min(tmp$症例.基本情報.年齢, na.rm = TRUE),
              max = max(tmp$症例.基本情報.年齢, na.rm = TRUE),
              step = 1)
})

output$select_PS = renderUI({
  pickerInput("PS", "Filter by performance status",
              choices = sort(unique(Data_case_raw()$症例.背景情報.ECOG.PS.名称.)),
              selected = sort(unique(Data_case_raw()$症例.背景情報.ECOG.PS.名称.)), multiple = TRUE)
})
output$select_smoking = renderUI({
  pickerInput("smoking", "Filter by smoking status",
              choices = sort(unique(Data_case_raw()$症例.背景情報.喫煙歴有無.名称.)),
              selected = sort(unique(Data_case_raw()$症例.背景情報.喫煙歴有無.名称.)), multiple = TRUE)
})
output$select_gene = renderUI({
  req(Data_report_raw())
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("gene", "Genes of interest",
              candidate,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = unique(names(sort(table(Data_report_raw()$Hugo_Symbol),
                                           decreasing = T)))[1],
              multiple = TRUE)
})
output$select_oncoprint_gene = renderUI({
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("oncoprint_gene", "Genes to be preferentially displayed",
              candidate,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})
output$select_drug_multi_gene = renderUI({
  req(Data_report_raw())
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("drug_multi_gene", "Genes for multivariable analyses",
              candidate,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = unique(names(sort(table(Data_report_raw()$Hugo_Symbol),
                                           decreasing = T)))[1:5],
              multiple = TRUE)
})
output$select_lolliplot_gene = renderUI({
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  candidate = candidate[!candidate %in% c("TMB", "MSI")]
  gene_candidate = unique(names(sort(table(Data_report_raw()$Hugo_Symbol),
                                     decreasing = T)))
  gene_candidate = gene_candidate[gene_candidate %in% candidate]
  gene_candidate = gene_candidate[!gene_candidate %in% c("TMB", "MSI")]

  pickerInput("lolliplot_gene", "Gene for lolliplot",
              options = list(`live-search`=TRUE),
              choices = candidate,
              selected = gene_candidate[1],
              multiple = FALSE)
})
output$select_lolliplot_no = renderUI({
  tmp = Data_case_raw() %>% dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, 症例.基本情報.年齢)
  sliderInput(inputId = "lolliplot_no", label = "Threshold mutation count for lolliplot",
              value = 5,
              min = 1,
              max = 500,
              step = 1)
})


output$select_gene_group_1 = renderUI({
  req(input$gene)
  pickerInput("gene_group_1", "Gene-set 1 of interest",
              input$gene,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})
output$select_gene_group_2 = renderUI({
  req(input$gene_group_1)
  pickerInput("gene_group_2", "Gene-set 2 of interest",
              input$gene[-which(input$gene %in% input$gene_group_1)],
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})
output$select_table_summary_1_gene_group_1 = renderUI({
  req(Data_report_raw())
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("table_summary_1_gene_group_1", "Gene-set 1 of interest",
              candidate,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})
output$select_table_summary_1_gene_group_2 = renderUI({
  req(Data_report_raw())
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("table_summary_1_gene_group_2", "Gene-set 2 of interest",
              candidate[-which(candidate %in% input$table_summary_1_gene_group_1)],
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})

output$select_special_gene = renderUI({
  req(Data_report_raw())
  pickerInput("special_gene", "Gene to analyze",
              c("", sort(unique(Data_report_raw()$Hugo_Symbol))),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = FALSE)
})
output$select_special_gene_mutation_1 = renderUI({
  req(input$special_gene)
  pickerInput("special_gene_mutation_1", "Variants 1", sort(unique(
    (Data_report_raw() %>% dplyr::filter(Hugo_Symbol == input$special_gene | (input$special_gene != "" & str_detect(Hugo_Symbol, paste0(input$special_gene, "_")))))$amino.acid.change)),
    options = list(`actions-box` = TRUE, `live-search`=TRUE),
    multiple = TRUE)
})
output$select_special_gene_mutation_1_name = renderUI({
  req(input$special_gene_mutation_1)
  textInput(inputId = "special_gene_mutation_1_name", label = "Name for variants 1",
            value = input$special_gene_mutation_1[[1]])
})
output$select_special_gene_mutation_2 = renderUI({
  req(input$special_gene_mutation_1)
  tmp_mut = sort(unique(
    (Data_report_raw() %>% dplyr::filter(Hugo_Symbol == input$special_gene | (input$special_gene != "" & str_detect(Hugo_Symbol, paste0(input$special_gene, "_")))))$amino.acid.change))
  pickerInput("special_gene_mutation_2", "Variants 2", tmp_mut[!tmp_mut %in% input$special_gene_mutation_1],
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              multiple = TRUE)
})
output$select_special_gene_mutation_2_name = renderUI({
  req(input$special_gene_mutation_2)
  textInput(inputId = "special_gene_mutation_2_name", label = "Name for variants 2",
            value = input$special_gene_mutation_2[[1]])
})
output$select_gene_no = renderUI({
  req(input$gene)
  numericInput(inputId = "gene_no", label = "Gene number to be automatically analyzed",
               value = min(GENE_NO_THRESHOLD, length(sort(unique(Data_report_raw()$Hugo_Symbol)))),
               min = max(1, length(input$gene)),
               max = length(sort(unique(Data_report_raw()$Hugo_Symbol))),
               step = 1)
})
output$select_oncoprint_gene_no = renderUI({
  numericInput(inputId = "oncoprint_gene_no", label = "Gene number for oncoprint",
               value = min(GENE_NO_THRESHOLD, length(sort(unique(Data_report_raw()$Hugo_Symbol)))),
               min = max(1, length(input$gene)),
               max = length(sort(unique(Data_report_raw()$Hugo_Symbol))),
               step = 1)
})


output$select_patho = renderUI({
  radioButtons(
    "patho", "Variants for analysis",
    choices = c("Only pathogenic muts", "All muts"),
    #selected = "All muts")
    selected = "Only pathogenic muts")
})
output$select_oncoprint_option = renderUI({
  radioButtons(
    "oncoprint_option", "Oncoprint setting",
    choices = c("Sort by mutation frequency", "Sort by pathology", "Sort by age", "Sort by metastasis status at CGP"),
    selected = "Sort by mutation frequency")
})
output$select_target_line = renderUI({
  pickerInput("target_line", "CTx lines to analyze",
              choices = c("1","2","3","4", "5~"),
              selected =  c("1","2","3","4", "5~"),
              options = list(`actions-box` = TRUE),
              multiple = TRUE)
})
output$select_gene_group_analysis = renderUI({
  radioButtons(
    "gene_group_analysis", "Case selection based on the mutations",
    choices = c("Only cases with mutations in the gene set are analyzed",
                "Only cases without mutations in the gene set are analyzed",
                "All cases"),
    selected = "All cases")
})
output$select_special_gene_annotation = renderUI({
  radioButtons(
    "special_gene_annotation", "Significance of the genes",
    choices = c("Only pathogenic mutations are included in the analysis",
                "This gene is also analyzed for VUS and Benign"),
    selected = "Only pathogenic mutations are included in the analysis")
})
output$select_special_gene_independent = renderUI({
  radioButtons(
    "special_gene_independent", "Treat specified variants independently?",
    choices = c("Treat as independent genes in analyses",
                "No"),
    selected = "Treat as independent genes in analyses")
})
output$select_RMST_CGP = renderUI({
  numericInput("RMST_CGP", "Timing for RMST measurement in survival analysis after CGP (years)",
               value = 2,
               min = 1,
               max = 5,
               step = 1)
})
output$select_RMST_CTx = renderUI({
  numericInput("RMST_CTx", "Timing for RMST measurement in survival analysis after CTx (years)",
               value = 2,
               min = 1,
               max = 10,
               step = 1)
})
output$select_RMST_drug = renderUI({
  numericInput("RMST_drug", "Timing for RMST measurement in drug analysis (months)",
               value = 12,
               min = 1,
               max = 60,
               step = 1)
})
# output$select_URL = renderUI({
#   textInput(inputId = "URL", "Hash ID list to be excluded, URL", value = "http://google.com", width = NULL, placeholder = NULL)
# })
output$select_new_analysis = renderUI({
  radioButtons(
    "new_analysis", "Analyze with new dataset",
    choices = c("Yes, input new csv files",
                "No, use the previous dataset"),
  # selected = "Yes, input new csv files")
    selected = ifelse(CCAT_FLAG, "No, use the previous dataset", "Yes, input new csv files"))
})
output$select_intermediate_file = renderUI({
  radioButtons(
    "intermediate_file", "Save intermediate files",
    choices = c("Yes (not recommended)",
                "No"),
    selected = "No")
})
output$select_histology_detail = renderUI({
  radioButtons(
    "histology_detail", "Analyze without detailed histology",
    choices = c("Yes, use oncotree 1st level",
                "No"),
    selected = "No")
})
output$select_MSI = renderUI({
  radioButtons(
    "MSI", "Analyze MSI-PCR test results",
    choices = c("Yes, additional analysis in survival after CGP and drug response",
                "No"),
    selected = "No")
})
output$select_MMR = renderUI({
  radioButtons(
    "MMR", "Analyze MMR-IHC test results",
    choices = c("Yes, additional analysis in survival after CGP and drug response",
                "No"),
    selected = "No")
})
output$select_HER2 = renderUI({
  radioButtons(
    "HER2", "Analyze HER2-IHC test results",
    choices = c("Yes, additional analysis in survival after CGP and drug response",
                "No"),
    selected = "No")
})
output$select_PD_L1 = renderUI({
  radioButtons(
    "PD_L1", "Analyze PD_L1-IHC test results, lung only",
    choices = c("Yes, additional analysis in survival after CGP and drug response",
                "No"),
    selected = "No")
})
output$select_AE_analysis = renderUI({
  radioButtons(
    "AE_analysis", "Analyze drug effect and adverse effect",
    choices = c("Yes", "No"),
    selected = "No")
})
output$select_fusion = renderUI({
  radioButtons(
    "fusion", "How to analyze fusion genes",
    choices = c("Split into two and treat as mutation in each gene",
                "Split into two and treat like EML-fusion, ALK-fusion",
                "All fusion genes are named as 'fusion_genes'",
                "As is"),
    selected = "Split into two and treat as mutation in each gene")
})
output$select_amplification = renderUI({
  radioButtons(
    "amplification", "How to analyze gene amplification/loss",
    choices = c("Treat as indipendent gene",
                "Amplification/loss is a kind of mutation"),
    selected = "Amplification/loss is a kind of mutation")
})
output$select_T_N = renderUI({
  radioButtons(
    "T_N", "Analysis for germline/somatic origin",
    choices = c("Only somatic for T/N panel", "All muts", "Only germline"),
    selected = "All muts")
})
output$download_ID_histology = downloadHandler(
  filename = "ID_histology.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/ID_histology.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_test_clinical_data = downloadHandler(
  filename = "sample_clinical.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/sample_clinical.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
# output$download_treatment_reach_data = downloadHandler(
#   filename = "treatment_reach_summary.csv",
#   content = function(file) {
#     csv_file = read.csv(header = TRUE, file(app_path("source/treatment_reach_summary.csv"),
#                                             encoding='UTF-8-BOM'))
#     write_excel_csv(csv_file, file)
#   }
# )
output$download_test_report_data = downloadHandler(
  filename = "sample_report.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/sample_report.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_ID_drug_pre_CGP = downloadHandler(
  filename = "ID_drug_pre_CGP.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/ID_drug_pre_CGP.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_ID_drug_post_CGP = downloadHandler(
  filename = "ID_drug_post_CGP.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/ID_drug_post_CGP.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_drug_rename = downloadHandler(
  filename = "table_drug_rename.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_drug_rename.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_drug_combination_rename = downloadHandler(
  filename = "table_drug_combination_rename.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_drug_rename.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_regimen_rename = downloadHandler(
  filename = "table_regimen_rename.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_regimen_rename.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_histology_rename = downloadHandler(
  filename = "table_histology_rename.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_histology_rename.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_mutation_rename = downloadHandler(
  filename = "table_mutation_rename.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_mutation_rename.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$download_reannotation = downloadHandler(
  filename = "table_reannotation.csv",
  content = function(file) {
    csv_file = read.csv(header = TRUE, file(app_path("source/table_reannotation.csv"),
                                            encoding='UTF-8-BOM'))
    write_excel_csv(csv_file, file)
  }
)
output$select_drug_volcano = renderUI({
  radioButtons(
    "drug_volcano", "Only volcano plot analysis",
    choices = c("Yes",
                "No, all analyses"),
    selected = "No, all analyses")
})
output$select_ToT_TTF = renderUI({
  radioButtons(
    "ToT_TTF", "Drug treatment period to analyze",
    choices = c("Time on treatment" = "ToT",
                "Time to progression (not reliable)" = "TTF",
                "Time to death" = "TTD",
                "Time to adverse effect" = "TTAE"),
    selected = "ToT")
})
output$select_volcano_pathology = renderUI({
  radioButtons(
    "volcano_pathology", "Considering pathology in volcano plot analysis",
    choices = c("Yes",
                "No"),
    selected = "Yes")
})
output$select_machine_learning = renderUI({
  radioButtons(
    "machine_learning", "Perform machine learning in CGP benefit prediction, time consuming",
    choices = c("Yes",
                "No"),
    selected = "No")
})
output$select_importance = renderUI({
  radioButtons(
    "importance", "Importance of factors in machine learning",
    choices = c("SHAP value (heavy analysis)",
                "Importance (very fast)"),
    selected = "Importance (very fast)")
})
output$select_clustering = renderUI({
  radioButtons(
    "clustering", "Clustering cohort",
    choices = c("Fixed to the whole C-CAT cohort",
                "Fixed to the 1st analysis",
                "Vabiable according to selected cohort"),
    selected = "Fixed to the whole C-CAT cohort")
})
output$select_drug_other = renderUI({
  radioButtons(
    "drug_other", "Rename unselected drugs as others",
    choices = c("Yes",
                "No"),
    selected = "No")
})

output$select_nomogram_max_samples = renderUI({
  numericInput("nomogram_max_samples", "Max nomogram cases (before filtering)",
               value = 10000,
               min = 5000,
               max = 200000,
               step = 5000)
})
output$select_figure_survival_CGP_5_1_max_samples = renderUI({
  numericInput("figure_survival_CGP_5_1_max_samples", "Max regression analysis cases (before filtering)",
               value = 5000,
               min = 5000,
               max = 200000,
               step = 5000)
})
output$select_figure_Bayes_max_samples = renderUI({
  numericInput("figure_Bayes_max_samples", "Max Bayesian analysis cases (before filtering)",
               value = 2000,
               min = 0,
               max = 200000,
               step = 1000)
})
output$select_figure_Drug_Bayes_max_samples = renderUI({
  numericInput("figure_Drug_Bayes_max_samples", "Max Bayesian analysis cases (before filtering)",
               value = 2000,
               min = 1000,
               max = 200000,
               step = 1000)
})
output$select_oncoprint_max_samples = renderUI({
  numericInput("oncoprint_max_samples", "Max oncoplot cases",
               value = 5000,
               min = 5000,
               max = 200000,
               step = 5000)
})
output$select_AE_max_samples = renderUI({
  numericInput("AE_max_samples", "Max adverse effect analyzed cources",
               value = 5000,
               min = 5000,
               max = 200000,
               step = 5000)
})
output$select_TMB_threshold = renderUI({
  numericInput("TMB_threshold", "Threshold for TMB",
               value = 10,
               min = 1,
               max = 30,
               step = 1)
})
output$select_color_var_surv_CTx_1 = renderUI({
  radioButtons(
    inputId = "color_var_surv_CTx_1",
    label = "Survival analysis start date",
    choices = c("1st-line CTx start date" = "1L_CTx",
                "2nd-line CTx start date" = "2L_CTx",
                "Diagnosis date" = "diagnosis",
                "CGP test date" = "CGP"
    ),
    selected = "1L_CTx")
})
output$select_color_var_surv_CTx_3 = renderUI({
  radioButtons(
    inputId = "color_var_surv_CTx_3",
    label = "Risk-set adjustment for left-truncation bias",
    choices = c("Yes", "No"),
    selected = "Yes")
})

output$select_table_summary = renderUI({
  radioButtons(
    inputId = "table_summary",
    label = "Split table columns by",
    choices = c("No split, all patients in one cohort" = "all",
                "Gene variant: Please specify Gene-set 1/2 below." = "gene",
                "Panel" = "panel",
                "Diagnosis" = "diagnosis"
    ),
    selected = "all")
})
output$select_table_summary_no = renderUI({
  numericInput("table_summary_no", "Max histology subtype number (if exceeded, OncoTree first level displayed)",
               value = 20,
               min = 1,
               max = 100,
               step = 1)
})
output$select_table_summary_layout = renderUI({
  radioButtons(
    "table_summary_layout", "Display percentage",
    choices = c("Yes", "No"),
    selected = "No")
})


output$select_multiple_used_drug = renderUI({
  radioButtons(
    inputId = "multiple_used_drug",
    label = "If the same drug was administered multiple times to a single patient",
    choices = c("All treatments were included in the analysis",
                "Only the first administration of the treatment was included in the analysis"),
    selected = "All treatments were included in the analysis")
})
output$input_shape_llogis = renderUI({
  numericInput("shape_llogis", "Shape of log-logistic distribution followed by survival after chemotherapy induction",
               value = 2,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_scale_llogis = renderUI({
  numericInput("scale_llogis", "Scale of log-logistic distribution followed by survival after chemotherapy induction",
               value = 3,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_shape_CGP = renderUI({
  numericInput("shape_CGP", "Shape of log-logistic distribution followed by survival after CGP testing",
               value = 3,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_scale_CGP = renderUI({
  numericInput("scale_CGP", "Scale of log-logistic distribution followed by survival after CGP testing",
               value = 1,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_shape_weibull = renderUI({
  numericInput("shape_weibull", "Shape of Weibull distribution followed by the timing of CGP",
               value = 2,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_scale_weibull = renderUI({
  numericInput("scale_weibull", "Scale of Weibull distribution followed by the timing of CGP",
               value = 3,
               min = 0.1,
               max = 10,
               step = 0.1)
})
output$input_test_patients = renderUI({
  numericInput("test_patients", "Test patients",
               value = 1000,
               min = 100,
               max = 5000,
               step = 100)
})
output$input_censor_rate = renderUI({
  numericInput("censor_rate", "Censoring rate after CGP",
               value = 0.5,
               min = 0,
               max = 1,
               step = 0.1)
})
output$select_figure_drug_evidence_var = renderUI({
  radioButtons(
    inputId = "figure_drug_evidence_var",
    label = "Grouped by",
    choices = c("Cluster" = "cluster",
                "Histology" = "Cancers"),
    selected = "Cancers")
})
output$select_color_var_cluster = renderUI({
  radioButtons(
    inputId = "color_var_cluster",
    label = "Color by",
    choices = c("Age" = "YoungOld",
                "Histology" = "Cancers",
                "Cluster" = "cluster",
                "Treatment option recommended" = "EP_option",
                "Recommended treatment done" = "EP_treat"),
    selected = "cluster")
})


output$select_Bayes_basic_var = renderUI({
  radioButtons(
    inputId = "Bayes_basic_var",
    label = "Analysis period to final observation from",
    choices = c("Diagnosis date" = "diagnosis",
                "1L-CTx initiation" = "1L",
                "2L-CTx initiation" = "2L"),
    selected = "1L")
})
output$select_Bayes_prediction_var = renderUI({
  radioButtons(
    inputId = "Bayes_prediction_var",
    label = "Analysis period to final observation from",
    choices = c("Diagnosis date" = "diagnosis",
                "1L-CTx initiation" = "1L",
                "2L-CTx initiation" = "2L"),
    selected = "1L")
})
output$select_Bayes_diagnosis_var = renderUI({
  radioButtons(
    inputId = "Bayes_diagnosis_var",
    label = "Analysis period to final observation from",
    choices = c("Diagnosis date" = "diagnosis",
                "1L-CTx initiation" = "1L",
                "2L-CTx initiation" = "2L"),
    selected = "1L")
})
output$select_Bayes_diagnosis_no = renderUI({
  numericInput("Bayes_diagnosis_no",
               label = "No of diagnoses to be analyzed",
               value = 3,
               min = 1,
               max = 12,
               step = 1)
})
output$select_Bayes_mutation_var = renderUI({
  radioButtons(
    inputId = "Bayes_mutation_var",
    label = "Analysis period to final observation from",
    choices = c("Diagnosis date" = "diagnosis",
                "1L-CTx initiation" = "1L",
                "2L-CTx initiation" = "2L"),
    selected = "1L")
})
output$select_Bayes_mutation_no = renderUI({
  numericInput("Bayes_mutation_no",
               label = "No of genes to be analyzed",
               value = 5,
               min = 1,
               max = 30,
               step = 1)
})




output$select_color_var_surv_CGP = renderUI({
  radioButtons(
    inputId = "color_var_surv_CGP",
    label = "Grouped by",
    choices = c("Entire cohort" = "entire",
                "Treatment option" = "treat_option",
                "Only pationts with treatment option" = "treat_option_pos",
                "Treatment option and treatment after CGP" = "treat_group",
                "Treatment after CGP in detail" = "treat_group_2",
                "Treatment after CGP, simple" = "treat_group_3",
                "Treatment lines before CGP" = "treat_group_4",
                "Best treatment effect before CGP" = "treat_group_5",
                "Recommended vs not-recommended treatment" = "treat_group_6",
                "ECOG Performance status" = "treat_group_7",
                "ECOG Performance status (unknown patients excluded)" = "treat_group_8",
                "Diagnosis" = "treat_group_9",
                "Panel" = "treat_group_10",
                "Best evidence level" = "treat_group_11",
                "Test year" = "treat_group_12"
    ),
    selected = "entire")
})

output$select_propensity_survival_CGP = renderUI({
  pickerInput(
    inputId = "propensity_survival_CGP",
    label = "Propensiry score matching with clinical factors (if any)",
    choices = c("Panel" = "Panel",
                "ECOG Performance status" = "PS",
                "Diagnosis" = "Cancers",
                "Mutational clustering" = "cluster",
                "Treatment option" = "EP_option",
                "Sex" = "Sex",
                "Age" = "YoungOld",
                "Best evidence level" = "Best_Evidence_Level",
                "Treatment lines before CGP" = "CTx_lines_before_CGP",
                "Best treatment effect before CGP" = "pre_CGP_best_RECIST"
    ),
    selected = NULL,
    multiple = TRUE)
})
output$select_IPW_survival_CGP = renderUI({
  pickerInput(
    inputId = "IPW_survival_CGP",
    label = "Inverse probability weighting with clinical factors (if any)",
    choices = c("Panel" = "Panel",
                "ECOG Performance status" = "PS",
                "Diagnosis" = "Cancers",
                "Mutational clustering" = "cluster",
                "Treatment option" = "EP_option",
                "Sex" = "Sex",
                "Age" = "YoungOld",
                "Best evidence level" = "Best_Evidence_Level",
                "Treatment lines before CGP" = "CTx_lines_before_CGP",
                "Best treatment effect before CGP" = "pre_CGP_best_RECIST"
    ),
    selected = NULL,
    multiple = TRUE)
})
output$select_IPW_threshold = renderUI({
  numericInput("IPW_threshold",
               label = "Exclusion threshold for very heavy weight",
               value = 10,
               min = 1,
               max = 100,
               step = 1)
})


output$select_gene_survival_reach_1_Best_Evidence_Level = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level)
  pickerInput("gene_survival_reach_1_Best_Evidence_Level", "Best Evidence Level",
              OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_Year = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Year)
  pickerInput("gene_survival_reach_1_Year", "Test year",
              OUTPUT_DATA$figure_surv_CGP_candidate_Year,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_M = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_meta)
  pickerInput("gene_survival_reach_1_M", "Metastatic site",
              OUTPUT_DATA$figure_surv_CGP_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_P = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_PS)
  pickerInput("gene_survival_reach_1_P", "Performance status at CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_PS,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_C = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_cluster)
  pickerInput("gene_survival_reach_1_C", "Mutation cluster",
              OUTPUT_DATA$figure_surv_CGP_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_H = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Histology)
  pickerInput("gene_survival_reach_1_H", "Histology",
              OUTPUT_DATA$figure_surv_CGP_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_S = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Sex)
  pickerInput("gene_survival_reach_1_S", "Sex",
              OUTPUT_DATA$figure_surv_CGP_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_A = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Age)
  pickerInput("gene_survival_reach_1_A", "Age",
              OUTPUT_DATA$figure_surv_CGP_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_R = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_RECIST)
  pickerInput("gene_survival_reach_1_R", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_L = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_lines)
  pickerInput("gene_survival_reach_1_L", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_P_1 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_reach_1_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_P_2 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_reach_1_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_W = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_reach_1_W", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_D = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_reach_1_D", "Drugs used before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_ND = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_reach_1_ND", "Drugs not used before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_EP_option = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_EP_option)
  pickerInput("gene_survival_reach_1_EP_option", "Genome-matched therapy recommendation",
              OUTPUT_DATA$figure_surv_CGP_candidate_EP_option,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_reach_1_Panel = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Panel)
  pickerInput("gene_survival_reach_1_Panel", "Panel",
              OUTPUT_DATA$figure_surv_CGP_candidate_Panel,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})

output$select_gene_survival_CGP_1_M = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_meta)
  pickerInput("gene_survival_CGP_1_M", "Metastatic site",
              OUTPUT_DATA$figure_surv_CGP_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_M = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_meta)
  pickerInput("gene_survival_CGP_2_M", "Metastatic site",
              OUTPUT_DATA$figure_surv_CGP_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_Panel = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Panel)
  pickerInput("gene_survival_CGP_1_Panel", "CGP panel",
              OUTPUT_DATA$figure_surv_CGP_candidate_Panel,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_Panel = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Panel)
  pickerInput("gene_survival_CGP_2_Panel", "CGP panel",
              OUTPUT_DATA$figure_surv_CGP_candidate_Panel,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_Best_Evidence_Level = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level)
  pickerInput("gene_survival_CGP_1_Best_Evidence_Level", "Best Evidence Level",
              OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_Best_Evidence_Level = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level)
  pickerInput("gene_survival_CGP_2_Best_Evidence_Level", "Best Evidence Level",
              OUTPUT_DATA$figure_surv_CGP_candidate_Best_Evidence_Level,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_Year = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Year)
  pickerInput("gene_survival_CGP_1_Year", "Test year",
              OUTPUT_DATA$figure_surv_CGP_candidate_Year,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_Year = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Year)
  pickerInput("gene_survival_CGP_2_Year", "Test year",
              OUTPUT_DATA$figure_surv_CGP_candidate_Year,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})


output$select_gene_survival_CGP_1_P = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_PS)
  pickerInput("gene_survival_CGP_1_P", "Performance status at CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_PS,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_P = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_PS)
  pickerInput("gene_survival_CGP_2_P", "Performance status at CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_PS,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_C = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_cluster)
  pickerInput("gene_survival_CGP_1_C", "Mutation cluster",
              OUTPUT_DATA$figure_surv_CGP_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_C = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_cluster)
  pickerInput("gene_survival_CGP_2_C", "Mutation cluster",
              OUTPUT_DATA$figure_surv_CGP_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_H = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Histology)
  pickerInput("gene_survival_CGP_1_H", "Histology",
              OUTPUT_DATA$figure_surv_CGP_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_H = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Histology)
  pickerInput("gene_survival_CGP_2_H", "Histology",
              OUTPUT_DATA$figure_surv_CGP_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_S = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Sex)
  pickerInput("gene_survival_CGP_1_S", "Sex",
              OUTPUT_DATA$figure_surv_CGP_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_S = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Sex)
  pickerInput("gene_survival_CGP_2_S", "Sex",
              OUTPUT_DATA$figure_surv_CGP_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_A = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Age)
  pickerInput("gene_survival_CGP_1_A", "Age",
              OUTPUT_DATA$figure_surv_CGP_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_A = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_Age)
  pickerInput("gene_survival_CGP_2_A", "Age",
              OUTPUT_DATA$figure_surv_CGP_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_R = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_RECIST)
  pickerInput("gene_survival_CGP_1_R", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_R = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_RECIST)
  pickerInput("gene_survival_CGP_2_R", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_EP_option = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_EP_option)
  pickerInput("gene_survival_CGP_1_EP_option", "Genome-matched therapy recommendation",
              OUTPUT_DATA$figure_surv_CGP_candidate_EP_option,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_EP_option = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_EP_option)
  pickerInput("gene_survival_CGP_2_EP_option", "Genome-matched therapy recommendation",
              OUTPUT_DATA$figure_surv_CGP_candidate_EP_option,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_EP_treat = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_EP_treat)
  pickerInput("gene_survival_CGP_1_EP_treat", "Treatment status after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_EP_treat,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_EP_treat = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_EP_treat)
  pickerInput("gene_survival_CGP_2_EP_treat", "Treatment status after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_EP_treat,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})

output$select_gene_survival_CGP_1_L = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_lines)
  pickerInput("gene_survival_CGP_1_L", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_L = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_lines)
  pickerInput("gene_survival_CGP_2_L", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_P_1 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_1_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = OUTPUT_DATA$figure_surv_CGP_Top_gene[1],
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_P_2 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_1_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_W = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_1_W", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_D = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_CGP_1_D", "Drugs used after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_1_ND = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_CGP_1_ND", "Drugs not used after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_P_1 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_2_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_P_2 = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_2_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_W = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_genes)
  pickerInput("gene_survival_CGP_2_W", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_CGP_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_D = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_CGP_2_D", "Drugs used after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_CGP_2_ND = renderUI({
  req(OUTPUT_DATA$figure_surv_CGP_candidate_drugs)
  pickerInput("gene_survival_CGP_2_ND", "Drugs not used after CGP",
              OUTPUT_DATA$figure_surv_CGP_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_figure_survival_CGP_5_1_var = renderUI({
  radioButtons(
    inputId = "figure_survival_CGP_5_1_var",
    label = "Multivariable logistic regression analysis, considering gene mutations",
    choices = c("None" = "None",
                "Frequently mutated genes" = "gene",
                "Mutation-based cluster" = "cluster"),
    selected = NULL)
})
output$select_prediction_var = renderUI({
  req(input$evidence_level)
  tmp_text = paste0("Variants with evidence level ", paste(input$evidence_level, collapse = "/"))
  choices <- c(
    "Genome-matched therapy reach, pre-CGP information" = "",
    "Treatment option recommendation, pre-CGP information" = "_rec",
    setNames("_GMT", tmp_text)
  )
  radioButtons(
    inputId = "prediction_var",
    label = "Prediction target",
    choices = choices,
    selected = "")
})

output$input_age = renderUI({
  numericInput("predict_age", "Age",
               value = 60,
               min = 1,
               max = 120,
               step = 1)
})
output$input_sex = renderUI({
  radioButtons(
    "predict_Sex", "Sex",
    choices = c("Male", "Female"),
    selected = "Male")
})
output$input_Smoking_history = renderUI({
  radioButtons(
    "predict_Smoking_history", "Smoking history",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_Alcoholic_history = renderUI({
  radioButtons(
    "predict_Alcoholic_history", "Alcoholic history",
    choices = c("Yes", "No"),
    selected = "No")
})
# output$input_Enroll_date = renderUI({
#   if(file.exists(paste0(tempdir(), "/Enroll_date_data", input$prediction_var, ".rda"))){
#     load(paste0(tempdir(), "/Enroll_date_data", input$prediction_var, ".rda"))
#   }
#   object_name <- paste0("Enroll_date_data", input$prediction_var)
#   req(exists(object_name))
#   dynamic_organ_data <- get(object_name)
#   pickerInput("predict_Enroll_date", "Enroll year",
#               choices = dynamic_organ_data$Enroll_date_type,
#               selected = dynamic_organ_data$Enroll_date_type[1], multiple = FALSE)
# })
# output$input_Enroll_date_raw = renderUI({
#   pickerInput(
#     "predict_Enroll_date", "Enroll year, machine learning",
#     choices = c("2019", "2020", "2021", "2022", "2023", "2024", "2025"),
#     selected = "2022", multiple = FALSE)
# })
output$input_time_diagnosis_enroll = renderUI({
  radioButtons(
    "predict_time_diagnosis_enroll", "Diagnosis to CGP test",
    choices = c(">1-year", "<=1-year"),
    selected = ">1-year")
})
output$input_Lymph_met = renderUI({
  radioButtons(
    "predict_Lymph_met", "Lymph node metastasis",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_Brain_met = renderUI({
  radioButtons(
    "predict_Brain_met", "Brain metastasis",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_Lung_met = renderUI({
  radioButtons(
    "predict_Lung_met", "Lung metastasis",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_Bone_met = renderUI({
  radioButtons(
    "predict_Bone_met", "Bone metastasis",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_Liver_met = renderUI({
  radioButtons(
    "predict_Liver_met", "Liver metastasis",
    choices = c("Yes", "No"),
    selected = "No")
})
output$input_organ = renderUI({
  if(file.exists(paste0(tempdir(), "/Organ_data", input$prediction_var, ".rda"))){
    load(paste0(tempdir(), "/Organ_data", input$prediction_var, ".rda"))
  }
  object_name <- paste0("Organ_data", input$prediction_var)
  req(exists(object_name))
  dynamic_organ_data <- get(object_name)
  pickerInput("predict_organ", "Histology",
              choices = dynamic_organ_data$Organ_type,
              selected = dynamic_organ_data$Organ_type[1], multiple = FALSE)
})
output$input_organ_raw = renderUI({
  if(file.exists(paste0(tempdir(), "/Organ_data", input$prediction_var, ".rda"))){
    load(paste0(tempdir(), "/Organ_data", input$prediction_var, ".rda"))
  }
  object_name <- paste0("Organ_data", input$prediction_var)
  req(exists(object_name))
  dynamic_data <- get(object_name)
  pickerInput("predict_organ_raw", "Histology, machine learning",
              choices = dynamic_data$Organ_type_raw,
              selected = dynamic_data$Organ_type_raw[1], multiple = FALSE)
})
output$input_PS = renderUI({
  if(file.exists(paste0(tempdir(), "/PS_data", input$prediction_var, ".rda"))){
    load(paste0(tempdir(), "/PS_data", input$prediction_var, ".rda"))
    object_name <- paste0("PS_data", input$prediction_var)
    req(exists(object_name))
    dynamic_data <- get(object_name)
    PS_choice = dynamic_data$PS_type
  } else {
    PS_choice = c("0", "1", "2_4")
  }
  radioButtons(
    "predict_PS", "Performance status",
    choices = PS_choice,
    selected = PS_choice[1])
})
output$input_PS_raw = renderUI({
  radioButtons(
    "predict_PS_raw", "Performance status, machine learning",
    choices = c("0", "1", "2_4"),
    selected = c("0", "1", "2_4")[1])
})
output$input_Lines = renderUI({
  if(file.exists(paste0(tempdir(), "/Lines_data", input$prediction_var, ".rda"))){
    load(paste0(tempdir(), "/Lines_data", input$prediction_var, ".rda"))
    object_name <- paste0("Lines_data", input$prediction_var)
    req(exists(object_name))
    dynamic_data <- get(object_name)
    Lines_choice = dynamic_data$Lines_type
  } else {
    Lines_choice = c("0", "1", "2~")
  }
  radioButtons(
    "predict_Lines", "CTx lines before CGP",
    choices = Lines_choice,
    selected = Lines_choice[1])
})
output$input_Lines_raw = renderUI({
  radioButtons(
    "predict_Lines_raw", "CTx lines before CGP, machine learning",
    choices = c("0", "1", "2~"),
    selected = "1")
})
output$input_Best_effect = renderUI({
  radioButtons(
    "predict_Best_effect", "Best CTx effect in previous lines",
    choices = c("CR", "PR", "SD", "PD", "Unknown", "CTx_naive"),
    selected = "SD")
})
output$input_Panel = renderUI({
  if(file.exists(paste0(tempdir(), "/Panel_data", input$prediction_var, ".rda"))){
    load(paste0(tempdir(), "/Panel_data", input$prediction_var, ".rda"))
    object_name <- paste0("Panel_data", input$prediction_var)
    req(exists(object_name))
    dynamic_data <- get(object_name)
    Panel_choice = dynamic_data$Panel_type
  } else {
    Panel_choice = c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP", "Liquid")
  }
  radioButtons(
    "predict_Panel", "Panel",
    choices = Panel_choice,
    selected = Panel_choice[1])
})
output$input_Panel_raw = renderUI({
  radioButtons(
    "predict_Panel_raw", "Panel, machine learning",
    choices = c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP", "Liquid"),
    selected = c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP", "Liquid")[1])
})

output$select_color_var_surv_CTx_2 = renderUI({
  radioButtons(
    inputId = "color_var_surv_CTx_2",
    label = "Grouped by",
    choices = c("Entire cohort" = "entire",
                "Treatment option" = "treat_option",
                "Only pationts with treatment option" = "treat_option_pos",
                "Treatment option and treatment after CGP" = "treat_group",
                "Treatment after CGP in detail" = "treat_group_2",
                "Treatment after CGP, simple" = "treat_group_3",
                "Treatment lines before CGP" = "treat_group_4",
                "Best treatment effect before CGP" = "treat_group_5",
                "Recommended vs not-recommended treatment" = "treat_group_6",
                "ECOG Performance status" = "treat_group_7",
                "ECOG Performance status (unknown patients excluded)" = "treat_group_8",
                "Diagnosis" = "treat_group_9"
    ),
    selected = "entire")
})
output$select_gene_survival_interactive_1_M = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_meta)
  pickerInput("gene_survival_interactive_1_M", "Metastatic site",
              OUTPUT_DATA$figure_surv_interactive_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_M = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_meta)
  pickerInput("gene_survival_interactive_2_M", "Metastatic site",
              OUTPUT_DATA$figure_surv_interactive_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_C = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_cluster)
  pickerInput("gene_survival_interactive_1_C", "Mutation cluster",
              OUTPUT_DATA$figure_surv_interactive_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_C = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_cluster)
  pickerInput("gene_survival_interactive_2_C", "Mutation cluster",
              OUTPUT_DATA$figure_surv_interactive_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_H = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Histology)
  pickerInput("gene_survival_interactive_1_H", "Histology",
              OUTPUT_DATA$figure_surv_interactive_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_H = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Histology)
  pickerInput("gene_survival_interactive_2_H", "Histology",
              OUTPUT_DATA$figure_surv_interactive_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_S = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Sex)
  pickerInput("gene_survival_interactive_1_S", "Sex",
              OUTPUT_DATA$figure_surv_interactive_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_S = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Sex)
  pickerInput("gene_survival_interactive_2_S", "Sex",
              OUTPUT_DATA$figure_surv_interactive_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_A = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Age)
  pickerInput("gene_survival_interactive_1_A", "Age",
              OUTPUT_DATA$figure_surv_interactive_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_A = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Age)
  pickerInput("gene_survival_interactive_2_A", "Age",
              OUTPUT_DATA$figure_surv_interactive_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_R = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_RECIST)
  pickerInput("gene_survival_interactive_1_R", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_R = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_RECIST)
  pickerInput("gene_survival_interactive_2_R", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_L = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_lines)
  pickerInput("gene_survival_interactive_1_L", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_L = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_lines)
  pickerInput("gene_survival_interactive_2_L", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_P_1 = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_1_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = OUTPUT_DATA$figure_surv_interactive_Top_gene[1],
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_P_2 = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_1_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_W = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_1_W", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_D = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_drugs)
  pickerInput("gene_survival_interactive_1_D", "Drugs used in any lines",
              OUTPUT_DATA$figure_surv_interactive_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_P_1 = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_2_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_P_2 = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_2_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_W = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes)
  pickerInput("gene_survival_interactive_2_W", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_D = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_drugs)
  pickerInput("gene_survival_interactive_2_D", "Drugs used in any lines",
              OUTPUT_DATA$figure_surv_interactive_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})




output$select_gene_survival_interactive_1_M_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_meta_Bayes)
  pickerInput("gene_survival_interactive_1_M_Bayes", "Metastatic site",
              OUTPUT_DATA$figure_surv_interactive_candidate_meta_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_M_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_meta_Bayes)
  pickerInput("gene_survival_interactive_2_M_Bayes", "Metastatic site",
              OUTPUT_DATA$figure_surv_interactive_candidate_meta_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_C_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_cluster_Bayes)
  pickerInput("gene_survival_interactive_1_C_Bayes", "Mutation cluster",
              OUTPUT_DATA$figure_surv_interactive_candidate_cluster_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_C_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_cluster_Bayes)
  pickerInput("gene_survival_interactive_2_C_Bayes", "Mutation cluster",
              OUTPUT_DATA$figure_surv_interactive_candidate_cluster_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_H_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Histology_Bayes)
  pickerInput("gene_survival_interactive_1_H_Bayes", "Histology",
              OUTPUT_DATA$figure_surv_interactive_candidate_Histology_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_H_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Histology_Bayes)
  pickerInput("gene_survival_interactive_2_H_Bayes", "Histology",
              OUTPUT_DATA$figure_surv_interactive_candidate_Histology_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_S_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Sex_Bayes)
  pickerInput("gene_survival_interactive_1_S_Bayes", "Sex",
              OUTPUT_DATA$figure_surv_interactive_candidate_Sex_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_S_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Sex_Bayes)
  pickerInput("gene_survival_interactive_2_S_Bayes", "Sex",
              OUTPUT_DATA$figure_surv_interactive_candidate_Sex_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_A_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Age_Bayes)
  pickerInput("gene_survival_interactive_1_A_Bayes", "Age",
              OUTPUT_DATA$figure_surv_interactive_candidate_Age_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_A_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_Age_Bayes)
  pickerInput("gene_survival_interactive_2_A_Bayes", "Age",
              OUTPUT_DATA$figure_surv_interactive_candidate_Age_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_R_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_RECIST_Bayes)
  pickerInput("gene_survival_interactive_1_R_Bayes", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_RECIST_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_R_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_RECIST_Bayes)
  pickerInput("gene_survival_interactive_2_R_Bayes", "Best CTx effect before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_RECIST_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_L_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_lines_Bayes)
  pickerInput("gene_survival_interactive_1_L_Bayes", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_lines_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_L_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_lines_Bayes)
  pickerInput("gene_survival_interactive_2_L_Bayes", "CTx lines before CGP",
              OUTPUT_DATA$figure_surv_interactive_candidate_lines_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_P_1_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_1_P_1_Bayes", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = OUTPUT_DATA$figure_surv_interactive_Top_gene_Bayes[1],
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_P_2_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_1_P_2_Bayes", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_W_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_1_W_Bayes", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_1_D_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_drugs_Bayes)
  pickerInput("gene_survival_interactive_1_D_Bayes", "Drugs used in any lines",
              OUTPUT_DATA$figure_surv_interactive_candidate_drugs_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_P_1_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_2_P_1_Bayes", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_P_2_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_2_P_2_Bayes", "Pathogenic variant in any genes",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_W_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes)
  pickerInput("gene_survival_interactive_2_W_Bayes", "Genes without pathogenic variant",
              OUTPUT_DATA$figure_surv_interactive_candidate_genes_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_interactive_2_D_Bayes = renderUI({
  req(OUTPUT_DATA$figure_surv_interactive_candidate_drugs_Bayes)
  pickerInput("gene_survival_interactive_2_D_Bayes", "Drugs used in any lines",
              OUTPUT_DATA$figure_surv_interactive_candidate_drugs_Bayes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})


output$select_nomogram_type = renderUI({
  req(input$evidence_level)
  tmp_text = paste0("Variants with evidence level ", paste(input$evidence_level, collapse = "/"))
  choices <- c(
    "Genome-matched therapy reach, pre-CGP information" = "EP_treat",
    "Treatment option recommendation, pre-CGP information" = "EP_option",
    setNames("_GMT", tmp_text)
  )
  radioButtons(
    inputId = "nomogram_type",
    label = "Nomogram for",
    choices = choices,
    selected = NULL)
})
output$select_nomogram_type_2 = renderUI({
  req(input$evidence_level)
  tmp_text = paste0("Variants with evidence level ", paste(input$evidence_level, collapse = "/"))
  choices <- c(
    "Genome-matched therapy reach, pre-CGP information" = "EP_treat",
    "Treatment option recommendation, pre-CGP information" = "EP_option",
    setNames("_GMT", tmp_text)
  )
  radioButtons(
    inputId = "nomogram_type_2",
    label = "Logistic regression analysis for",
    choices = choices,
    selected = NULL)
})
output$select_nomogram_type_3 = renderUI({
  req(input$evidence_level)
  tmp_text = paste0("Variants with evidence level ", paste(input$evidence_level, collapse = "/"))
  choices <- c(
    "Genome-matched therapy reach, pre-CGP information" = "EP_treat",
    "Treatment option recommendation, pre-CGP information" = "EP_option",
    setNames("_GMT", tmp_text)
  )
  radioButtons(
    inputId = "nomogram_type_3",
    label = "Decision curve analysis for",
    choices = choices,
    selected = NULL)
})
output$select_nomogram_type_4 = renderUI({
  req(input$evidence_level)
  tmp_text = paste0("Variants with evidence level ", paste(input$evidence_level, collapse = "/"))
  choices <- c(
    "Genome-matched therapy reach, pre-CGP information" = "EP_treat",
    "Treatment option recommendation, pre-CGP information" = "EP_option",
    setNames("_GMT", tmp_text)
  )
  radioButtons(
    inputId = "nomogram_type_4",
    label = "Decision curve analysis for",
    choices = choices,
    selected = NULL)
})

output$select_figure_survival_drug_var = renderUI({
  radioButtons(
    inputId = "figure_survival_drug_var",
    label = "Survival analysis from drug initiation, with selected lines and selected regimens, grouped by",
    choices = c("Mutation status, with selected drugs in the selected lines" = "Mutation_OS",
                "With or without dverse event (G3-5), with selected drugs in the selected lines" = "AE_OS",
                "Regimens in the selected lines, whole selected regimens vs other drugs" = "regimen_OS",
                "Regimens in the selected lines, regimen set-1 vs regimen set-2" = "Drug_OS"
    ),
    selected = NULL)
})

output$select_figure_survival_drug_var_2 = renderUI({
  radioButtons(
    inputId = "figure_survival_drug_var_2",
    label = "Left truncation bias adjustment with",
    choices = c("Bayesian inference" = "Bayes",
                "Risk-set adjustment" = "Number at risk",
                "No adjustment" = "No"
    ),
    selected = "No")
})

output$select_drug_survival_gene = renderUI({
  req(Data_report_raw())
  candidate = sort(unique(c(Data_report_raw()$Hugo_Symbol,
                            paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
                            paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
                            paste0(input$special_gene, "_NOS"))))
  candidate = candidate[!candidate %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
  candidate = candidate[!is.na(candidate)]
  pickerInput("drug_survival_gene", "Genes of interest",
              candidate,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = unique(names(sort(table(Data_report_raw()$Hugo_Symbol),
                                           decreasing = T)))[1],
              multiple = TRUE)
})

output$select_drug_ORR_table_var = renderUI({
  radioButtons(
    inputId = "drug_ORR_table_var",
    label = "Objective response/disease control rate table for",
    choices = c("Mutated genes" = "Mutation",
                "Mutation-based clustering" = "Cluster",
                "Disease types" = "Histology",
                "About the 'Gene to analyze' classified in detail in Settings" = "About the gene classified in detail in Settings"
    ),
    selected = "No")
})

output$select_table_var_volcano_2_AE = renderUI({
  radioButtons(
    inputId = "table_var_volcano_2_AE",
    label = "Multivariable logistic regression analysis, considering gene mutations",
    choices = c("None" = "None",
                "Frequently mutated genes" = "gene",
                "Mutation-based cluster" = "cluster"),
    selected = NULL)
})

output$select_ToT_var_1 = renderUI({
  pickerInput(
    inputId = "ToT_var_1",
    label = "Figures with preset parameters",
    choices = c(#"Concurrent side effects correlation matrix of the selected regimen" = "AE_matrix",
                "Comparison of the treatment duration of the selected drug and the preceding treatment line" = "previous",
                "Treatment duration of all drugs in the selected line (all drugs)" = "selected_line",
                "Treatment duration of all drugs in the selected line, grouped by presence of gene mutation" = "mutation",
                "Treatment duration in the selected regimen, grouped by regimen" = "regimen",
                "Treatment duration in the selected regimen, grouped by detailed regimen" = "regimen_detail",
                "Treatment duration in the selected regimen, grouped by mutation and regimen" = "mutation_regimen",
                "Treatment duration in the selected regimen, grouped by detailed mutation" = "mutation_detail",
                "Treatment duration in the selected regimen, grouped by lymph node metastasis" = "lymph",
                "Treatment duration of all drugs in the selected line, grouped by diagnosis" = "gc_diagnosis",
                "Treatment duration of all drugs in the selected line, grouped by diagnosis with total cohort data" = "gc_diagnosis_all",
                "Treatment duration in the selected regimen, grouped by diagnosis" = "gc_regimen",
                "Treatment duration of all drugs in the selected line, grouped by presence of adverse effect" = "AE",
                "Treatment duration of all drugs in the selected line, grouped by adverse effect within 1-month, 3-month, or later" = "Early_AE",
                "Treatment duration in the selected regimen, grouped by adverse effect" = "gc_AE",
                "Treatment duration in the selected regimen, grouped by adverse effect within 1-month, 3-month, or later" = "gc_EAE"
    ),
    options = list(`actions-box` = TRUE, `live-search`=TRUE),
    multiple = FALSE,
    selected = "previous")
})

output$select_table_var_drug_dataset = renderUI({
  radioButtons(
    inputId = "table_var_drug_dataset",
    label = "Dataset for analysis",
    choices = c("Whole" = "all",
                "Treatment time" = "time",
                "Objective response" = "RECIST",
                "Adverse effect" = "AE"
    ),
    selected = "time")
})

output$select_gene_survival_drug_interactive_1_PD_L1 = renderUI({
  req(OUTPUT_DATA$Data_drug_Data_case_target)
  if("LUNG" %in% OUTPUT_DATA$Data_drug_Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
    candidate_PD_L1 = sort(unique(OUTPUT_DATA$Data_drug_Data_case_target$PD_L1))
    pickerInput("gene_survival_drug_interactive_1_PD_L1", "PD-L1 (lung only)",
              candidate_PD_L1,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
  }
})
output$select_gene_survival_drug_interactive_2_PD_L1 = renderUI({
  req(OUTPUT_DATA$Data_drug_Data_case_target)
  if("LUNG" %in% OUTPUT_DATA$Data_drug_Data_case_target$症例.基本情報.がん種.OncoTree.LEVEL1. & input$PD_L1 != "No"){
    candidate_PD_L1 = sort(unique(OUTPUT_DATA$Data_drug_Data_case_target$PD_L1))
    pickerInput("gene_survival_drug_interactive_2_PD_L1", "PD-L1 (lung only)",
                candidate_PD_L1,
                options = list(`actions-box` = TRUE, `live-search`=TRUE),
                selected = NULL,
                multiple = TRUE)
  }
})
output$select_gene_survival_drug_interactive_1_AEG = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Overall_AE_grage)
  pickerInput("gene_survival_drug_interactive_1_AEG", "Max AE grade",
              OUTPUT_DATA$drug_analysis_candidate_Overall_AE_grage,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_AEG = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Overall_AE_grage)
  pickerInput("gene_survival_drug_interactive_2_AEG", "Max AE grade",
              OUTPUT_DATA$drug_analysis_candidate_Overall_AE_grage,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_AED = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_AE)
  pickerInput("gene_survival_drug_interactive_1_AED", "Detailed adverse effect",
              OUTPUT_DATA$drug_analysis_candidate_AE,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_AED = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_AE)
  pickerInput("gene_survival_drug_interactive_2_AED", "Detailed adverse effect",
              OUTPUT_DATA$drug_analysis_candidate_AE,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_R = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_RECIST)
  pickerInput("gene_survival_drug_interactive_1_R", "RECIST",
              OUTPUT_DATA$drug_analysis_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_R = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_RECIST)
  pickerInput("gene_survival_drug_interactive_2_R", "RECIST",
              OUTPUT_DATA$drug_analysis_candidate_RECIST,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_M = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_meta)
  pickerInput("gene_survival_drug_interactive_1_M", "Metastatic site",
              OUTPUT_DATA$drug_analysis_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_M = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_meta)
  pickerInput("gene_survival_drug_interactive_2_M", "Metastatic site",
              OUTPUT_DATA$drug_analysis_candidate_meta,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_C = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_cluster)
  pickerInput("gene_survival_drug_interactive_1_C", "Mutation cluster",
              OUTPUT_DATA$drug_analysis_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_C = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_cluster)
  pickerInput("gene_survival_drug_interactive_2_C", "Mutation cluster",
              OUTPUT_DATA$drug_analysis_candidate_cluster,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_H = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Histology)
  pickerInput("gene_survival_drug_interactive_1_H", "Histology",
              OUTPUT_DATA$drug_analysis_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_H = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Histology)
  pickerInput("gene_survival_drug_interactive_2_H", "Histology",
              OUTPUT_DATA$drug_analysis_candidate_Histology,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_S = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Sex)
  pickerInput("gene_survival_drug_interactive_1_S", "Sex",
              OUTPUT_DATA$drug_analysis_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_S = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Sex)
  pickerInput("gene_survival_drug_interactive_2_S", "Sex",
              OUTPUT_DATA$drug_analysis_candidate_Sex,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_A = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Age)
  pickerInput("gene_survival_drug_interactive_1_A", "Age",
              OUTPUT_DATA$drug_analysis_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_A = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_Age)
  pickerInput("gene_survival_drug_interactive_2_A", "Age",
              OUTPUT_DATA$drug_analysis_candidate_Age,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_P_1 = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_1_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = OUTPUT_DATA$drug_analysis_Top_gene[1],
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_P_2 = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_1_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_W = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_1_W", "Genes without pathogenic variant",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_D = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_drugs)
  pickerInput("gene_survival_drug_interactive_1_D", "Drugs used in any lines",
              OUTPUT_DATA$drug_analysis_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_AE = renderUI({
  pickerInput("gene_survival_drug_interactive_1_AE", "Adverse effect (G3-G5) with the drugs",
              c("AE (+)", "AE (-)"),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_EAE = renderUI({
  pickerInput("gene_survival_drug_interactive_1_EAE", "Adverse effect within 1-month, 3-month, or later",
              c("Early AE (+)", "Middle AE (+)", "Late AE (+)", "Unknown AE (+)", "AE (-)"),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_P_1 = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_2_P_1", "Pathogenic variant in any genes",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_P_2 = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_2_P_2", "Pathogenic variant in any genes",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_W = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_genes)
  pickerInput("gene_survival_drug_interactive_2_W", "Genes without pathogenic variant",
              OUTPUT_DATA$drug_analysis_candidate_genes,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_D = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_drugs)
  pickerInput("gene_survival_drug_interactive_2_D", "Drugs used in any lines",
              OUTPUT_DATA$drug_analysis_candidate_drugs,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_AE = renderUI({
  pickerInput("gene_survival_drug_interactive_2_AE", "Adverse effect (G3-G5) with the drugs",
              c("AE (+)", "AE (-)"),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_EAE = renderUI({
  pickerInput("gene_survival_drug_interactive_2_EAE", "Adverse effect within 1-month, 3-month, or later",
              c("Early AE (+)", "Middle AE (+)", "Late AE (+)", "Unknown AE (+)", "AE (-)"),
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_1_L = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_lines)
  pickerInput("gene_survival_drug_interactive_1_L", "CTx lines in which drugs used",
              OUTPUT_DATA$drug_analysis_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
output$select_gene_survival_drug_interactive_2_L = renderUI({
  req(OUTPUT_DATA$drug_analysis_candidate_lines)
  pickerInput("gene_survival_drug_interactive_2_L", "CTx line in which drugs used",
              OUTPUT_DATA$drug_analysis_candidate_lines,
              options = list(`actions-box` = TRUE, `live-search`=TRUE),
              selected = NULL,
              multiple = TRUE)
})
