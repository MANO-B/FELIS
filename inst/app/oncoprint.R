oncoprint_logic <- function() {
  analysis_env <- new.env()
  clear_reactive_data(OUTPUT_DATA)
  gc()
  with(analysis_env, {
    withProgress(message = sample(nietzsche)[1], {
      req(input$gene_group_analysis, input$patho)
      Data_case_target = Data_case()
      if(input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed", "Only cases without mutations in the gene set are analyzed")){
        Data_case_target = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in%
                          Data_report()$Tumor_Sample_Barcode)
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
      if(nrow(Data_MAF[Data_MAF$TMB > 30,]) > 0){
        Data_MAF[Data_MAF$TMB > 30,]$TMB = 30
      }

      incProgress(1 / 8)

      Data_oncoprint = Data_case_target %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                      症例.基本情報.がん種.OncoTree.,
                      YoungOld,
                      Lymph_met,
                      Brain_met,
                      Lung_met,
                      Bone_met,
                      Liver_met,
                      Other_met) %>%
        dplyr::distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all=TRUE) %>%
        dplyr::arrange(C.CAT調査結果.基本項目.ハッシュID)
      Data_MAF_target = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in%
                        Data_oncoprint$C.CAT調査結果.基本項目.ハッシュID)
      if(input$patho == "Only pathogenic muts"){
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Evidence_level == "F")
      } else {
        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(!Hugo_Symbol %in% c("TMB", "MSI") | Evidence_level == "F")
      }
      Data_MAF_target_tmp = Data_MAF_target %>%
        dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change, .keep_all = T) %>%
        dplyr::arrange(Tumor_Sample_Barcode, Hugo_Symbol, amino.acid.change)
      incProgress(1 / 8)

      OUTPUT_DATA$oncoprint_Data_oncoprint = Data_oncoprint
      OUTPUT_DATA$oncoprint_Data_MAF_target = Data_MAF_target_tmp
      incProgress(1 / 8)
      incProgress(1 / 8)

      Data_case_target = Data_case_target %>%
        dplyr::select(C.CAT調査結果.基本項目.ハッシュID,
                      症例.基本情報.年齢,
                      #YoungOld,
                      Lymph_met,
                      Brain_met,
                      Lung_met,
                      Bone_met,
                      Liver_met,
                      Other_met,
                      EP_option,
                      EP_treat,
                      症例.基本情報.性別.名称.,
                      症例.背景情報.初回治療前のステージ分類.名称.,
                      症例.背景情報.病理診断名,
                      症例.背景情報.臨床診断名,
                      症例.検体情報.病理診断名,
                      症例.基本情報.がん種.OncoTree.,
                      症例.基本情報.がん種.OncoTree..名称.,
                      症例.基本情報.がん種.OncoTree.LEVEL1.,
                      症例.検体情報.パネル.名称.,
                      症例.背景情報.ECOG.PS.名称.,
                      症例.背景情報.喫煙歴有無.名称.,
                      症例.背景情報.アルコール多飲有無.名称.,
                      症例.背景情報.重複がん有無.異なる臓器..名称.,
                      症例.背景情報.多発がん有無.同一臓器..名称.,
                      症例.がん種情報.登録時転移部位.名称.,
                      症例.検体情報.腫瘍細胞含有割合,
                      症例.検体情報.検体採取部位.名称.,
                      症例.背景情報.診断日,
                      症例.管理情報.登録日,
                      症例.転帰情報.転帰.名称.,
                      症例.転帰情報.最終生存確認日,
                      症例.転帰情報.死亡日,
                      症例.検体情報.検体採取日.腫瘍組織.
        )
      Data_drug = Data_drug_raw()
      filtered_data <- Data_drug %>%
        dplyr::select(ID, 投与開始日, Drug)
      filtered_data = filtered_data %>%
        dplyr::left_join((Data_case_target %>%
                            dplyr::select(C.CAT調査結果.基本項目.ハッシュID, 症例.検体情報.検体採取日.腫瘍組織.) %>%
                            dplyr::distinct( C.CAT調査結果.基本項目.ハッシュID,.keep_all = T)),
                         by=c("ID" = "C.CAT調査結果.基本項目.ハッシュID")) %>%
        dplyr::filter(投与開始日 < 症例.検体情報.検体採取日.腫瘍組織.)
      filtered_data$Drug[is.na(filtered_data$Drug)] = ""
      filtered_data <- filtered_data %>%
        dplyr::filter(Drug != "") %>%
        dplyr::select(ID, Drug)
      summarized_data <- filtered_data %>%
        group_by(ID) %>%
        summarise(
          検査前使用薬剤 = paste(unique(Drug), collapse = ","),
          .groups = "drop"
        )
      Data_case_target <- Data_case_target %>%
        distinct(C.CAT調査結果.基本項目.ハッシュID, .keep_all = TRUE) %>%
        left_join(summarized_data, by = c("C.CAT調査結果.基本項目.ハッシュID" = "ID")) %>%
        dplyr::arrange(C.CAT調査結果.基本項目.ハッシュID)
      Data_case_target$検査前使用薬剤[is.na(Data_case_target$検査前使用薬剤)] <- ""
      incProgress(1 / 8)
      Data_TMB = (Data_MAF %>%
                    dplyr::distinct(Tumor_Sample_Barcode, TMB))
      colnames(Data_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB_raw")
      Data_case_target = left_join(Data_case_target,
                                   Data_TMB,
                                   by = "C.CAT調査結果.基本項目.ハッシュID")
      Data_case_target$TMB_raw[is.na(Data_case_target$TMB_raw)] <- 0
      Data_MAF_target_tmp = Data_MAF_target_tmp %>%
        dplyr::mutate(
          Evidence = case_when(
            Evidence_level == "F" ~ "Oncogenic",
            TRUE ~ ""
          ))
      gene_list = unique(c(input$gene,
                           names(sort(table(Data_MAF_target$Hugo_Symbol),
                                      decreasing = T))))
      gene_list = gene_list[gene_list %in% unique(Data_MAF_target$Hugo_Symbol)]
      gene_list = gene_list[1:min(input$gene_no, length(gene_list))]
      # ② Data_MAF_target_tmp から gene_list に含まれる遺伝子について、
      #    サンプルごとに amino.acid.change をカンマ区切りでまとめる
      gene_data <- Data_MAF_target_tmp %>%
        filter(Hugo_Symbol %in% gene_list) %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        summarise(Variant_mut = paste(amino.acid.change, collapse = ","), .groups = "drop") %>%
        pivot_wider(names_from = Hugo_Symbol,
                    values_from = Variant_mut,
                    values_fill = list(Variant_mut = ""))

      # ③ gene_table の初期化：各サンプル（Data_case_target のハッシュID）を行とし、
      #     gene_list の各遺伝子を列とする。gene_table$target_gene_sequence は別途追加。
      gene_table <- Data_case_target %>%
        select(Tumor_Sample_Barcode = `C.CAT調査結果.基本項目.ハッシュID`)

      # gene_data のサンプルIDと gene_table の ID をキーに left_join
      gene_table <- gene_table %>%
        left_join(gene_data, by = c("Tumor_Sample_Barcode" = "Tumor_Sample_Barcode"))

      # gene_table に不足している列（gene_list に存在するが、gene_data にないサンプル）を空文字に補完
      if (!is.null(gene_list) && length(gene_list) > 0 && !all(is.na(gene_list))) {
        gene_table <- gene_table %>%
          mutate(across(all_of(gene_list), ~ replace_na(.x, ""))) %>%
          mutate(target_gene_sequence = "")
      } else {
        gene_table <- gene_table %>%
          mutate(target_gene_sequence = "")
      }
      # 1. 最頻出の Hugo_Symbol を取得
      if(!is.null(input$special_gene) && input$special_gene != ""){
        gene_name = input$special_gene[[1]]
      } else {
        gene_name = names(sort(table(Data_MAF_target_tmp$Hugo_Symbol),decreasing = T))
        gene_name = gene_name[!is.na(gene_name)]
        gene_name = gene_name[gene_name %in% gene_list][[1]]
      }
      # 2. 該当遺伝子のデータを抽出し Variant_mut を作成
      tmp_mut <- Data_MAF_target_tmp %>%
        filter(Hugo_Symbol == gene_name) %>%
        mutate(Variant_mut = paste0("chr", Chromosome, ":", Start_Position, "_",
                                    Reference_Allele, ">", Tumor_Seq_Allele2, "_", Evidence))
      # 3. サンプルごとに Variant_mut を結合（カンマ区切り）
      mutations_summary <- tmp_mut %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarise(variants = paste(Variant_mut, collapse = ","), .groups = "drop")
      # 4. gene_table に、すでに target_gene_sequence が空なら variants を、そうでなければ既存と結合した文字列に更新
      gene_table <- gene_table %>%
        left_join(mutations_summary, by = "Tumor_Sample_Barcode") %>%
        mutate(target_gene_sequence = if_else(
          is.na(target_gene_sequence) | target_gene_sequence == "",
          variants,
          paste(target_gene_sequence, variants, sep = ",")
        )) %>%
        select(-variants)  # 不要な列は削除
      incProgress(1 / 8)

      Data_case_target = left_join(Data_case_target,
                                   gene_table,
                                   by = c("C.CAT調査結果.基本項目.ハッシュID" = "Tumor_Sample_Barcode"))
      Data_case_target_table = Data_case_target %>% dplyr::select(-症例.基本情報.がん種.OncoTree.)
      replacements <- c(
        "症例\\.検体情報\\.病理診断名" = "症例.検体情報.提出検体の病理診断名",
        "C.CAT調査結果\\.基本項目\\." = "",
        "症例\\.EP後レジメン情報\\." = "",
        "症例\\.基本情報\\." = "",
        "症例\\.検体情報\\." = "",
        "症例\\.背景情報\\." = "",
        "症例\\.がん種情報\\." = "",
        "症例\\.転帰情報\\." = "",
        "症例\\.管理情報\\." = "",
        "\\.名称\\." = "",
        "target_gene_sequence" = paste0(gene_name, "_sequence (set Gene to analyze)"),
        "ハッシュID" = "ID",
        "年齢" = "Age",
        "性別" = "Sex",
        "初回治療前のステージ分類" = "Stage at diagnosis",
        "提出検体の病理診断名" = "Diagnosis with CGP test specimen",
        "病理診断名" = "Pathological diagnosis",
        "臨床診断名" = "Clinical diagnosis",
        "がん種.OncoTree." = "OncoTree classification",
        "がん種.OncoTree.LEVEL1." = "OncoTree LEVEL 1",
        "パネル" = "Panel",
        "喫煙歴有無" = "Smoking history",
        "アルコール多飲有無" = "Heavy alcohol drinking",
        "重複がん有無.異なる臓器." = "Multiple cancer in multiple organ",
        "多発がん有無.同一臓器." = "Multiple tumor nodules in a organ",
        "登録時転移部位" = "Metastasis site at CGP test timing",
        "腫瘍細胞含有割合" = "Tumor content rate",
        "検体採取部位" = "Sampling site",
        "診断日" = "Diagnosis date",
        "登録日" = "Enrollment date",
        "転帰" = "Final survival status",
        "最終生存確認日" = "Last survival confirmation date",
        "死亡日" = "Date of death",
        "検体採取日.腫瘍組織." = "Sampling date of tumor tissue",
        "検査前使用薬剤" = "Drugs used before CGP test"
      )

      # 一括置換
      new_names <- str_replace_all(colnames(Data_case_target_table), replacements)
      colnames(Data_case_target_table) <- new_names
      OUTPUT_DATA$oncoprint_Data_case_target_table = Data_case_target_table
      incProgress(1 / 30)
    })
  })
  rm(analysis_env)
  gc()
}

output$download_mutation_data = downloadHandler(
  filename = "table_mutation_data.csv",
  content = function(file) {
    req(OUTPUT_DATA$oncoprint_Data_case_target_table)
    write_excel_csv(OUTPUT_DATA$oncoprint_Data_case_target_table, file)
  }
)
output$table_patient = DT::renderDataTable(server = TRUE,{
  req(OUTPUT_DATA$oncoprint_Data_case_target_table)
  create_datatable_with_confirm(OUTPUT_DATA$oncoprint_Data_case_target_table, "FELIS downloaded raw data in Table of clinical and mutation information per patient tab")
})

output$figure_lolliplot2 = renderGirafe({
  req(input$lolliplot_gene, input$lolliplot_no,
      OUTPUT_DATA$oncoprint_Data_MAF_target)
  withProgress(message = "Internet access may be necessary", {
    lolliplot_gene = input$lolliplot_gene
    if(lolliplot_gene == ""){
      lolliplot_gene = names(sort(table(OUTPUT_DATA$oncoprint_Data_MAF_target$Hugo_Symbol),decreasing = T))[[1]]
    }
    data_lolliplot = OUTPUT_DATA$oncoprint_Data_MAF_target %>%
      dplyr::select(Chromosome,
                    Start_Position,
                    End_Position,
                    Reference_Allele,
                    Tumor_Seq_Allele2,
                    Tumor_Sample_Barcode,
                    Hugo_Symbol,
                    Variant_Classification,
                    amino.acid.change) %>%
      dplyr::filter(!Variant_Classification %in% c("rearrangement",
                                                   "truncation",
                                                   "amplification",
                                                   "fusion",
                                                   "TMB_MSI_high",
                                                   "duplication",
                                                   "exon_skipping",
                                                   "deletion",
                                                   "inversion",
                                                   "loss")) %>%
      dplyr::filter(!str_detect(amino.acid.change, "c.")) %>%
      dplyr::filter(!str_detect(amino.acid.change, "_")) %>%
      dplyr::filter(Hugo_Symbol == lolliplot_gene)
    if(nrow(data_lolliplot) == 0){
      g = gg_empty()
    } else {
      var<-unique(data_lolliplot[!duplicated(data_lolliplot),])#remove duplicates
      var[order(var$Hugo_Symbol,gsub("([A-Z]+)([0-9]+)","\\2",var$amino.acid.change)),]#order
      var.freq<-plyr::count(var,c("Hugo_Symbol","amino.acid.change","Variant_Classification"))
      var.aanum<-unlist(lapply(regmatches(var.freq$amino.acid.change,gregexpr(pattern="*(\\d+)",var.freq$amino.acid.change)),function(x) x[[1]]))
      var.plot.data<-cbind(var.freq,as.numeric(as.character(var.aanum)))#hard way to remove factor
      colnames(var.plot.data)[5]<-"var.aanum"
      var.plot.data = var.plot.data[var.plot.data$var.aanum >=1 & is.finite(var.plot.data$var.aanum) & !is.na(var.plot.data$var.aanum),]

      var.plot<-isolate(var.plot.data)
      var.plot<-var.plot[var.plot$Hugo_Symbol==lolliplot_gene,]
      var.plot<-var.plot[order(var.plot$var.aanum),]
      var.plot$amino.acid.change<-as.character(var.plot$amino.acid.change)
      var.plot$tooltip_text <- paste0(
        "Amino acid: ", var.plot$amino.acid.change, "\n",
        "Frequency: ", var.plot$freq
      )
      incProgress(1 / 4)

      g<-ggplot(var.plot,aes(var.aanum,freq,color=Variant_Classification,tooltip = tooltip_text,label=amino.acid.change)) +
        geom_segment_interactive(aes(x=var.aanum,y=0,xend=var.aanum,yend=freq),
                                 color=ifelse(var.plot$freq>=input$lolliplot_no,"orange","grey50"),
                                 linewidth=ifelse(var.plot$freq>=input$lolliplot_no,1.3,0.7)) +
        geom_point_interactive(aes(color = Variant_Classification,
                                   tooltip = tooltip_text,
                                   data_id = amino.acid.change), size = 3) +
        geom_text_repel_interactive(nudge_y=0.2,
                                    color=ifelse(var.plot$freq>=input$lolliplot_no,"orange",NA),
                                    size=ifelse(var.plot$freq>=input$lolliplot_no,3,3),
                                    fontface="bold", max.overlaps=Inf) +
        theme_classic() +
        theme(plot.title=element_text(face="bold",size=(15),hjust=0),axis.text.x=element_text(angle=90,hjust=1)) +
        labs(y="Mutation frequency in our dataset",
             title=paste(nrow(data_lolliplot), " small scale variants of ",lolliplot_gene," in the dataset",sep="")) +
        scale_y_continuous(breaks=seq(0,max(var.plot$freq),max(floor(max(var.plot$freq)/100)*10,5)),name="Frequency")
      proteins_acc<-read.table(file("source/UniProt.txt",
                                    encoding='UTF-8-BOM'),header=T)
      if(!lolliplot_gene %in% c("NKX2_1", "NKX3_1","H1_2", "H3_4","H3_5","H3_3A","H3_3B")){
        lolliplot_gene_raw = str_split(lolliplot_gene, "_",simplify = T)[[1]]
      } else {
        lolliplot_gene_raw = lolliplot_gene
      }
      if(length(proteins_acc[proteins_acc$Hugo_Symbol==lolliplot_gene_raw,2]) == 1){
        proteins_acc<-proteins_acc[proteins_acc$Hugo_Symbol==lolliplot_gene_raw,2]
        proteins_acc_url<-gsub(" ","%2C",proteins_acc)
        baseurl<-"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
        URL<-paste0(baseurl,proteins_acc_url)
        incProgress(1 / 6)
        flag_lolliplot = 0
        if(file.exists(paste0("source/lolliplot/", lolliplot_gene_raw, ".rda"))){
          load(paste0("source/lolliplot/", lolliplot_gene_raw, ".rda"))
          flag_lolliplot = 1
        } else if(file.exists(paste0(tempdir(), "/", lolliplot_gene_raw, ".rda"))){
          load(paste0("source/lolliplot/", lolliplot_gene_raw, ".rda"))
          flag_lolliplot = 1
        } else {
          prots_feat<-try(GET(URL,accept_json()), silent = FALSE)
          if (class(prots_feat) != "try-error"){
            flag_lolliplot = 1
            save(prots_feat, file=paste0(tempdir(), "/", lolliplot_gene_raw, ".rda"))
          }
        }
        if (flag_lolliplot == 1){
          prots_feat_red<-(httr::content(prots_feat))
          incProgress(1 / 6)
          features_total_plot<-NULL
          for(i in 1:length(prots_feat_red)){
            features_temp<-drawProteins::extract_feat_acc(prots_feat_red[[i]])
            features_temp$order<-i
            features_total_plot<-rbind(features_total_plot,features_temp)
          }
          plot_start<-0
          plot_end<-max(features_total_plot$end)
          if(nrow(features_total_plot)>0){
            g <- g +
              geom_rect_interactive(data=features_total_plot[1,],
                                    mapping=aes(xmin=plot_start,xmax=plot_end,ymin=-0.08*max(var.plot$freq),ymax=-0.06*max(var.plot$freq)),
                                    colour="grey",
                                    fill="grey",
                                    inherit.aes=F)
          }
          if("CHAIN"%in%(unique(features_total_plot$type))){
            g <- g +
              geom_rect_interactive(data=features_total_plot[features_total_plot$type=="CHAIN",],
                                    mapping=aes(xmin=begin,xmax=end,ymin=-0.1*max(var.plot$freq),ymax=-0.04*max(var.plot$freq)),
                                    colour="grey",
                                    fill="grey",
                                    inherit.aes=F)
          }
          incProgress(1 / 6)
          if("DNA_BIND"%in%(unique(features_total_plot$type))) {
            features_total_plot[features_total_plot$type=="DNA_BIND","description"] <- features_total_plot[features_total_plot$type=="DNA_BIND","type"]
            g <- g +
              geom_rect_interactive(data=features_total_plot[features_total_plot$type=="DNA_BIND",],
                                    mapping=aes(xmin=begin,xmax=end,ymin=-0.12*max(var.plot$freq),ymax=-0.02*max(var.plot$freq),fill=description),
                                    inherit.aes=F)
          }
          if("DOMAIN"%in%(unique(features_total_plot$type))){
            g <- g +
              geom_rect_interactive(data=features_total_plot[features_total_plot$type=="DOMAIN",],
                                    mapping=aes(xmin=begin,xmax=end,ymin=-0.12*max(var.plot$freq),ymax=-0.02*max(var.plot$freq),fill=description),
                                    inherit.aes=F)
          }
          if("MOTIF"%in%(unique(features_total_plot$type))){
            g <- g +
              geom_rect_interactive(data=features_total_plot[features_total_plot$type=="MOTIF",],
                                    mapping=aes(xmin=begin,xmax=end,ymin=-0.12*max(var.plot$freq),ymax=-0.02*max(var.plot$freq),fill=description),
                                    inherit.aes=F)
          }
          if("REPEAT"%in%(unique(features_total_plot$type))){
            g <- g +
              geom_rect_interactive(data=features_total_plot[features_total_plot$type=="REPEAT",],
                                    mapping=aes(xmin=begin,xmax=end,ymin=-0.12*max(var.plot$freq),ymax=-0.02*max(var.plot$freq),fill=description),
                                    inherit.aes=F)
          }
          g <- g +
            scale_x_continuous(breaks=round(seq(plot_start,plot_end,by=50), digits = 0),name="Amino acid number")
        }
      }
      incProgress(1 / 6)
    }
  })
  ggiraph::girafe(ggobj = g,width_svg = 15,height_svg = 5, options = list(opts_sizing(rescale = FALSE)))
})


output$figure_oncoprint = renderPlot({
  req(OUTPUT_DATA$oncoprint_Data_oncoprint,
      OUTPUT_DATA$oncoprint_Data_MAF_target,
      oncoprint_trigger(),
      oncoprint_trigger_2())
  max_samples = isolate(input$oncoprint_max_samples)
  max_samples = ifelse(is.null(max_samples), 5000, max_samples)
  if (nrow(OUTPUT_DATA$oncoprint_Data_oncoprint) > max_samples){
    sample_indices = sample(nrow(OUTPUT_DATA$oncoprint_Data_oncoprint),max_samples)
    Data_oncoprint = OUTPUT_DATA$oncoprint_Data_oncoprint[sample_indices, ]
    ENTIRE = "subsampled"
  } else {
    Data_oncoprint = OUTPUT_DATA$oncoprint_Data_oncoprint
    ENTIRE = "whole"
  }

  Data_MAF_target_oncoprint = OUTPUT_DATA$oncoprint_Data_MAF_target %>%
    dplyr::filter(Tumor_Sample_Barcode %in%
                    Data_oncoprint$C.CAT調査結果.基本項目.ハッシュID &
                    Variant_Classification %in% c("rearrangement",
                                                  "small_scale_variant",
                                                  "nonsense_frameshift",
                                                  "truncation",
                                                  "amplification",
                                                  "fusion",
                                                  "TMB_MSI_high",
                                                  "duplication",
                                                  "exon_skipping",
                                                  "deletion",
                                                  "inversion",
                                                  "loss"))
  gene = isolate(get_ordered_genes())
  gene_no = length(gene)
  oncogenic_genes = gene
  oncoprint_option = ifelse(is.null(input$oncoprint_option), "Null", input$oncoprint_option)
  mid_age = input$mid_age
  withProgress(message = sample(nietzsche)[1], {
    mat_oncoprint = matrix("",
                           nrow = nrow(Data_oncoprint),
                           ncol = gene_no)
    rownames(mat_oncoprint) = Data_oncoprint$C.CAT調査結果.基本項目.ハッシュID
    colnames(mat_oncoprint) = gene
    Data_TMB = (Data_MAF_target_oncoprint %>%
                  dplyr::distinct(Tumor_Sample_Barcode, TMB))
    colnames(Data_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB")
    Data_oncoprint = left_join(Data_oncoprint,
                               Data_TMB,
                               by = "C.CAT調査結果.基本項目.ハッシュID")
    Data_TMB = Data_oncoprint$TMB
    Data_TMB[is.na(Data_TMB)] <- 0
    incProgress(1 / 4)
    Data_major_gene = Data_MAF_target_oncoprint %>%
      dplyr::filter(Hugo_Symbol %in% colnames(mat_oncoprint)) %>%
      dplyr::arrange(Tumor_Sample_Barcode)

    Data_major_gene_agg <- Data_major_gene %>%
      group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
      summarise(Variant_Classification = paste(unique(Variant_Classification), collapse = ";"), .groups = "drop")
    mat_oncoprint_update <- acast(Data_major_gene_agg, Tumor_Sample_Barcode ~ Hugo_Symbol, value.var = "Variant_Classification", fill = "")
    mat_oncoprint[rownames(mat_oncoprint_update), colnames(mat_oncoprint_update)] <- mat_oncoprint_update
    incProgress(1 / 4)
    rownames(mat_oncoprint) = Data_oncoprint$症例.基本情報.がん種.OncoTree.
    mat_oncoprint = t(mat_oncoprint)
    column_title = paste0("OncoPrint for frequent gene mutations, ", nrow(Data_oncoprint), " ", ENTIRE, " patients")
    tmp = mat_oncoprint
    tmp[tmp != ""] = 1
    tmp[tmp == ""] = 0
    tmp = matrix(as.numeric(tmp), ncol = ncol(tmp))
    tmp2 = gene_no + 1 - seq_len(gene_no)
    tmp = data.frame(tmp)
    colnames(tmp) = 1:ncol(tmp)
    sample_order = 1:ncol(tmp)
    colPal3 <- colorRampPalette(brewer.pal(11, "Spectral"))
    anno_diag = colnames(mat_oncoprint)
    anno_age = Data_oncoprint$YoungOld
    Data_oncoprint = Data_oncoprint %>%
      dplyr::mutate(anno_meta = case_when(
        Lymph_met == "Yes" ~ "(+)",
        Brain_met == "Yes" ~ "(+)",
        Lung_met == "Yes" ~ "(+)",
        Bone_met == "Yes" ~ "(+)",
        Liver_met == "Yes" ~ "(+)",
        Other_met == "Yes" ~ "(+)",
        TRUE ~ "(-)"))
    anno_meta = Data_oncoprint$anno_meta
    diag_list = names(sort(table(anno_diag),decreasing = TRUE))
    diag_list = factor(diag_list, levels = diag_list)
    diag_col = colPal3(length(diag_list))
    meta_list = unique(anno_meta)
    age_list = unique(anno_age)
    names(diag_col) = diag_list
    diag_list2 = sort(unique(anno_diag))
    diag_col2 = diag_col
    names(diag_col2) = diag_list2
    tmp3 = factor(mat_oncoprint)
    tmp3_ref = sort(mat_oncoprint, decreasing = F)[[1]]
    tmp3 = relevel(tmp3, ref = tmp3_ref)
    tmp3 = matrix(as.numeric(as.logical(as.numeric(tmp3) - 1)),
                  ncol = ncol(mat_oncoprint))
    tmp3 = data.frame(tmp3)
    colnames(tmp3) = 1:ncol(tmp3)
    ord <- do.call(order, c(as.data.frame(t(tmp3[rev(tmp2), ])), list(decreasing = TRUE)))
    # 並び替え適用
    tmp3 <- tmp3[, ord]
    # 各列の移動先を取得
    sample_order <- rank(ord, ties.method = "first")
    sample_order_simple = sample_order
    if(oncoprint_option != "Null"){
      if(oncoprint_option == "Sort by pathology"){
        sample_order_final = NULL
        for(l in diag_list){
          sample_order_final = c(sample_order_final,
                                 sample_order[sample_order %in% (1:length(sample_order))[(anno_diag == l)]])
        }
      } else if(oncoprint_option == "Sort by age"){
        sample_order_final_age = NULL
        for(l in age_list){
          sample_order_final_age = c(sample_order_final_age,
                                     sample_order[sample_order %in% (1:length(sample_order))[(anno_age == l)]])
        }
      } else if(oncoprint_option == "Sort by metastasis status at CGP"){
        sample_order_final_meta = NULL
        for(l in meta_list){
          sample_order_final_meta = c(sample_order_final_meta,
                                      sample_order[sample_order %in% (1:length(sample_order))[(anno_meta == l)]])
        }
      }
      incProgress(1 / 4)
      # 縦に並べる凡例: Mutation Types, Age, Metastasis
      lgd_mut <- Legend(
        at = names(color_mut),
        labels = names(color_mut),
        title = "Mutation Types",
        legend_gp = gpar(fill = unlist(color_mut)),
        direction = "vertical"
      )
      lgd_age <- Legend(
        at = c("Older", "Younger"),
        labels = c("Older", "Younger"),
        title = paste0("Age <= ", mid_age, ", or older"),
        legend_gp = gpar(fill = c("green", "blue")),
        direction = "vertical"
      )
      lgd_meta <- Legend(
        at = c("(+)", "(-)"),
        labels = c("Positive", "Negative"),
        title = "Metastasis",
        legend_gp = gpar(fill = c("purple1", "orange")),
        direction = "vertical"
      )
      if(length(diag_col2) > 100){
        wrap_every <- max(14, ceil(length(diag_col2)/7))
        lgd_diag_chunks <- split(names(diag_col2), ceiling(seq_along(names(diag_col2)) / wrap_every))
        lgd_diag_rows <- lapply(seq_along(lgd_diag_chunks), function(i) {
          chunk <- lgd_diag_chunks[[i]]
          if (i == 1) {
            # 1つ目にタイトルを追加
            Legend(
              at = chunk,                          # 凡例項目
              labels = chunk,                      # ラベル
              legend_gp = gpar(fill = diag_col2[chunk]),  # 色設定
              title = "Diagnosis"                  # タイトル
            )
          } else {
            # 2つ目以降はタイトルなし
            Legend(
              at = chunk,                          # 凡例項目
              labels = chunk,                      # ラベル
              legend_gp = gpar(fill = diag_col2[chunk]),   # 色設定
              title = "\u200B"
            )
          }
        })
        # 横方向で分割された凡例をパッキング
        final_lgd_diag_with_title <- do.call(packLegend, c(lgd_diag_rows, list(direction = "horizontal", gap = unit(3, "mm"))))
        # 縦方向に並べる凡例を統合
        vertical_legends <- packLegend(
          lgd_mut,
          lgd_age,
          lgd_meta,
          direction = "vertical",  # 縦に並べる
          gap = unit(3, "mm")
        )
        final_legends = list(vertical_legends, final_lgd_diag_with_title)
      } else {
        wrap_every <- max(10, ceil(length(diag_col2)/6))
        lgd_diag_chunks <- split(names(diag_col2), ceiling(seq_along(names(diag_col2)) / wrap_every))
        lgd_diag_rows <- lapply(seq_along(lgd_diag_chunks), function(i) {
          chunk <- lgd_diag_chunks[[i]]
          if (i == 1) {
            # 1つ目にタイトルを追加
            Legend(
              at = chunk,                          # 凡例項目
              labels = chunk,                      # ラベル
              legend_gp = gpar(fill = diag_col2[chunk]),  # 色設定
              title = "Diagnosis"                  # タイトル
            )
          } else {
            # 2つ目以降はタイトルなし
            Legend(
              at = chunk,                          # 凡例項目
              labels = chunk,                      # ラベル
              legend_gp = gpar(fill = diag_col2[chunk]),   # 色設定
              title = "\u200B"
            )
          }
        })
        # 横方向で分割された凡例をパッキング
        final_lgd_diag_with_title <- do.call(packLegend, c(lgd_diag_rows, list(direction = "horizontal", gap = unit(3, "mm"))))
        # 縦方向に並べる凡例を統合
        vertical_legends <- packLegend(
          lgd_age,
          lgd_meta,
          direction = "vertical",  # 縦に並べる
          gap = unit(3, "mm")
        )
        final_legends = list(lgd_mut, vertical_legends, final_lgd_diag_with_title)
      }
      ha = HeatmapAnnotation(
        Diagnosis = anno_diag,
        Age = anno_age,
        Metastasis = anno_meta,
        col = list(Diagnosis = diag_col2,
                   Age = c("Older" = "green", "Younger" = "blue"),
                   Metastasis = c("(+)" = "purple1", "(-)" = "orange")),
        annotation_height = unit(15, "mm"),
        show_legend = FALSE
      )
      hb = HeatmapAnnotation(TMB = anno_points(Data_TMB, size = unit(0.6, "mm"),
                                               ylim = c(0, 30),
                                               axis_param = list(
                                                 side = "right",
                                                 at = c(0, 10, 20, 30),
                                                 labels = c("0", "10", "20", ">30")
                                               )),
                             annotation_height = unit(15, "mm"),
                             show_legend = FALSE)
      if(oncoprint_option == "Sort by mutation frequency"){
        ht = oncoPrint(mat_oncoprint,
                       get_type = get_type_optimized,
                       column_order = sample_order_simple,
                       row_order = rev(tmp2),
                       alter_fun = alter_fun,
                       col = color_mut,
                       pct_gp = gpar(fontsize = 10),
                       column_title = column_title,
                       remove_empty_columns = FALSE,
                       remove_empty_rows = FALSE,
                       bottom_annotation = ha,
                       top_annotation = hb,
                       use_raster = FALSE,
                       raster_quality = 2,
                       alter_fun_is_vectorized = TRUE,
                       show_column_names = FALSE,
                       show_row_names = TRUE)
      } else if(oncoprint_option == "Sort by age"){
        ht = oncoPrint(mat_oncoprint,
                       get_type = get_type_optimized,
                       column_order = sample_order_final_age,
                       row_order = rev(tmp2),
                       alter_fun = alter_fun,
                       col = color_mut,
                       pct_gp = gpar(fontsize = 10),
                       column_title = column_title,
                       remove_empty_columns = FALSE,
                       remove_empty_rows = FALSE,
                       bottom_annotation = ha,
                       top_annotation = hb,
                       use_raster = FALSE,
                       raster_quality = 2,
                       alter_fun_is_vectorized = TRUE,
                       show_column_names = FALSE,
                       show_row_names = TRUE)
      } else if(oncoprint_option == "Sort by metastasis status at CGP"){
        ht = oncoPrint(mat_oncoprint,
                       get_type = get_type_optimized,
                       column_order = sample_order_final_meta,
                       row_order = rev(tmp2),
                       alter_fun = alter_fun,
                       col = color_mut,
                       pct_gp = gpar(fontsize = 10),
                       column_title = column_title,
                       remove_empty_columns = FALSE,
                       remove_empty_rows = FALSE,
                       bottom_annotation = ha,
                       top_annotation = hb,
                       use_raster = FALSE,
                       raster_quality = 2,
                       alter_fun_is_vectorized = TRUE,
                       show_column_names = FALSE,
                       show_row_names = TRUE)
      } else{
        ht = oncoPrint(mat_oncoprint,
                       get_type = get_type_optimized,
                       column_order = sample_order_final,
                       row_order = rev(tmp2),
                       alter_fun = alter_fun,
                       col = color_mut,
                       pct_gp = gpar(fontsize = 10),
                       column_title = column_title,
                       remove_empty_columns = FALSE,
                       remove_empty_rows = FALSE,
                       bottom_annotation = ha,
                       top_annotation = hb,
                       use_raster = FALSE,
                       raster_quality = 2,
                       alter_fun_is_vectorized = TRUE,
                       show_column_names = FALSE,
                       show_row_names = TRUE)
      }
      draw(
        ht,
        show_heatmap_legend = FALSE,
        annotation_legend_list = final_legends,  # 凡例を縦・横に配置
        annotation_legend_side = "bottom",  # 凡例を下部に配置
        heatmap_legend_side = "bottom"      # ヒートマップの凡例も下部に配置
      )
    }
    incProgress(1 / 4)
  })
})


# orderInputを動的に描画
output$dynamic_order <- renderUI({
  # 必要な入力値の存在確認
  req(input$oncoprint_gene_no, input$oncoprint_max_samples,
      OUTPUT_DATA$oncoprint_Data_MAF_target)

  gene_count <- as.numeric(input$oncoprint_gene_no)

  # 全体の腫瘍遺伝子リストを頻度順で取得
  oncogenic_genes <- names(sort(table(OUTPUT_DATA$oncoprint_Data_MAF_target$Hugo_Symbol), decreasing = TRUE))

  # ユーザー選択遺伝子の処理
  if (!is.null(input$oncoprint_gene) && length(input$oncoprint_gene) > 0) {
    valid_selected_genes <- input$oncoprint_gene[input$oncoprint_gene %in% OUTPUT_DATA$oncoprint_Data_MAF_target$Hugo_Symbol]
    if (length(valid_selected_genes) > 0) {
      oncogenic_genes <- unique(c(valid_selected_genes, oncogenic_genes))
    }
  }

  # 表示する遺伝子数を制限
  final_gene_count <- min(length(oncogenic_genes), gene_count)
  display_genes <- oncogenic_genes[1:final_gene_count]

  if (length(display_genes) == 0) {
    return(div(class = "alert alert-warning", icon("exclamation-triangle"), " no genes can be displayed."))
  }

  tagList(
    h4("Gene Display Order Settings"),
    p(paste0("Displayed Genes: ", final_gene_count, " / ", length(oncogenic_genes))),

    # カスタムCSS
    tags$style(HTML("
        .gene-list {
          min-height: 200px;
          border: 2px dashed #ccc;
          border-radius: 8px;
          padding: 15px;
          background-color: #fafafa;
        }
        .gene-item {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          border-radius: 20px;
          padding: 8px 16px;
          margin: 5px;
          display: inline-block;
          cursor: move;
          user-select: none;
          transition: all 0.3s ease;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .gene-item:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        .gene-item.dragging {
          opacity: 0.5;
          transform: rotate(5deg);
        }
      ")),

    # JavaScript for drag and drop
    tags$script(HTML("
        $(document).ready(function() {
          let draggedElement = null;
          let draggedIndex = null;

          function makeItemsDraggable() {
            $('.gene-item').attr('draggable', true);

            $('.gene-item').off('dragstart dragover drop dragend');

            $('.gene-item').on('dragstart', function(e) {
              draggedElement = this;
              draggedIndex = $(this).index();
              $(this).addClass('dragging');
              e.originalEvent.dataTransfer.effectAllowed = 'move';
            });

            $('.gene-item').on('dragover', function(e) {
              e.preventDefault();
              e.originalEvent.dataTransfer.dropEffect = 'move';
            });

            $('.gene-item').on('drop', function(e) {
              e.preventDefault();
              if (this !== draggedElement) {
                let targetIndex = $(this).index();
                let container = $(this).parent();

                if (draggedIndex < targetIndex) {
                  $(this).after(draggedElement);
                } else {
                  $(this).before(draggedElement);
                }

                updateGeneOrder();
              }
            });

            $('.gene-item').on('dragend', function(e) {
              $(this).removeClass('dragging');
              draggedElement = null;
              draggedIndex = null;
            });
          }

          function updateGeneOrder() {
            let geneOrder = [];
            $('.gene-item').each(function() {
              geneOrder.push($(this).text().trim());
            });
            Shiny.setInputValue('gene_order', geneOrder);
          }

          // 初期化
          setTimeout(makeItemsDraggable, 100);

          // Shinyの再描画後に再初期化
          $(document).on('shiny:value', function(event) {
            if (event.target.id === 'dynamic_order') {
              setTimeout(makeItemsDraggable, 100);
            }
          });
        });
      ")),

    div(
      class = "gene-list",
      id = "gene_container",
      lapply(display_genes, function(gene) {
        div(class = "gene-item", gene)
      })
    ),
    br(),
    div(
      class = "well well-sm",
      h5("Instructions:"),
      tags$ul(
        tags$li("Drag and drop gene names to change their order."),
        tags$li("The order will be reflected in the Oncoprint."),
        tags$li("Genes at the top will appear in higher-priority positions.")
      )
    )
  )
})

# オンコプリント生成用の遺伝子順序を取得する関数
get_ordered_genes <- reactive({
  req(input$oncoprint_gene_no, OUTPUT_DATA$oncoprint_Data_MAF_target)

  if (is.null(input$gene_order) || length(input$gene_order) == 0) {
    # デフォルトの順序（頻度順）を返す
    oncogenic_genes <- names(sort(table(OUTPUT_DATA$oncoprint_Data_MAF_target$Hugo_Symbol), decreasing = TRUE))
    gene_count <- min(input$oncoprint_gene_no, length(oncogenic_genes))
    return(oncogenic_genes[1:gene_count])
  }

  # ユーザー定義の順序を返す
  return(input$gene_order[input$gene_order %in% OUTPUT_DATA$oncoprint_Data_MAF_target$Hugo_Symbol])
})

color_mut = c("amplification" = "red",
              "loss" = "blue",
              "nonsense_frameshift" = "green4",
              "TMB_MSI_high" = "orange",
              "exon_skipping" = "red4",
              "small_scale_variant" = "green",
              "duplication" = "pink",
              "rearrangement" = "purple",
              "fusion" = "gray",
              "truncation" = "black",
              "inversion" = "coral2",
              "deletion" = "yellow")
height_ratios <- c(
  amplification = 1.0,
  loss = 0.95,
  nonsense_frameshift = 0.9,
  TMB_MSI_high = 0.85,
  exon_skipping = 0.8,
  small_scale_variant = 0.75,
  duplication = 0.7,
  rearrangement = 0.65,
  fusion = 0.6,
  truncation = 0.55,
  inversion = 0.45,
  deletion = 0.5
)
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#F4F4F4", col = NA))
  },
  amplification = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["amplification"],
              gp = gpar(fill = color_mut["amplification"], col = NA))
  },
  loss = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["loss"],
              gp = gpar(fill = color_mut["loss"], col = NA))
  },
  nonsense_frameshift = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["nonsense_frameshift"],
              gp = gpar(fill = color_mut["nonsense_frameshift"], col = NA))
  },
  TMB_MSI_high = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["TMB_MSI_high"],
              gp = gpar(fill = color_mut["TMB_MSI_high"], col = NA))
  },
  exon_skipping = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["exon_skipping"],
              gp = gpar(fill = color_mut["exon_skipping"], col = NA))
  },
  small_scale_variant = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["small_scale_variant"],
              gp = gpar(fill = color_mut["small_scale_variant"], col = NA))
  },
  duplication = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["duplication"],
              gp = gpar(fill = color_mut["duplication"], col = NA))
  },
  rearrangement = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["rearrangement"],
              gp = gpar(fill = color_mut["rearrangement"], col = NA))
  },
  fusion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["fusion"],
              gp = gpar(fill = color_mut["fusion"], col = NA))
  },
  truncation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["truncation"],
              gp = gpar(fill = color_mut["truncation"], col = NA))
  },
  inversion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["inversion"],
              gp = gpar(fill = color_mut["inversion"], col = NA))
  },
  deletion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "pt"), h * height_ratios["deletion"],
              gp = gpar(fill = color_mut["deletion"], col = NA))
  }
)
