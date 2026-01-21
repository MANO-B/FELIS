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
        if(ENV_ != "server")
          QS_SAVE(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), clin_tmp, file=file.path(tempdir(), "variant_data.qs"))
        return(clin_tmp)
      }
    })
  })
}
Data_report_raw = reactive({
  req(Data_report_raw_original())
  req(input$fusion)
  req(input$T_N)
  withProgress(message = "Mutation data loading.", {

    clin_tmp = Data_report_raw_original()

    # --- 1. Re-annotation (高速化: forループ廃止とleft_join利用) ---
    if(!is.null(input$reaanotation)){
      # ファイルを一括読み込み
      reaanotation_list <- lapply(input$reaanotation$datapath, function(x){
        read.csv(x, header = TRUE, encoding = 'UTF-8-BOM')
      }) %>% dplyr::bind_rows()

      # キーを作成して結合 (match関数によるループ処理を排除)
      clin_tmp <- clin_tmp %>%
        dplyr::mutate(Gene_alt = paste0(C.CAT調査結果.変異情報.マーカー, "_", C.CAT調査結果.変異情報.変異内容)) %>%
        dplyr::left_join(reaanotation_list %>% dplyr::select(Gene_alt, Level), by = "Gene_alt") %>%
        dplyr::mutate(
          C.CAT調査結果.エビデンス.エビデンスレベル = case_when(
            !is.na(Level) ~ Level,
            TRUE ~ C.CAT調査結果.エビデンス.エビデンスレベル
          )
        ) %>%
        dplyr::select(-Gene_alt, -Level)
    }

    # --- 2. 文字列置換 (高速化: str_replace_allで一括処理) ---
    # 置換リストの定義
    rep_patterns <- c(
      "NKX2\\-1" = "NKX2XXXX1", "NKX3\\-1" = "NKX3XXXX1",
      "H1\\-2" = "H1XXXX2", "H3\\-3A" = "H3XXXX3A", "H3\\-3B" = "H3XXXX3B",
      "H3\\-4" = "H3XXXX4", "H3\\-5" = "H3XXXX5",
      "HLA\\-A" = "HLAXXXXA", "HLA\\-DRB1" = "HLAXXXXDRB1",
      "-" = "_" # 最後にハイフンをアンダースコアに
    )

    clin_tmp$C.CAT調査結果.変異情報.マーカー <- stringr::str_replace_all(
      clin_tmp$C.CAT調査結果.変異情報.マーカー, rep_patterns
    )

    incProgress(1 / 4)

    # --- 3. Fusion Gene Handling (既存ロジックを維持しつつ整理) ---
    # Fusion（"_"を含むもの）の有無を確認
    has_fusion <- any(stringr::str_detect(clin_tmp$C.CAT調査結果.変異情報.マーカー, "_"))

    if(has_fusion) {
      if(input$fusion == "Split into two and treat as mutation in each gene"){
        # Split logic
        clin_base = clin_tmp %>% dplyr::filter(!str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_fusion = clin_tmp %>% dplyr::filter(str_detect(C.CAT調査結果.変異情報.マーカー, "_"))

        # Split marker strings efficiently
        splits <- stringr::str_split(clin_fusion$C.CAT調査結果.変異情報.マーカー, "_", simplify = TRUE)

        clin_f1 <- clin_fusion %>% dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー, C.CAT調査結果.変異情報.マーカー = splits[,1])
        clin_f2 <- clin_fusion %>% dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー, C.CAT調査結果.変異情報.マーカー = splits[,2])

        clin_tmp <- dplyr::bind_rows(clin_base, clin_f1, clin_f2)

      } else if(input$fusion == "Split into two and treat like EML-fusion, ALK-fusion"){
        clin_base = clin_tmp %>% dplyr::filter(!str_detect(C.CAT調査結果.変異情報.マーカー, "_"))
        clin_fusion = clin_tmp %>% dplyr::filter(str_detect(C.CAT調査結果.変異情報.マーカー, "_"))

        splits <- stringr::str_split(clin_fusion$C.CAT調査結果.変異情報.マーカー, "_", simplify = TRUE)

        clin_f1 <- clin_fusion %>% dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー, C.CAT調査結果.変異情報.マーカー = paste0(splits[,1], "_fusion"))
        clin_f2 <- clin_fusion %>% dplyr::mutate(C.CAT調査結果.変異情報.変異内容 = C.CAT調査結果.変異情報.マーカー, C.CAT調査結果.変異情報.マーカー = paste0(splits[,2], "_fusion"))

        clin_tmp <- dplyr::bind_rows(clin_base, clin_f1, clin_f2)

      } else if(input$fusion == "All fusion genes are named as 'fusion_genes'"){
        clin_tmp <- clin_tmp %>%
          dplyr::mutate(
            C.CAT調査結果.変異情報.変異内容 = ifelse(str_detect(C.CAT調査結果.変異情報.マーカー, "_"), C.CAT調査結果.変異情報.マーカー, C.CAT調査結果.変異情報.変異内容),
            C.CAT調査結果.変異情報.マーカー = ifelse(str_detect(C.CAT調査結果.変異情報.マーカー, "_"), "fusion_genes", C.CAT調査結果.変異情報.マーカー)
          )
      }
    }

    # --- 4. Amplification Handling ---
    if(input$amplification == "Treat as indipendent gene"){
      # mutateとcase_whenで一括処理（filter -> rbindよりも高速）
      clin_tmp <- clin_tmp %>%
        dplyr::mutate(C.CAT調査結果.変異情報.マーカー = case_when(
          C.CAT調査結果.変異情報.変異内容 == "amplification" ~ paste0(C.CAT調査結果.変異情報.マーカー, "_amp"),
          C.CAT調査結果.変異情報.変異内容 == "loss" ~ paste0(C.CAT調査結果.変異情報.マーカー, "_loss"),
          TRUE ~ C.CAT調査結果.変異情報.マーカー
        ))
    }

    # Restore special markers (vectorized replace)
    restore_patterns <- c(
      "NKX2XXXX1" = "NKX2_1", "NKX3XXXX1" = "NKX3_1", "H1XXXX2" = "H1_2",
      "H3XXXX4" = "H3_4", "H3XXXX5" = "H3_5", "H3XXXX3B" = "H3_3B",
      "HLAXXXXA" = "HLA_A", "HLAXXXXDRB1" = "HLA_DRB1", "H3XXXX3A" = "H3_3A"
    )
    clin_tmp$C.CAT調査結果.変異情報.マーカー <- stringr::str_replace_all(
      clin_tmp$C.CAT調査結果.変異情報.マーカー, restore_patterns
    )

    incProgress(1 / 4)

    # --- 5. Mutation Rename (高速化: forループ廃止とleft_join利用) ---
    if(!is.null(input$mutation_rename)){
      # ファイル一括読み込み
      mutation_rename_list <- lapply(input$mutation_rename$datapath, function(x){
        read.csv(x, header = TRUE, encoding = 'UTF-8-BOM')
      }) %>% dplyr::bind_rows()

      # ターゲットとなる遺伝子リスト
      target_genes <- unique(mutation_rename_list$Gene)

      # left_joinでリネーム情報を付与
      # 結合キー: Gene(マーカー) と Mutation(変異内容)
      clin_tmp <- clin_tmp %>%
        dplyr::left_join(mutation_rename_list,
                         by = c("C.CAT調査結果.変異情報.マーカー" = "Gene",
                                "C.CAT調査結果.変異情報.変異内容" = "Mutation")) %>%
        dplyr::mutate(
          C.CAT調査結果.変異情報.変異内容 = case_when(
            # リネーム定義が見つかればそれを適用
            !is.na(Rename) ~ Rename,
            # リネーム定義がないが、ターゲット遺伝子である場合は "Other" にする
            C.CAT調査結果.変異情報.マーカー %in% target_genes ~ "Other",
            # それ以外は元のまま
            TRUE ~ C.CAT調査結果.変異情報.変異内容
          )
        ) %>%
        dplyr::select(-Rename) # 結合した列を削除
    }

    Data_MAF = clin_tmp

    # --- 6. 列名変更と整形 ---
    # 物理位置の複製
    Data_MAF$C.CAT調査結果.変異情報.物理位置2 = Data_MAF$C.CAT調査結果.変異情報.物理位置

    # 列選択とリネーム (selectでリネームも同時に行うと効率的)
    Data_MAF = Data_MAF %>%
      dplyr::select(
        Hugo_Symbol = C.CAT調査結果.変異情報.マーカー,
        Chromosome = C.CAT調査結果.変異情報.クロモソーム番号,
        Start_Position = C.CAT調査結果.変異情報.物理位置,
        End_Position = C.CAT調査結果.変異情報.物理位置2, # 元コードに合わせて物理位置2を使用
        Reference_Allele = C.CAT調査結果.変異情報.リファレンス塩基,
        Tumor_Seq_Allele2 = C.CAT調査結果.変異情報.塩基変化,
        Variant_Classification = C.CAT調査結果.変異情報.変異種類,
        Evidence_level = C.CAT調査結果.エビデンス.エビデンスレベル,
        Drug = C.CAT調査結果.エビデンス.薬剤,
        Tumor_Sample_Barcode = C.CAT調査結果.基本項目.ハッシュID,
        VAF = C.CAT調査結果.変異情報.変異アレル頻度,
        amino.acid.change = C.CAT調査結果.変異情報.変異内容,
        TMB = C.CAT調査結果.変異情報TMB.TMB,
        Somatic_germline = C.CAT調査結果.変異情報.変異由来
      )

    # --- 7. Variant Type Determination ---
    # 事前に文字数を計算しておくと少し速い
    n_ref <- nchar(Data_MAF$Reference_Allele)
    n_alt <- nchar(Data_MAF$Tumor_Seq_Allele2)

    Data_MAF = Data_MAF %>% dplyr::mutate(
      Variant_Type = case_when(
        Hugo_Symbol == "MSI" ~ "MSI",
        Hugo_Symbol == "TMB" ~ "TMB",
        n_ref == 1 & n_alt == 1 ~ "SNP",
        n_ref == 2 & n_alt == 2 ~ "DNP",
        n_ref == 3 & n_alt == 3 ~ "TNP",
        n_ref == 1 & Reference_Allele != "-" & n_alt > 1 ~ "INS",
        n_ref > 1 & Tumor_Seq_Allele2 != "-" ~ Variant_Classification,
        TRUE ~ "Other"
      )
    )

    incProgress(1 / 4)

    # --- 8. TMB Handling (高速化: readr::parse_number利用) ---
    # 文字列操作を減らすため parse_number を使用（数値以外の文字を除去して変換）
    # 元コード: "X Muts/Mb" -> "X" -> 数値化
    # parse_numberは "10 Muts/Mb" -> 10 を直接返す
    Data_MAF$TMB <- suppressWarnings(readr::parse_number(as.character(Data_MAF$TMB)))
    Data_MAF$TMB[is.na(Data_MAF$TMB)] <- 0

    # TMBレポート行の作成
    Data_report_TMB = Data_MAF %>%
      dplyr::filter(Hugo_Symbol == "TMB") %>%
      dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE) %>%
      dplyr::arrange(Tumor_Sample_Barcode) %>%
      dplyr::mutate(
        Evidence_level = ifelse(TMB >= input$TMB_threshold, "F", ""),
        amino.acid.change = ifelse(TMB >= input$TMB_threshold, "high", "")
      )

    # TMB以外の行を整形して結合
    Data_MAF = Data_MAF %>%
      dplyr::filter(Hugo_Symbol != "TMB") %>%
      dplyr::mutate(
        Evidence_level = case_when(
          Hugo_Symbol == "MSI" & amino.acid.change == "high" & Evidence_level == "" ~ "F",
          TRUE ~ Evidence_level
        )
      ) %>%
      dplyr::bind_rows(Data_report_TMB) %>%
      dplyr::arrange(Tumor_Sample_Barcode)

    # NAの処理
    Data_MAF$Evidence_level[is.na(Data_MAF$Evidence_level)] <- ""

    # --- 9. T/N Filtering ---
    if(!is.null(input$T_N)){
      if(input$T_N == "Only somatic for T/N panel"){
        Data_MAF = Data_MAF %>% dplyr::filter(Somatic_germline != "Germline")
      } else if(input$T_N == "Only germline"){
        Data_MAF = Data_MAF %>% dplyr::filter(Somatic_germline == "Germline")
      }
    }

    incProgress(1 / 4)
  })

  return(Data_MAF)
})

Data_report = reactive({
  #req(input$special_gene_annotation)
  # 生データの取得
  Data_MAF = Data_report_raw()

  # --- 1. Evidence Level の書き換え処理 ---
  if(!is.null(input$special_gene_annotation)){
    if(input$special_gene_annotation == "This gene is also analyzed for VUS and Benign" &
       !is.null(input$special_gene)){
      Data_MAF = Data_MAF %>% dplyr::mutate(
        Evidence_level = case_when(
          Hugo_Symbol == input$special_gene & Evidence_level == "" ~ "F",
          TRUE ~ Evidence_level
        )
      )
    }
  }

  # --- 2. 遺伝子シンボルの書き換え処理 ---
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

  # --- 3. 遺伝子グループ解析（フィルタリング処理） ---
  if(!is.null(input$gene_group_analysis) && !is.null(input$gene)){

    # 共通処理: ID_geneset の特定
    # 条件に応じて対象となるTumor_Sample_Barcodeを抽出
    target_maf <- Data_MAF %>% dplyr::filter(Hugo_Symbol %in% input$gene)

    if(input$patho == "Only pathogenic muts"){
      target_maf <- target_maf %>% dplyr::filter(Evidence_level == "F")
    }
    ID_geneset <- unique(target_maf$Tumor_Sample_Barcode)


    # A. 遺伝子セットに変異がある症例のみ残す場合
    if(input$gene_group_analysis == "Only cases with mutations in the gene set are analyzed"){
      Data_MAF = Data_MAF %>%
        dplyr::filter(Tumor_Sample_Barcode %in% ID_geneset)

      # B. 遺伝子セットに変異がある症例を除外する場合
    } else if(input$gene_group_analysis == "Only cases without mutations in the gene set are analyzed"){

      # まず除外を実行
      Data_MAF = Data_MAF %>%
        dplyr::filter(!Tumor_Sample_Barcode %in% ID_geneset)

      # 本来存在するはずの全症例IDを取得（除外対象以外）
      ALL_ID_vector = (Data_case_raw() %>% dplyr::filter(
        !C.CAT調査結果.基本項目.ハッシュID %in% ID_geneset
      ))$C.CAT調査結果.基本項目.ハッシュID

      # Data_MAFに残っていないID（変異が全くないため消えてしまった症例など）を特定
      # setdiff(A, B) は Aに含まれるがBに含まれない要素を返します
      DEFF_ID = setdiff(ALL_ID_vector, Data_MAF$Tumor_Sample_Barcode)

      # ★★★ 高速化の主要箇所 ★★★
      if(length(DEFF_ID) > 0){
        # テンプレートとなる行を抽出・作成
        TEST_base = Data_MAF %>%
          dplyr::filter(Hugo_Symbol == "TMB") %>%
          dplyr::slice(1) %>%
          dplyr::mutate(
            TMB = 0,
            Evidence_level = "",
            Drug = ""
          )

        # テンプレートが見つかった場合のみ実行
        if(nrow(TEST_base) > 0){
          # 1. テンプレートを行数分コピー (repを使ったインデックス参照で高速複製)
          TEST_ALL = TEST_base[rep(1, length(DEFF_ID)), ]

          # 2. IDを一括代入
          TEST_ALL$Tumor_Sample_Barcode = DEFF_ID

          # 3. 結合
          Data_MAF = dplyr::bind_rows(Data_MAF, TEST_ALL)
        }
      }
    }
  }

  # --- 最終整形 ---
  Data_MAF = Data_MAF %>%
    dplyr::distinct() %>%
    dplyr::arrange(Tumor_Sample_Barcode) %>%
    dplyr::select(-Somatic_germline)

  return(Data_MAF)
})
