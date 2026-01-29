files_to_remove <- c(
  "Data_mutation_cord.qs",
  "m2.qs", "lgb_m2.qs", "rf_m2.qs",
  "m2_rec.qs", "lgb_m2_rec.qs", "rf_m2_rec.qs",
  "m2_GMT.qs", "lgb_m2_GMT.qs", "rf_m2_GMT.qs",
  "PS_data.rda", "Organ_data.rda", "Lines_data.rda",
  "PS_distribution.pdf", "check_IPCW.pdf", "IPCW_weight.pdf", "IPW_weight.pdf",
  "love_plot_PSM.pdf"
)

paths <- file.path(tempdir(), files_to_remove)

# 存在するものだけ削除
file.remove(paths[file.exists(paths)])

if(Sys.getenv("SHINY_SERVER_VERSION") != ""){
  if(file.exists(file.path(tempdir(), "Data_mutation_cord.qs"))){
    file.remove(file.path(tempdir(), "Data_mutation_cord.qs"))
  }
  if(file.exists(file.path(tempdir(), "clinical_data.qs"))){
    file.remove(file.path(tempdir(), "clinical_data.qs"))
  }
  if(file.exists(file.path(tempdir(), "variant_data.qs"))){
    file.remove(file.path(tempdir(), "variant_data.qs"))
  }
  if(file.exists(file.path(tempdir(), "drug_data.qs"))){
    file.remove(file.path(tempdir(), "drug_data.qs"))
  }
}

if(AMAZON_FLAG){
  # 結果を見やすく出力
  Identity_text_1 = paste0(
    "System environment for AWS, etc.\n",
    "------------------\n",
    "AWS_ACCESS_KEY_ID:", Sys.getenv("AWS_ACCESS_KEY_ID"), "\n",
    "AWS_SESSION_TOKEN:", Sys.getenv("AWS_SESSION_TOKEN"), "\n",
    "AWS_ROLE_ARN:", Sys.getenv("AWS_ROLE_ARN"), "\n",
    "AWS_DEFAULT_REGION:", Sys.getenv("AWS_DEFAULT_REGION"), "\n",
    "AWS_REGION:", Sys.getenv("AWS_REGION"), "\n",
    "AWS_WEB_IDENTITY_TOKEN_FILE:", Sys.getenv("AWS_WEB_IDENTITY_TOKEN_FILE"), "\n",
    "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI:", Sys.getenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"), "\n",
    "AWS_CONTAINER_CREDENTIALS_FULL_URI:", Sys.getenv("AWS_CONTAINER_CREDENTIALS_FULL_URI"), "\n",
    "AWS_PROFILE:", Sys.getenv("AWS_PROFILE"), "\n\n",
    "# STSサービスクライアントを作成\n",
    "sts <- paws::sts(config = list(region = 'ap-northeast-1'))\n",
    "auth_info <- tryCatch({\n",
    "# GetCallerIdentityを実行\n",
    "identity <- sts$get_caller_identity()\n",
    "info_text <- paste0(...略...)\n",
    "return(info_text)\n",
    "}, error = function(e) {\n",
    "return(paste0('認証エラー: AWSに接続できませんでした。\n詳細: ', e$message))\n",
    "})\n\n",
    "========= 実行結果 ========="
  )

  # STSサービスクライアントを作成
  sts <- paws::sts(config = list(region = "ap-northeast-1"))

  auth_info <- tryCatch({
    # GetCallerIdentityを実行
    identity <- sts$get_caller_identity()

    # 文字列を結合して返す (catではなくpasteやsprintfを使う)
    info_text <- paste0(
      "--- AWS Identity Info ---\n",
      "Account: ", identity$Account, "\n",
      "Arn:     ", identity$Arn, "\n",
      "UserId:  ", identity$UserId, "\n"
    )

    return(info_text)

  }, error = function(e) {
    # エラー時も文字列として返す
    return(paste0("認証エラー: AWSに接続できませんでした。\n詳細: ", e$message))
  })

  Identity_text_1 = paste0(Identity_text_1, "\n", auth_info)

  # ログ保存用変数の初期化
  EXCLUDED_TEXT_1 = ""

  download_withdrawn_list_safe <- function(bucket, key,
                                           local_path = file.path(tempdir(), "ic_withdrawn_list.txt")) {

    if(file.exists(local_path)){
      file.remove(local_path)
    }

    # ログバッファの準備
    log_buffer <- character(0)

    # ログ記録用関数
    add_log <- function(...) {
      msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
      log_buffer <<- c(log_buffer, msg)
    }

    add_log("=== S3ダウンロード開始 (東京リージョン固定) ===")
    add_log("Target Bucket: ", bucket)
    add_log("Target Key:    ", key)
    add_log("Target Region: ap-northeast-1")

    success <- FALSE

    connect_timeout_sec = 5
    total_timeout_sec = 10

    tryCatch({
      # タイムアウト設定
      httr_opts <- list(
        httr::config(connecttimeout = connect_timeout_sec),
        httr::timeout(total_timeout_sec)
      )

      # ダウンロード実行
      aws.s3::save_object(
        object = key,
        bucket = bucket,
        file   = local_path,
        region = "ap-northeast-1",
        check_region = FALSE,
        config = httr_opts,
        overwrite = TRUE
      )

      success <- TRUE
      add_log("  -> ダウンロード成功！ (Success)")

      if (file.exists(local_path)) {
        f_size <- file.info(local_path)$size
        add_log("  -> 保存先: ", local_path)
        add_log("  -> サイズ: ", f_size, " bytes")
      }

    }, error = function(e) {
      err_msg <- conditionMessage(e)
      add_log("  -> 失敗 (Error): ", err_msg)

      # --- 追加: 詳細なオブジェクト構造をログに取り込む ---
      add_log("  -> [詳細エラーオブジェクト構造]")
      # str(e) の出力を文字列ベクトルとしてキャプチャ
      detailed_structure <- capture.output(str(e))
      for(line in detailed_structure) {
        add_log("     ", line)
      }
      # -----------------------------------------------

      # エラー内容に応じたヒント
      if (grepl("403", err_msg)) {
        add_log("     [ヒント] 403 Forbidden: 場所は合っていますが、アクセス権限がありません。")
        add_log("     AWSアクセスキーまたはIAMロールの権限を確認してください。")
      } else if (grepl("404", err_msg)) {
        add_log("     [ヒント] 404 Not Found: バケットまたはファイル名が間違っています。")
      }
    })

    add_log("=== 処理終了 ===")

    # ログをグローバル変数に保存
    EXCLUDED_TEXT_1 <<- paste(log_buffer, collapse = "\n")

    # 失敗時は空ファイル作成
    if (!success && !file.exists(local_path)) {
      file.create(local_path)
      add_log("  -> (安全策) 空ファイルを作成しました")
    }

    return(local_path)
  }

  # 実行
  withdrawn_file <- download_withdrawn_list_safe(bucket = bucket, key = key)

  ID_exclude <- tryCatch({
    # ファイルサイズ0または存在しない場合は空文字を返す
    if (!file.exists(withdrawn_file) || is.na(file.info(withdrawn_file)$size) || file.info(withdrawn_file)$size == 0) {
      character(0)
    } else {
      # エンコーディング判定
      Encode <- tryCatch(
        readr::guess_encoding(withdrawn_file, n_max = 1000)[[1]]$encoding,
        error = function(e) "UTF-8"
      )
      # readr::guess_encoding は Shift_JIS を返すことがあるが、locale では CP932 推奨
      if (!is.null(Encode) && identical(Encode, "Shift_JIS")) Encode <- "CP932"
      if (is.null(Encode)) Encode <- "UTF-8" # 念のため

      # CSV読み込み (headerなし前提)
      df <- readr::read_csv(
        file = withdrawn_file,
        col_names = "X1", # ヘッダーなしとして読み込み、列名をX1に固定
        locale = readr::locale(encoding = Encode),
        progress = FALSE,
        show_col_types = FALSE
      )

      # データフレームが空でないか確認してベクトル化
      if(nrow(df) > 0){
        if(stringr::str_detect(df$X1[1], "xml version")){
          character(0)
        } else if(nrow(df) > 1){
          if(stringr::str_detect(df$X1[2], "<Error>")) character(0)
        } else {
          as.character(df$X1)
        }
      }else character(0)
    }
  }, error = function(e) {
    # 読み込みエラー時は安全に空ベクトルを返す
    if(interactive()) message("Error reading csv: ", e$message)
    character(0)
  })

  # 結果確認
  EXCLUDED_TEXT_2 = (paste("除外リスト:", withdrawn_file, ", 件数:", length(ID_exclude)))

  capture_s3_get_probe <- function(bucket, key, region = "ap-northeast-1") {
    tmp <- tempfile(fileext = ".txt")

    txt_lines <- capture.output({
      cat("=== S3 GetObject probe (aws.s3::save_object) ===\n")
      cat("time_utc: ", format(Sys.time(), tz = "UTC"), "\n", sep = "")
      cat("bucket  : ", bucket, "\n", sep = "")
      cat("key     : ", key, "\n", sep = "")
      cat("region  : ", region, "\n", sep = "")
      cat("tmpfile : ", tmp, "\n", sep = "")
      cat("\n--- code ---\n")
      cat('tmp <- tempfile(fileext = ".txt")\n')
      cat('res_get <- tryCatch({\n')
      cat('  aws.s3::save_object(object = key, bucket = bucket, file = tmp,\n')
      cat('                      region = "ap-northeast-1", check_region = FALSE, overwrite = TRUE)\n')
      cat('  paste("OK: GetObject/save_object allowed. downloaded to", tmp)\n')
      cat('}, error = function(e) paste("NG: GetObject/save_object failed:", conditionMessage(e)))\n')
      cat("\n--- output ---\n")

      # 実行
      res_get <- tryCatch({
        aws.s3::save_object(
          object = key, bucket = bucket, file = tmp,
          region = region, check_region = FALSE, overwrite = TRUE
        )
        paste("OK: GetObject/save_object allowed. downloaded to", tmp)
      }, error = function(e) {
        paste("NG: GetObject/save_object failed:", conditionMessage(e))
      })

      # cat の部分（あなたがやったのと同じ）
      cat(res_get, "\n")
      cat("aws.s3::save_objectを実行するたびにHostIdが変化するらしい\n")

      # 追加：aws.s3 の aws_error オブジェクトも取れるように（tryCatchで“捕捉”）
      cat("\n--- structured (tryCatch error object) ---\n")
      err_obj <- tryCatch({
        aws.s3::save_object(
          object = key, bucket = bucket, file = tmp,
          region = region, check_region = FALSE, overwrite = TRUE
        )
        NULL
      }, error = function(e) e)

      if (is.null(err_obj)) {
        cat("No error object (success)\n")
      } else {
        # print結果をそのまま
        print(err_obj)
      }
    })

    # 1本の文字列にまとめて返す（監視用途向け）
    paste(txt_lines, collapse = "\n")
  }

  # 使い方：
  log_text_S3_1 <- capture_s3_get_probe(bucket = bucket, key = key)

  log_text_S3_2 <- paste(
    "aws.signature::locate_credentials()",
    "---------------",
    capture.output({
      creds <- aws.signature::locate_credentials()

      print(list(
        key     = creds$key,
        secret  = "********",
        region  = creds$region,
        file    = creds$file,
        profile = creds$profile
      ))
    }), collapse = "\n")
} else {
  ID_exclude = character(0)
}
