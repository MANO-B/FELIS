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
  files_to_remove <- c(
    "clinical_data.qs", "variant_data.qs", "drug_data.qs")
  file.remove(paths[file.exists(paths)])
}


# ==============================================================================
# AWS ECS Credential Retrieval and S3 Download Script
# ==============================================================================

bucket <- "3e040086-data-production-felisanalysisportal"
key    <- "ic_withdrawn_list.txt"
connect_timeout_sec = 5
total_timeout_sec = 10
AMAZON_FLAG <- file.exists("/srv/shiny-server/felis-ccat/AMAZON")
# AMAZON_FLAG <- TRUE # for debugging

if(AMAZON_FLAG){
  library(jsonlite)
  library(httr)
  library(aws.s3)
  library(readr)
  library(stringr)

  # ----------------------------------------------------------------------------
  # 1. Fetch Credentials from ECS Task Metadata Endpoint
  # ----------------------------------------------------------------------------

  # Get the relative URI from the environment variable
  ecs_creds_uri <- Sys.getenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI")
  Identity_text_1 <- "【Log: Credential Retrieval】"

  # Check if the URI is defined (Fallbacks for local testing if needed)
  if (ecs_creds_uri == "") {
    # If empty, we can't fetch credentials. Log this state.
    # Note: In a real ECS environment, this should not be empty.
    Identity_text_1 <- paste0(Identity_text_1, "\nAWS_CONTAINER_CREDENTIALS_RELATIVE_URI is not defined.")
  } else {
    Identity_text_1 <- paste0(Identity_text_1, "\nAWS_CONTAINER_CREDENTIALS_RELATIVE_URI: http://169.254.170.2", ecs_creds_uri)
  }

  # Construct the full URL
  # The IP 169.254.170.2 is static for ECS Task Metadata
  credential_url <- paste0("http://169.254.170.2", ecs_creds_uri)

  # Execute HTTP GET request with a timeout
  response <- tryCatch({
    httr::GET(credential_url, timeout(5))
  }, error = function(e) {
    return(NULL) # Return NULL on connection failure
  })

  # Initialize credential list
  creds_list <- NULL

  # Process the response
  if (is.null(response)) {
    Identity_text_1 <- paste0(Identity_text_1, "\nError: Failed to connect to credential endpoint (Connection Error).")
  } else if (httr::status_code(response) != 200) {
    Identity_text_1 <- paste0(Identity_text_1, "\nError: HTTP Status ", as.character(httr::status_code(response)))
  } else {
    Identity_text_1 <- paste0(Identity_text_1, "\nSuccess: Connected to credential endpoint. Status 200.")

    # Parse JSON content
    content_text <- httr::content(response, as = "text", encoding = "UTF-8")

    creds_list <- tryCatch({
      parsed <- jsonlite::fromJSON(content_text)

      # Verify required keys exist
      required_keys <- c("AccessKeyId", "SecretAccessKey", "Token")
      if (!all(required_keys %in% names(parsed))) {
        stop("JSON missing required credential keys.")
      }

      Identity_text_1 <<- paste0(Identity_text_1, "\nJSON validation: OK")
      parsed
    }, error = function(e) {
      Identity_text_1 <<- paste0(Identity_text_1, "\nWarning: JSON parsing failed. Using dummy credentials.")
      return(NULL)
    })

    if(!is.null(creds_list)) {
      Identity_text_1 <- paste0(Identity_text_1, "\nExtracted Token Expiration: ", creds_list$Expiration)
    }
  }

  # ----------------------------------------------------------------------------
  # 2. Set Credentials to Variables (Handling Failures)
  # ----------------------------------------------------------------------------

  if(is.null(creds_list)){
    # Fallback to dummy credentials if retrieval failed
    # (Operations will fail, but the app won't crash immediately)
    my_access_key    <- "DUMMY_ACCESS_KEY"
    my_secret_key    <- "DUMMY_SECRET_KEY"
    my_session_token <- "DUMMY_TOKEN"
    my_expiration    <- "N/A"
  } else {
    my_access_key    <- creds_list$AccessKeyId
    my_secret_key    <- creds_list$SecretAccessKey
    my_session_token <- creds_list$Token
    my_expiration    <- creds_list$Expiration
  }

  message("AccessKeyId retrieved: ", my_access_key)

  # ============================================================================
  # [IMPORTANT] Store keys in global variables for explicit passing
  # ============================================================================
  MY_AWS_ACCESS_KEY_ID     <- my_access_key
  MY_AWS_SECRET_ACCESS_KEY <- my_secret_key
  MY_AWS_SESSION_TOKEN     <- my_session_token
  MY_AWS_REGION            <- "ap-northeast-1"

  # Also set as environment variables for compatibility with other libs/logging
  Sys.setenv("AWS_ACCESS_KEY_ID"     = MY_AWS_ACCESS_KEY_ID)
  Sys.setenv("AWS_SECRET_ACCESS_KEY" = MY_AWS_SECRET_ACCESS_KEY)
  Sys.setenv("AWS_SESSION_TOKEN"     = MY_AWS_SESSION_TOKEN)
  Sys.setenv("AWS_DEFAULT_REGION"    = MY_AWS_REGION)
  Sys.setenv("AWS_REGION"            = MY_AWS_REGION)

  # ----------------------------------------------------------------------------
  # 3. Log Environment Details
  # ----------------------------------------------------------------------------

  Identity_text_1 <- paste0(Identity_text_1, "\n\n",
                            "System environment for AWS\n",
                            "------------------\n",
                            "AWS_ACCESS_KEY_ID: ", Sys.getenv("AWS_ACCESS_KEY_ID"), "\n",
                            "AWS_SESSION_TOKEN: ", Sys.getenv("AWS_SESSION_TOKEN"), "\n", # Token is visible here
                            "AWS_DEFAULT_REGION: ", Sys.getenv("AWS_DEFAULT_REGION"), "\n",
                            "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI: ", Sys.getenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"), "\n\n",
                            "========= STS Caller Identity Check =========\n"
  )

  # Check identity using the explicitly retrieved credentials
  paws_config <- list(
    region = MY_AWS_REGION,
    credentials = list(
      creds = list(
        access_key_id = MY_AWS_ACCESS_KEY_ID,
        secret_access_key = MY_AWS_SECRET_ACCESS_KEY,
        session_token = if(MY_AWS_SESSION_TOKEN != "") MY_AWS_SESSION_TOKEN else NULL
      )
    )
  )

  auth_info <- tryCatch({
    sts <- paws::sts(config = paws_config)
    identity <- sts$get_caller_identity()

    paste0(
      "--- AWS Identity Info ---\n",
      "Account: ", identity$Account, "\n",
      "Arn:     ", identity$Arn, "\n",
      "UserId:  ", identity$UserId, "\n"
    )
  }, error = function(e) {
    paste0("Authentication Error: Could not connect to AWS STS.\nDetail: ", e$message)
  })

  Identity_text_1 <- paste0(Identity_text_1, auth_info)

  # ----------------------------------------------------------------------------
  # 4. Define Download Function with Explicit Credentials
  # ----------------------------------------------------------------------------

  EXCLUDED_TEXT_1 <- "" # Initialize log variable

  download_withdrawn_list_safe <- function(bucket, key,
                                           local_path = file.path(tempdir(), "ic_withdrawn_list.txt")) {

    if(file.exists(local_path)){
      file.remove(local_path)
    }

    log_buffer <- character(0)
    add_log <- function(...) {
      msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
      log_buffer <<- c(log_buffer, msg)
    }

    add_log("=== S3 Download Start (Explicit Credentials) ===")
    add_log("Target Bucket: ", bucket)
    add_log("Target Key:    ", key)
    add_log("Target Region: ", MY_AWS_REGION)

    success <- FALSE
    connect_timeout_sec <- 5
    total_timeout_sec <- 10

    tryCatch({
      httr_opts <- list(
        httr::config(connecttimeout = connect_timeout_sec),
        httr::timeout(total_timeout_sec)
      )

      # [CRITICAL] Passing credentials explicitly here
      aws.s3::save_object(
        object = key,
        bucket = bucket,
        file   = local_path,
        region = MY_AWS_REGION,
        key    = MY_AWS_ACCESS_KEY_ID,      # Explicit
        secret = MY_AWS_SECRET_ACCESS_KEY,  # Explicit
        session_token = if(MY_AWS_SESSION_TOKEN != "") MY_AWS_SESSION_TOKEN else NULL, # Explicit
        check_region = FALSE,
        config = httr_opts,
        overwrite = TRUE
      )

      success <- TRUE
      add_log(" -> Download Success!")

      if (file.exists(local_path)) {
        f_size <- file.info(local_path)$size
        add_log(" -> Path: ", local_path)
        add_log(" -> Size: ", f_size, " bytes")
      }

    }, error = function(e) {
      err_msg <- conditionMessage(e)
      add_log(" -> Failed (Error): ", err_msg)

      # Log detailed error structure
      add_log(" -> [Error Details]")
      detailed_structure <- capture.output(str(e))
      for(line in detailed_structure) {
        add_log("    ", line)
      }

      if (grepl("403", err_msg)) {
        add_log("    [Hint] 403 Forbidden: Check IAM Role permissions.")
      } else if (grepl("404", err_msg)) {
        add_log("    [Hint] 404 Not Found: Check Bucket/Key name.")
      }
    })

    add_log("=== Process End ===")
    EXCLUDED_TEXT_1 <<- paste(log_buffer, collapse = "\n")

    # Safe fallback: create empty file if failed
    if (!success && !file.exists(local_path)) {
      file.create(local_path)
      add_log(" -> Created empty file as fallback.")
    }

    return(local_path)
  }

  # ----------------------------------------------------------------------------
  # 5. Execute Download and Read File
  # ----------------------------------------------------------------------------

  withdrawn_file <- download_withdrawn_list_safe(bucket = bucket, key = key)

  ID_exclude <- tryCatch({
    if (!file.exists(withdrawn_file) || is.na(file.info(withdrawn_file)$size) || file.info(withdrawn_file)$size == 0) {
      character(0)
    } else {
      # Encoding detection
      Encode <- tryCatch(
        readr::guess_encoding(withdrawn_file, n_max = 1000)[[1]]$encoding,
        error = function(e) "UTF-8"
      )
      if (!is.null(Encode) && identical(Encode, "Shift_JIS")) Encode <- "CP932"
      if (is.null(Encode)) Encode <- "UTF-8"

      # Read CSV (assuming no header)
      df <- readr::read_csv(
        file = withdrawn_file,
        col_names = "X1",
        locale = readr::locale(encoding = Encode),
        progress = FALSE,
        show_col_types = FALSE
      )

      # Validation and Vectorization
      if(nrow(df) > 0){
        if(stringr::str_detect(df$X1[1], "xml version")){
          character(0)
        } else if(nrow(df) > 1 && stringr::str_detect(df$X1[2], "<Error>")){
          character(0)
        } else {
          as.character(df$X1)
        }
      } else {
        character(0)
      }
    }
  }, error = function(e) {
    if(interactive()) message("Error reading csv: ", e$message)
    character(0)
  })

  EXCLUDED_TEXT_2 <- paste("Withdrawn List:", withdrawn_file, ", Count:", length(ID_exclude))

  # ----------------------------------------------------------------------------
  # 6. Capture Probe Function (Debugging Tool)
  # ----------------------------------------------------------------------------

  # Define timeouts globally for the probe function
  httr_opts_probe <- list(
    httr::config(connecttimeout = 5),
    httr::timeout(10)
  )

  capture_s3_get_probe <- function(bucket, key, region = "ap-northeast-1") {
    local_path_probe <- file.path(tempdir(), "ic_withdrawn_list_probe.txt")

    txt_lines <- capture.output({
      cat("=== S3 GetObject Probe (aws.s3::save_object) ===\n")
      cat("Time (UTC) : ", format(Sys.time(), tz = "UTC"), "\n", sep = "")
      cat("Bucket     : ", bucket, "\n", sep = "")
      cat("Key        : ", key, "\n", sep = "")
      cat("Region     : ", region, "\n", sep = "")
      cat("Local Path : ", local_path_probe, "\n", sep = "")
      cat("\n--- Execution Log ---\n")

      # Explicitly use global credential variables
      res_get <- tryCatch({
        aws.s3::save_object(
          object = key,
          bucket = bucket,
          file = local_path_probe,
          region = MY_AWS_REGION,
          key = MY_AWS_ACCESS_KEY_ID,
          secret = MY_AWS_SECRET_ACCESS_KEY,
          session_token = if(MY_AWS_SESSION_TOKEN != "") MY_AWS_SESSION_TOKEN else NULL,
          check_region = FALSE,
          config = httr_opts_probe,
          overwrite = TRUE
        )
        paste("OK: GetObject/save_object allowed. Downloaded to", local_path_probe)
      }, error = function(e) {
        paste("NG: GetObject/save_object failed:", conditionMessage(e))
      })

      cat(res_get, "\n")
      cat("---------------------\n")
    })

    paste(txt_lines, collapse = "\n")
  }

  # Execute Probe
  log_text_S3_1 <- capture_s3_get_probe(bucket = bucket, key = key)

} else {
  # When AMAZON_FLAG is FALSE
  ID_exclude <- character(0)
  # Define empty log variables to prevent renderText errors
  EXCLUDED_TEXT_1 <- "AWS Flag is OFF"
  EXCLUDED_TEXT_2 <- ""
  Identity_text_1 <- ""
  log_text_S3_1   <- ""
}

# ------------------------------------------------------------------------------
# 7. Render Output
# ------------------------------------------------------------------------------

output$log_text_output <- renderText({
  if(AMAZON_FLAG){
    log_text <- paste(EXCLUDED_TEXT_1, EXCLUDED_TEXT_2, Identity_text_1, log_text_S3_1, sep = "\n\n")
    return(log_text)
  }
})
