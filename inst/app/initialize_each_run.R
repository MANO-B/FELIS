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
object_key    <- "ic_withdrawn_list.txt"
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
  if (!is.null(response) && httr::status_code(response) == 200) {

    # Parse JSON content
    content_text <- httr::content(response, as = "text", encoding = "UTF-8")

    creds_list <- tryCatch({
      parsed <- jsonlite::fromJSON(content_text)

      # Verify required keys exist
      required_keys <- c("AccessKeyId", "SecretAccessKey", "Token")
      if (!all(required_keys %in% names(parsed))) {
        stop("JSON missing required credential keys.")
      }

      parsed
    }, error = function(e) {
      return(NULL)
    })
  }

  # ----------------------------------------------------------------------------
  # 2. Set Credentials to Variables (Handling Failures)
  # ----------------------------------------------------------------------------

  if(is.null(creds_list)){
    # Fallback to dummy credentials if retrieval failed
    # (Operations will fail, but the app won't crash immediately)
    MY_AWS_ACCESS_KEY_ID    <- "DUMMY_ACCESS_KEY"
    MY_AWS_SECRET_ACCESS_KEY    <- "DUMMY_SECRET_KEY"
    MY_AWS_SESSION_TOKEN <- "DUMMY_TOKEN"
    MY_AWS_REGION            <- "ap-northeast-1"
  } else {
    MY_AWS_ACCESS_KEY_ID    <- creds_list$AccessKeyId
    MY_AWS_SECRET_ACCESS_KEY    <- creds_list$SecretAccessKey
    MY_AWS_SESSION_TOKEN <- creds_list$Token
    MY_AWS_REGION            <- "ap-northeast-1"
  }

  # Also set as environment variables for compatibility with other libs/logging
  Sys.setenv("AWS_ACCESS_KEY_ID"     = MY_AWS_ACCESS_KEY_ID)
  Sys.setenv("AWS_SECRET_ACCESS_KEY" = MY_AWS_SECRET_ACCESS_KEY)
  Sys.setenv("AWS_SESSION_TOKEN"     = MY_AWS_SESSION_TOKEN)
  Sys.setenv("AWS_DEFAULT_REGION"    = MY_AWS_REGION)
  Sys.setenv("AWS_REGION"            = MY_AWS_REGION)

  # ----------------------------------------------------------------------------
  # 4. Define Download Function with Explicit Credentials
  # ----------------------------------------------------------------------------

  download_withdrawn_list_safe <- function(bucket, object_key,
                                           local_path = file.path(tempdir(), "ic_withdrawn_list.txt")) {

    if(file.exists(local_path)){
      file.remove(local_path)
    }

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
        object = object_key,
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
    }, error = function(e) {
      character(0)
    })

    # Safe fallback: create empty file if failed
    if (!success && !file.exists(local_path)) {
      file.create(local_path)
    }

    return(local_path)
  }

  # ----------------------------------------------------------------------------
  # 5. Execute Download and Read File
  # ----------------------------------------------------------------------------

  withdrawn_file <- download_withdrawn_list_safe(bucket = bucket, object_key = object_key)

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
    character(0)
  })
} else {
  # When AMAZON_FLAG is FALSE
  ID_exclude <- character(0)
}
