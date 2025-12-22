SOFT_VERSION = "4.2.6"
DATA_VERSION = "20251024"
bucket <- "3e040086-data-production-felisanalysisportal"
key    <- "ic_withdrawn_list.txt"
Sys.setenv(AWS_EC2_METADATA_DISABLED = "true")
app_dir <- Sys.getenv("FELIS_DATA_ROOT", unset = getOption("felis_data_root", tempdir()))
CCAT_FLAG <- TRUE
# CCAT_FLAG <- file.exists("ccat")

app_path <- function(...) {
  file.path(APP_DIR, ...)
}

shiny::addResourcePath("felis-src", file.path(APP_DIR, "source"))

if(file.exists(file.path(tempdir(), "Data_mutation_cord.qs"))){
  file.remove(file.path(tempdir(), "Data_mutation_cord.qs"))
}
if(file.exists(file.path(tempdir(), "m2.qs"))){
  file.remove(file.path(tempdir(), "m2.qs"))
}
if(file.exists(file.path(tempdir(), "lgb_m2.qs"))){
  file.remove(file.path(tempdir(), "lgb_m2.qs"))
}
if(file.exists(file.path(tempdir(), "rf_m2.qs"))){
  file.remove(file.path(tempdir(), "rf_m2.qs"))
}
if(file.exists(file.path(tempdir(), "m2_rec.qs"))){
  file.remove(file.path(tempdir(), "m2_rec.qs"))
}
if(file.exists(file.path(tempdir(), "lgb_m2_rec.qs"))){
  file.remove(file.path(tempdir(), "lgb_m2_rec.qs"))
}
if(file.exists(file.path(tempdir(), "rf_m2_rec.qs"))){
  file.remove(file.path(tempdir(), "rf_m2_rec.qs"))
}
if(file.exists(file.path(tempdir(), "m2_GMT.qs"))){
  file.remove(file.path(tempdir(), "m2_GMT.qs"))
}
if(file.exists(file.path(tempdir(), "lgb_m2_GMT.qs"))){
  file.remove(file.path(tempdir(), "lgb_m2_GMT.qs"))
}
if(file.exists(file.path(tempdir(), "rf_m2_GMT.qs"))){
  file.remove(file.path(tempdir(), "rf_m2_GMT.qs"))
}
if(file.exists(file.path(tempdir(), "PS_data.rda"))){
  file.remove(file.path(tempdir(), "PS_data.rda"))
}
if(file.exists(file.path(tempdir(), "Organ_data.rda"))){
  file.remove(file.path(tempdir(), "Organ_data.rda"))
}
if(file.exists(file.path(tempdir(), "Lines_data.rda"))){
  file.remove(file.path(tempdir(), "Lines_data.rda"))
}

stan_model_factor_code_path <- "source/Stan_code_loglog_weibull_factor.stan"
stan_model_factor_save_path_rstan <- file.path(app_dir, "compiled_model_factor_rstan.rds")
stan_model_factor_save_path_cmdstan <- file.path(app_dir, "compiled_model_factor_cmdstan.rds")

stan_model_simple_code_path <- "source/Stan_code_loglog_weibull.stan"
stan_model_simple_save_path_rstan <- file.path(app_dir, "compiled_model_rstan.rds")
stan_model_simple_save_path_cmdstan <- file.path(app_dir, "compiled_model_cmdstan.rds")

if(Sys.getenv("SHINY_SERVER_VERSION") != ""){
  DOCKER = FALSE
  PARALLEL = min(6, max(1, parallel::detectCores(), na.rm = TRUE))
  MAX_CORE = 1
  NO_IMPROVE = 5
  ITER = 1000
  BootNoVar=3
  if(!file.exists(stan_model_factor_save_path_rstan)){
    mod_factor <- stan_model(file = stan_model_factor_code_path)
    saveRDS(mod_factor, stan_model_factor_save_path_rstan)
    rm(mod_factor)
  }

  if(!file.exists(stan_model_simple_save_path_rstan)){
    mod_simple <- stan_model(file = stan_model_simple_code_path)
    saveRDS(mod_simple, stan_model_simple_save_path_rstan)
    rm(mod_simple)
  }
  ENV_ = "server"
} else if(file.exists("docker")){
  library(cmdstanr)
  DOCKER = TRUE
  PARALLEL = min(6, max(1, parallel::detectCores(), na.rm = TRUE))
  MAX_CORE = PARALLEL
  NO_IMPROVE = 5
  ITER = 1000
  BootNoVar=1

  if(!file.exists(stan_model_factor_save_path_cmdstan)){
    stan_model_compiled <- cmdstan_model(stan_file = stan_model_factor_code_path)
    saveRDS(stan_model_compiled, stan_model_factor_save_path_cmdstan)
    rm(stan_model_compiled)
  }

  if(!file.exists(stan_model_simple_save_path_cmdstan)){
    stan_model_compiled <- cmdstan_model(stan_file = stan_model_simple_code_path)
    saveRDS(stan_model_compiled, stan_model_simple_save_path_cmdstan)
    rm(stan_model_compiled)
  }
  ENV_ = "docker"
} else {
  library(cmdstanr)
  DOCKER = TRUE
  PARALLEL = min(6, max(1, parallel::detectCores(), na.rm = TRUE))
  MAX_CORE = PARALLEL
  NO_IMPROVE = 5
  ITER = 1000
  BootNoVar=1

  hash_file_path <- paste0(stan_model_factor_save_path_cmdstan, ".md5")
  current_hash <- md5sum(stan_model_factor_code_path)
  needs_recompile <- !file.exists(stan_model_factor_save_path_cmdstan) ||
    !file.exists(hash_file_path) ||
    readLines(hash_file_path) != current_hash
  if (needs_recompile) {
    stan_model_compiled <- cmdstan_model(stan_file = stan_model_factor_code_path)
    saveRDS(stan_model_compiled, stan_model_factor_save_path_cmdstan)
    rm(stan_model_compiled)
    writeLines(current_hash, hash_file_path)
  }
  hash_file_path <- paste0(stan_model_simple_save_path_cmdstan, ".md5")
  current_hash <- md5sum(stan_model_simple_code_path)
  needs_recompile <- !file.exists(stan_model_simple_save_path_cmdstan) ||
    !file.exists(hash_file_path) ||
    readLines(hash_file_path) != current_hash
  if (needs_recompile) {
    stan_model_compiled <- cmdstan_model(stan_file = stan_model_simple_code_path)
    saveRDS(stan_model_compiled, stan_model_simple_save_path_cmdstan)
    rm(stan_model_compiled)
    writeLines(current_hash, hash_file_path)
  }
  ENV_ = "local"
}

tidymodels_prefer()
options(pillar.advice = FALSE, pillar.min_title_chars = Inf)

nietzsche = read.csv(header = TRUE,
                     file("source/nietzsche.csv",
                          encoding='UTF-8-BOM'))$text
download_withdrawn_list_safe <- function(local_path = file.path(tempdir(), "ic_withdrawn_list.txt"),
                                         bucket, key,
                                         quiet = TRUE,
                                         connect_timeout_sec = 0.5,
                                         total_timeout_sec   = 1) {
  dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)
  ensure_empty <- function() {
    writeLines(character(0), con = local_path, useBytes = TRUE)
  }

  cfg <- httr::config(connecttimeout = connect_timeout_sec)
  to  <- httr::timeout(total_timeout_sec)

  exists <- tryCatch({
    !is.null(head_object(object = key, bucket = bucket, config = cfg, to))
  }, error = function(e) {
    if (!quiet) message("head_object failed; using empty withdrawn list: ", conditionMessage(e))
    FALSE
  })

  if (exists) {
    ok <- tryCatch({
      save_object(object = key, bucket = bucket, file = local_path, config = cfg, to)
      TRUE
    }, error = function(e) {
      if (!quiet) message("save_object failed; using empty withdrawn list: ", conditionMessage(e))
      FALSE
    })
    if (!ok) ensure_empty()
  } else {
    ensure_empty()
  }

  local_path
}
withdrawn_file <- download_withdrawn_list_safe()

ID_exclude <- tryCatch({
  if (!file.exists(withdrawn_file) || is.na(file.info(withdrawn_file)$size) || file.info(withdrawn_file)$size == 0) {
    character(0)
  } else {
    Encode <- tryCatch(
      readr::guess_encoding(withdrawn_file, n_max = 1000)[[1]]$encoding,
      error = function(e) "UTF-8"
    )
    if (identical(Encode, "Shift_JIS")) Encode <- "CP932"

    readr::read_csv(
      file = withdrawn_file,
      col_names = "X1",
      locale = readr::locale(encoding = Encode),
      progress = FALSE,
      show_col_types = FALSE
    )$X1
  }
}, error = function(e) character(0))

# Rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = PARALLEL)
options(expressions = 500000)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)

ht_opt$message = FALSE

Penal = 1
Toler = 10^-14
GENE_NO_THRESHOLD = 20

# 設定パラメータ
ANALYSIS_CONFIG <- list(
  seed = 1234,
  chains = 4,
  iter_sampling = ITER,
  iter_warmup = as.integer(ITER / 2),
  adapt_delta = 0.99,
  max_treedepth = 10,
  thin = 1,
  refresh = 500,
  parallel_chains = PARALLEL
)


options(shiny.maxRequestSize=16*1024^3)
conflicted::conflicts_prefer(.quiet = T, gtsummary::all_double)
conflicted::conflicts_prefer(.quiet = T, gtsummary::all_factor)
conflicted::conflicts_prefer(.quiet = T, gtsummary::all_integer)
conflicted::conflicts_prefer(.quiet = T, gtsummary::all_logical)
conflicted::conflicts_prefer(.quiet = T, gtsummary::all_numeric)
conflicted::conflicts_prefer(.quiet = T, patchwork::area)
conflicted::conflicts_prefer(.quiet = T, dplyr::arrange)
#conflicted::conflicts_prefer(.quiet = T, sjlabelled::as_factor)
conflicted::conflicts_prefer(.quiet = T, gtsummary::as_flextable)
#conflicted::conflicts_prefer(.quiet = T, sjlabelled::as_label)
conflicted::conflicts_prefer(.quiet = T, flextable::border)
conflicted::conflicts_prefer(.quiet = T, shinydashboard::box)
conflicted::conflicts_prefer(.quiet = T, readr::col_factor)
conflicted::conflicts_prefer(.quiet = T, purrr::compact)
conflicted::conflicts_prefer(.quiet = T, flextable::compose)
conflicted::conflicts_prefer(.quiet = T, gtsummary::continuous_summary)
conflicted::conflicts_prefer(.quiet = T, dplyr::count)
conflicted::conflicts_prefer(.quiet = T, DT::dataTableOutput)
conflicted::conflicts_prefer(.quiet = T, dplyr::desc)
conflicted::conflicts_prefer(.quiet = T, purrr::discard)
conflicted::conflicts_prefer(.quiet = T, dplyr::failwith)
conflicted::conflicts_prefer(.quiet = T, stringr::fixed)
conflicted::conflicts_prefer(.quiet = T, flextable::font)
conflicted::conflicts_prefer(.quiet = T, gt::html)
conflicted::conflicts_prefer(.quiet = T, dplyr::id)
#conflicted::conflicts_prefer(.quiet = T, plyr::is.discrete)
conflicted::conflicts_prefer(.quiet = T, dplyr::lag)
conflicted::conflicts_prefer(.quiet = T, dplyr::mutate)
conflicted::conflicts_prefer(.quiet = T, survival::myeloma)
conflicted::conflicts_prefer(.quiet = T, shiny::observe)
conflicted::conflicts_prefer(.quiet = T, tidyr::pack)
conflicted::conflicts_prefer(.quiet = T, recipes::prepare)
conflicted::conflicts_prefer(.quiet = T, dials::prune)
conflicted::conflicts_prefer(.quiet = T, dplyr::rename)
conflicted::conflicts_prefer(.quiet = T, DT::renderDataTable)
conflicted::conflicts_prefer(.quiet = T, flextable::rotate)
conflicted::conflicts_prefer(.quiet = T, dplyr::select)
conflicted::conflicts_prefer(.quiet = T, dials::smoothness)
conflicted::conflicts_prefer(.quiet = T, yardstick::spec)
conflicted::conflicts_prefer(.quiet = T, dplyr::src)
conflicted::conflicts_prefer(.quiet = T, recipes::step)
conflicted::conflicts_prefer(.quiet = T, dplyr::summarise)
conflicted::conflicts_prefer(.quiet = T, dplyr::summarize)
conflicted::conflicts_prefer(.quiet = T, parsnip::translate)
conflicted::conflicts_prefer(.quiet = T, tidyr::unpack)
conflicted::conflicts_prefer(.quiet = T, recipes::update)
conflicted::conflicts_prefer(.quiet = T, rms::validate)
conflicted::conflicts_prefer(.quiet = T, dplyr::between)
conflicted::conflicts_prefer(.quiet = T, dplyr::first)
conflicted::conflicts_prefer(.quiet = T, dplyr::last)
conflicted::conflicts_prefer(.quiet = T, purrr::transpose)
#conflicted::conflicts_prefer(.quiet = T, foreach::accumulate)
conflicted::conflicts_prefer(.quiet = T, fastshap::explain)
conflicted::conflicts_prefer(.quiet = T, tidyr::replace_na)
conflicted::conflicts_prefer(.quiet = T, bigstep::prepare_data)
conflicted::conflicts_prefer(.quiet = T, survival::cluster)
conflicted::conflicts_prefer(.quiet = T, data.table::dcast)
conflicted::conflicts_prefer(.quiet = T, yardstick::mcc)
conflicted::conflicts_prefer(.quiet = T, yardstick::rmse)
conflicted::conflicts_prefer(.quiet = T, yardstick::sensitivity)
conflicted::conflicts_prefer(.quiet = T, tidyr::smiths)
conflicted::conflicts_prefer(.quiet = T, purrr::when)
conflicted::conflicts_prefer(.quiet = T, shinyWidgets::alert)
conflicted::conflicts_prefer(.quiet = T, methods::removeClass)
conflicted::conflicts_prefer(.quiet = T, shiny::runExample)
conflicted::conflicts_prefer(.quiet = T, rstan::show)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten_chr)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten_dbl)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten_int)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten_lgl)
conflicted::conflicts_prefer(.quiet = T, purrr::flatten_raw)
conflicted::conflicts_prefer(.quiet = T, datasets::penguins)
conflicted::conflicts_prefer(.quiet = T, purrr::splice)
conflicted::conflicts_prefer(.quiet = T, ranger::timepoints)
conflicted::conflicts_prefer(.quiet = T, data.table::yearmon)
conflicted::conflicts_prefer(.quiet = T, data.table::yearqtr)
conflicted::conflicts_prefer(.quiet = T, data.table::`:=`)
