SOFT_VERSION = "4.2.8"
DATA_VERSION = "20251024"
bucket <- "3e040086-data-production-felisanalysisportal"
key    <- "ic_withdrawn_list.txt"
app_dir <- Sys.getenv("FELIS_DATA_ROOT", unset = getOption("felis_data_root", tempdir()))
CCAT_FLAG <- file.exists(file.path(app_dir, "ccat"))
AMAZON_FLAG <- file.exists("/srv/shiny-server/felis-ccat/AMAZON")
addResourcePath(prefix = "APP_DIR", directoryPath = APP_DIR)
app_path <- function(...) {
  file.path(APP_DIR, ...)
}

shiny::addResourcePath("felis-src", file.path(APP_DIR, "source"))

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
prefer <- list(
  gtsummary = c("all_double","all_factor","all_integer","all_logical","all_numeric","as_flextable","continuous_summary"),
  patchwork = c("area"),
  dplyr = c("arrange","count","desc","failwith","id","lag","mutate","rename","select","src","summarise","summarize","between","first","last"),
  flextable = c("border","compose","font","rotate"),
  shinydashboard = c("box"),
  readr = c("col_factor"),
  purrr = c("compact","discard","transpose","when","flatten","flatten_chr","flatten_dbl","flatten_int","flatten_lgl","flatten_raw","splice"),
  gt = c("html"),
  survival = c("myeloma","cluster"),
  shiny = c("observe","runExample"),
  tidyr = c("pack","unpack","replace_na","smiths"),
  recipes = c("prepare","step","update"),
  dials = c("prune","smoothness"),
  yardstick = c("spec","mcc","rmse","sensitivity"),
  parsnip = c("translate"),
  rms = c("validate"),
  fastshap = c("explain"),
  bigstep = c("prepare_data"),
  DT = c("dataTableOutput","renderDataTable"),
  data.table = c("dcast","yearmon","yearqtr",":="),
  shinyWidgets = c("alert"),
  methods = c("removeClass"),
  rstan = c("show"),
  datasets = c("penguins"),
  ranger = c("timepoints")
)
invisible(purrr::iwalk(prefer, \(funs, pkg) {
  purrr::walk(funs, \(fn) {
    conflicted::conflicts_prefer(.quiet = TRUE, !!rlang::parse_expr(paste0(pkg, "::`", fn, "`")))
  })
}))
