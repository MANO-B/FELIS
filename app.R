# Sys.setenv(FELIS_DATA_ROOT = "/srv/shiny-server/felis-ccat")
# Sys.setenv(FELIS_DATA_ROOT = getwd())
# FELIS::run_app()


#Sys.setenv(FELIS_DATA_ROOT = "/srv/shiny-server/felis-cs")
Sys.setenv(FELIS_DATA_ROOT = "/Users/ikegami/Desktop/length_bias/felis-ccat")
APP_DIR <- system.file("app", package = "FELIS")
stopifnot(nzchar(APP_DIR))

old <- getwd()
on.exit(setwd(old), add = TRUE)
setwd(APP_DIR)

source("app.R", local = TRUE, chdir = TRUE)
# Run the application
shinyApp(ui = ui, server = server)
