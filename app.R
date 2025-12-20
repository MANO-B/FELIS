# Sys.setenv(FELIS_DATA_ROOT = "/srv/shiny-server/felis-cs")
Sys.setenv(FELIS_DATA_ROOT = getwd())
APP_DIR <- system.file("app", package = "FELIS")
setwd(APP_DIR)
source("app.R", local = TRUE, chdir = TRUE)
# Run the application
shinyApp(ui = ui, server = server)
