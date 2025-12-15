#' Run FELIS Shiny app
#' @export
run_app <- function(...) {
  app_dir <- system.file("app", package = "FELIS")
  if (app_dir == "") stop("App directory not found. Is inst/app included in the package?")

  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(app_dir)

  shiny::runApp(app_dir, ...)
}
