#' Run FELIS Shiny app
#' @export
run_app <- function(...) {
  app_dir <- system.file("shiny", package = "FELIS")
  stopifnot(nzchar(app_dir))

  old <- options(felis_app_dir = app_dir)
  on.exit(options(old), add = TRUE)
  # old <- getwd()
  # on.exit(setwd(old), add = TRUE)
  # setwd(app_dir)

  shiny::runApp(app_dir, ...)
}
