#' Run FELIS Shiny app
#' @export

run_app <- function(data_root = NULL, ...) {
  app_dir <- system.file("app", package = "FELIS")
  shiny::runApp(app_dir, launch.browser = F, ...)
}
