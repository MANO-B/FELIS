#' Run FELIS Shiny app
#' @export

run_app <- function() {
  FELIS_dir <- system.file("app", package = "FELIS")
  shiny::runApp(FELIS_dir, launch.browser = F)
}
