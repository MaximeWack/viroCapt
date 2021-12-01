#' Shiny application to interactively visualize results from the pipeline
#'
#' @export
visu <- function(...) {
  appDir <- system.file("visu", package = "viroCapt")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `viroCapt`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", ...)
}
