#' Shiny application to interactively visualize results from the pipeline
#'
#' @export
visu <- function() {
  appDir <- system.file("visu", package = "HPVcap")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `HPVcap`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
