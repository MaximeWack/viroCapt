#' Replace invalid values in a (NA, NaN, NULL) with default value b
#'
#' @param a Any value
#' @param b A default return value
#' @return Either a, or b if a is any of NULL, NA, NaN
`%||%` <- function(a, b) if (length(a) == 0) b else if (a %>% is.na) b else a


#' Compose two functions (sequence order, not composition order)
#'
#' @param f A function a -> b
#' @param g A function b -> c
#' @return A function a -> c
compose2 <- function(f, g)
{
  function(...)
  {
    g(f(...))
  }
}


#' Compose a list of functions
#'
#' @param funs A list of functions, in sequence order
#' @return A composed function
compose <- function(funs)
{
  Reduce(compose2, funs)
}


#' Compose all combinations of the given named functions (left to right order)
#'
#' @param ... Named functions to compose-combine
#' @return A list of functions, result of the compositions of all ordered combinations of inputted functions
compose_combn <- function(...)
{
  funs <- lapply(list(...), match.fun)

  fun_names <- names(funs)
  sizes <- length(funs)

  1:sizes %>%
  lapply(utils::combn, x = fun_names) %>%
    lapply(apply, 2, function(x) {paste0(x, collapse = "")}) %>%
    unlist -> comb_names

  1:sizes %>%
    lapply(utils::combn, x = funs) %>%
    lapply(apply, 2, compose) %>%
    unlist %>%
    stats::setNames(comb_names)
}
