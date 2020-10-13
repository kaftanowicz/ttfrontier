##########################################################################################
# Auxiliary functions that are still useful enough to be exported
##########################################################################################

#' Create a list whose elements have the same names as the variables in input.
#'
#' @param ... objects
#'
#' @return a named list containing given objects
#' @export
#'
#' @examples
#' \dontrun{
#' x <- 10
#' y <- "letter y"
#' named_list(x, y) # same result as list(x = x, y = y)
#' }
named_list <- function(...) {
  var_list <- list(...)
  var_names <- sapply(substitute(list(...)), deparse)[-1] # -1 drops the list
  names(var_list) <- var_names
  return(var_list)
}

#' Get parameter scales for optimization
#'
#' Results are useful as an argument to the `parscale` component of `control` argument of `optim` function (and other optimization functions).
#'
#' @param start numeric vector of start values of coefficients
#'
#' @return numeric vector of scales, each being a power of 10
#' @export
#'
#' @examples
#' \dontrun{
#' get_parscale(c("(Intercept)" = 3.8, "educ" = 0.142, "exper" = 0.027, "tenure" = 0.021))
#' }
get_parscale <- function(start) {
  10^floor(log10(abs(start)))
}
