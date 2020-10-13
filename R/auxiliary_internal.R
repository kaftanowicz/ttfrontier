##########################################################################################
# Auxiliary functions (internal)
##########################################################################################

#' Compute dependent variable without inefficiencies
#'
#' @param x numeric matrix of values of dependent variables
#' @param coeffs numeric vector of coefficients, including the intercept as the first element
#'
#' @return numeric vector of dependent variable without inefficiencies
#' @noRd
get_mu <- function(x, coeffs) {
  coeffs[1L] + x %*% matrix(tail(coeffs, -1L), ncol = 1L)
}

#' Compute compound error terms
#'
#' @param variables vector of variables (coefficients and error parameters)
#' @param model_mx model matrix
#' @param y vector of dependent variables
#'
#' @return vector of compound error terms
#' @noRd
get_epsilon <- function(variables, model_mx, y) {
  var_sep <- separate_vecpar(variables = variables,
                             error_term_param_names = c("sigma_v", "sigma_u", "sigma_w"))
  return(y - as.vector(model_mx %*% var_sep$coeffs))
}

#' Check log-likelihood functions returns proper results for boundaries and start values.
#' Raises error on negative check result.
#'
#' @param boundaries named list of numeric vectors that will be proper arguments of LL function
#' @param neg_ll negative log-likelihood function
#' @param model_mx numeric model matrix (including intercept column of 1s)
#' @param y numeric vector of values of independent variable
#'
#' @return invisible(TRUE) if LL returns only finite numeric values
#' @noRd
check_variables_boundaries <- function(boundaries, neg_ll, model_mx, y) {
  for (values_kind in names(boundaries)) {
    values <- boundaries[[values_kind]]
    values_report <- paste(names(values), values, sep = ": ", collapse = ", ")
    neg_ll_val <- neg_ll(values, model_mx, y)
    if (!is.finite(neg_ll_val)) {
      stop(sprintf("Non-finite value (%s) of negative log likelihood function for %s values: %s",
                   as.character(neg_ll_val), values_kind, values_report))
    }
  }
  return(invisible(TRUE))
}

#' Separate a vector of variables into a list containing a vector of coefficients
#' and a named list of error term parameters.
#'
#' @param variables vector of variables of the form c(variables, error_term_params)
#' @param error_term_param_names vector of names of error term parameters
#'
#' @return 2-element list: a vector of coefficients and a named list of error term parameters
#' @noRd
separate_vecpar <- function(variables,
                            error_term_param_names = c("sigma_v", "sigma_u", "sigma_w")) {
  n_error_terms <- length(error_term_param_names)
  coeffs <- head(variables, -n_error_terms)
  error_term_params <- as.list(tail(variables, n_error_terms))
  names(error_term_params) <- error_term_param_names
  return(list(coeffs = coeffs, error_term_params = error_term_params))
}
