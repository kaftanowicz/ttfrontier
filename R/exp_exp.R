
#' Get start values and lower and upper boundaries for variables
#'
#' @param start_val numeric vector of start values
#' @param error_param_min numeric minimum value of error parameters
#' @param error_param_max numeric maximum value of error parameters
#' @param r non-negative numeric range parameter; the greater, the wider are boundaries for variables
#' @param r_intercept non-negative numeric range parameter; the greater, the wider are boundaries for intercept
#'
#' @return 3-element named list: `start` vector of start values, `lower` vector of minimum values and `upper` vector of maximum values of variables
#'
#' @noRd
get_variables_boundaries_exp_exp <- function(start_val,
                                             error_param_min = 0.01,
                                             error_param_max = 0.5,
                                             r = 0.3,
                                             r_intercept = 0.5) {
  res <- list()
  n_vals <- length(start_val)
  sigma_start <- c("sigma_v" = 0.1, "sigma_u" = 0.1, "sigma_w" = 0.1)
  n_sigma <- length(sigma_start)
  min_sigma <- error_param_min
  max_sigma <- error_param_max
  res$start <- c(start_val, sigma_start)
  smaller_abs = c((1 - r_intercept) * start_val[1],
                  (1 - r) * start_val[2:n_vals],
                  rep(min_sigma, n_sigma))
  greater_abs = c((1 + r_intercept) * start_val[1],
                  (1 + r) * start_val[2:n_vals],
                  rep(max_sigma, n_sigma))
  res$lower <- pmin(smaller_abs, greater_abs)
  res$upper <- pmax(smaller_abs, greater_abs)

  names(res$lower) <- names(res$start)
  names(res$upper) <- names(res$start)
  return(res)
}


#' Negative log-likelihood function for N +Exp -Exp compound error term
#'
#' @param variables numeric vector of coefficients, last 3 of which are parameters of error distributions
#' @param model_mx numeric matrix model
#' @param y numeric vector of values of dependent variable
#'
#' @return negative log likelihood of model given data
#' @noRd
neg_log_likelihood_exp_exp <- function(variables, model_mx, y) {
  # todo: make it only have 1 vector parameter for compatibility with mle2
  var_sep <- separate_vecpar(variables)
  sigma_v <- var_sep$error_term_params$sigma_v
  sigma_u <- var_sep$error_term_params$sigma_u
  sigma_w <- var_sep$error_term_params$sigma_w

  n <- nrow(model_mx)
  epsilon <- y - model_mx %*% var_sep$coeffs
  return(-1 * log_likelihood_exp_exp(n, epsilon, sigma_v, sigma_u, sigma_w))
}

#' Log-likelihood function of two-tier model
#'
#' Based on Kumbhakar, Subal & Parmeter, Christopher (2009).
#'
#' References:
#' Kumbhakar, Subal & Parmeter, Christopher. (2009). The effects of match uncertainty and bargaining on labor market outcomes: Evidence from firm and worker specific estimates. Journal of Productivity Analysis. 31. 1-14. 10.1007/s11123-008-0117-3.
#'
#' @param n integer number of observations
#' @param epsilon numeric error value
#' @param sigma_v numeric standard deviation of symmetric normal error term
#' @param sigma_u numeric mean (expected value) of negative exponentially distributed error term
#' @param sigma_w numeric mean (expected value) of positive exponentially distributed error term
#'
#' @return numeric value of log-likelihood
#' @noRd
log_likelihood_exp_exp <- function(n, epsilon,
                                   sigma_v, sigma_u, sigma_w) {
  par <- parametrize_exp_exp(epsilon, sigma_v, sigma_u, sigma_w)
  return(-1 * n * log(sigma_u + sigma_w) +
    sum(
      log(
        exp(par$alpha) * pnorm(par$beta) + exp(par$a) * pnorm(par$b)
      )
    ))
}

#' Calculate auxiliary parameters for N+Exp-Exp model
#'
#' @param epsilon numeric error value
#' @param sigma_v numeric standard deviation of symmetric normal error term
#' @param sigma_u numeric mean (expected value) of negative exponentially distributed error term
#' @param sigma_w numeric mean (expected value) of positive exponentially distributed error term
#'
#' @return named list of auxiliary parameters
#' @noRd
parametrize_exp_exp <- function(epsilon, sigma_v, sigma_u, sigma_w) {
  a <- sigma_v^2 / (2 * sigma_w^2) - epsilon / sigma_w
  b <- epsilon / sigma_v - sigma_v / sigma_w
  alpha <- epsilon / sigma_u + sigma_v^2 / (2 * sigma_u^2)
  beta <- -1 * (epsilon / sigma_v + sigma_v / sigma_u)
  lambda <- 1/sigma_u + 1/sigma_w
  chi_1 <- pnorm(b) + exp(alpha - a) * pnorm(beta)
  chi_2 <- pnorm(beta) + exp(a - alpha) * pnorm(b)
  return(named_list(a, b, alpha, beta, lambda, chi_1, chi_2))
}

#' Compute composite error term
#'
#' @param n integer number of values to return
#' @param sigma_v numeric standard deviation of symmetric normal error term
#' @param sigma_u numeric mean (expected value) of negative exponentially distributed error term
#' @param sigma_w numeric mean (expected value) of positive exponentially distributed error term
#'
#' @return numeric vector of error values
#' @noRd
get_epsilon_exp_exp <- function(n, sigma_v, sigma_u, sigma_w) {
  rnorm(n, mean = 0, sd = sigma_v) - rexp(n, rate = 1/sigma_u) + rexp(n, rate = 1/sigma_w)
}

#' Compute observation-specific extracted surplus
#'
#' Computes observation-specific point estimates of extracted surplus (u and w) from the composed
#' error term.
#'
#' @param variables vector of variables (coefficients and error parameters)
#' @param model_mx model matrix
#' @param y vector of dependent variables
#'
#' @return
#' \describe{
#'   \item{u}{negative extracted surplus = E(u|epsilon)}
#'   \item{w}{positive extracted surplus = E(w|epsilon)}
#'   \item{exp_neg_u}{percentage measure of price reduction = E(exp(-u)|epsilon), if log(price) is dependent variable}
#'   \item{exp_neg_w}{percentage measure of price increase = E(exp(-w)|epsilon), if log(price) is dependent variable}
#' }
#' @noRd
estimate_extracted_surplus_exp_exp <- function(variables, formula, data) {
  model_mx <- model.matrix(object = formula, data = data)
  y <- data[[head(all.vars(formula), 1)]] # values of dependent variable

  epsilon <- get_epsilon(variables, model_mx, y)
  sigma_v <- variables[["sigma_v"]]
  sigma_u <- variables[["sigma_u"]]
  sigma_w <- variables[["sigma_w"]]

  par <- parametrize_exp_exp(epsilon, sigma_v, sigma_u, sigma_w)

  u <- 1/par$lambda +
    (exp(par$alpha - par$a) * sigma_v * (dnorm(-par$beta) + par$beta * pnorm(par$beta))) / par$chi_1

  w <- 1/par$lambda +
    (sigma_v * (dnorm(-par$b) + par$b * pnorm(par$b))) / par$chi_1

  exp_neg_u <- (par$lambda / (1 + par$lambda)) *
    (1 / par$chi_2) * (
      pnorm(par$b) +
        exp(par$alpha - par$a) *
        exp(sigma_v ^ 2 / 2 - sigma_v * par$beta) *
        pnorm(par$beta - sigma_v)
    )
  exp_neg_w <- (par$lambda / (1 + par$lambda)) *
    (1 / par$chi_1) * (
      pnorm(par$beta) +
        exp(par$a - par$alpha) *
        exp(sigma_v ^ 2 / 2 - sigma_v * par$beta) *
        pnorm(par$b - sigma_v)
    )
  return(named_list(u, w, exp_neg_u, exp_neg_w))
}

