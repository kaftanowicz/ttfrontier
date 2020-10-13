##########################################################################################
# Main functions (exported)
##########################################################################################

#' Two-Tier Stochastic Frontier Analysis
#'
#' @param formula a symbolic description of the model to be estimated
#' @param data data.frame that contains cross-section data
#' @param start_val numerical vector of initial values for optimizer; if NULL (default), will use OLS estimates
#' @param ineff_distr 2-element character vector specifying the distributions of the inefficiency (asymmetric) error terms; currently only `c("exp", "exp")` is supported
#' @param trace	integer in range 0:6; higher values cause more information on the progress of the optimization to be printed
#' @return
#' \describe{
#'   \item{mle_fit}{an object of class `mle2`; see \code{\link[bbmle]{mle2}}}
#'   \item{ineff_distr}{equal to argument `ineff_distr`}
#'   \item{unexplained_variation}{numeric value of unexplained variation in dependent variable, in range [0, 1]}
#'   \item{unexplained_variation_from_inefficiency}{numeric value of ratio of unexplained variation that comes from inefficiency (asymmetric error terms), in range [0, 1]}
#'   \item{surplus}{list containing observation-specific surplus estimates from the composed error term:
#'   \describe{
#'      \item{u}{negative extracted surplus = E(u|epsilon)}
#'      \item{w}{positive extracted surplus = E(w|epsilon)}
#'      \item{exp_neg_u}{percentage measure of price reduction = E(exp(-u)|epsilon), if log(price) is dependent variable}
#'      \item{exp_neg_w}{percentage measure of price increase = E(exp(-w)|epsilon), if log(price) is dependent variable}
#'      }
#'   }
#'   \item{ols_fit}{an object of class `lm`; see \code{\link[stats]{lm}}}
#'   \item{formula}{equal to argument `formula`}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Replication of Kumbhakar, S. C., & Parmeter, C. F. (2009), <doi:10.1007/s11123-008-0117-3>
#' library(wooldridge) # for dataset wage2
#' data(wage2)
#' d <- na.omit(wage2)
#' d$educ_sq <- d$educ^2
#' d$exper_sq <- d$exper^2
#' d$tenure_sq <- d$tenure^2
#'
#' fit <- ttsfa(formula = lwage ~ educ + educ_sq + exper + exper_sq + tenure + tenure_sq +
#'                age + IQ + urban + married + south + sibs + brthord + meduc + feduc,
#'              data = d)
#' summary(fit)
#' }
ttsfa <- function(formula,
                  data,
                  start_val = NULL,
                  ineff_distr = c("exp", "exp"),
                  error_param_min = 0.01,
                  error_param_max = 0.5,
                  maxit = 1e5,
                  trace = 0L) {
  # todo: # ttsfa doesn't work properly when you specify function inside formula, e.g. log(y) ~ .
  model_mx <- model.matrix(object = formula, data = data)
  y <- data[[head(all.vars(formula), 1)]] # values of dependent variable

  # OLS estimation
  ols_fit <- lm(formula = formula, data = data)
  if (is.null(start_val)) {
    start_val <- ols_fit$coefficients
  } else {
    len_sv <- length(start_val)
    ncol_mm <- ncol(model_mx)
    if (len_sv != ncol_mm) {
      stop(sprintf("Length of provided start values vector (%d) is not equal to the number of columns of model matrix (%d)",
                   len_sv, ncol_mm))
    }
    names(start_val) <- colnames(model_mx)
  }
  # Match functions for specified distributions of error terms
  ineff_distr_suffix <- paste(ineff_distr, collapse = "_")
  boundaries_fun <- match.fun(sprintf("get_variables_boundaries_%s", ineff_distr_suffix))
  neg_log_likelihood_fun <- match.fun(sprintf("neg_log_likelihood_%s", ineff_distr_suffix))
  estimate_extracted_surplus_fun <- match.fun(sprintf("estimate_extracted_surplus_%s", ineff_distr_suffix))

  # Calculate boundaries for parameters
  boundaries <- boundaries_fun(start_val = start_val,
                               error_param_min,
                               error_param_max)
  # Check if negative log-likelihood function returns finite values
  # for initial and boundary values of arguments
  check_variables_boundaries(boundaries = boundaries,
                             neg_ll = neg_log_likelihood_fun,
                             model_mx = model_mx,
                             y = y)
  par_names <- names(boundaries$start)

  # MLE estimation
  require(bbmle)
  parnames(neg_log_likelihood_fun) <- par_names
  mle_fit <- mle2(minuslogl = neg_log_likelihood_fun,
                 start = boundaries$start,
                 method = "L-BFGS-B",
                 data = list(model_mx = model_mx, y = y),
                 vecpar = TRUE,
                 parnames = par_names,
                 trace = (trace > 0),
                 lower = boundaries$lower,
                 upper = boundaries$upper,
                 control = list(parscale = get_parscale(boundaries$start),
                                maxit = maxit,
                                trace = trace))

  unexplained_variation <- sum(mle_fit@coef[c("sigma_v", "sigma_u", "sigma_w")]^2)
  unexplained_variation_from_inefficiency <- sum(mle_fit@coef[c("sigma_u", "sigma_w")]^2) / unexplained_variation
  surplus <- estimate_extracted_surplus_fun(variables = mle_fit@coef,
                                        formula = formula,
                                        data = data)
  return(named_list(mle_fit,
                    unexplained_variation,
                    unexplained_variation_from_inefficiency,
                    surplus,
                    ols_fit,
                    formula))
}
