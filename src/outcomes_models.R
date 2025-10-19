#' Estimate outcome models for the alters
#' Obtain predicted conditional expectations for the alters (F=1 and F=0)
#' Also, obtain predicted conditional expectations for the egos (Z=0, F=1)
#'
#'
#' @param Y_a A numeric vector of outcomes for the main dataset (alters).
#'        Length is n_a.
#' @param X_a A numeric matrix of covariates for the main dataset, with n_a rows
#'        and p columns.
#' @param F_a A binary numeric vector (0 or 1) representing the observed
#'        exposure of alters. Length is n_a.
#' @param X_e A numeric matrix of covariates for a second dataset (egos),
#'        often a subset or different population, with n_e rows and p columns.
#'        The column space must match that of X_a.
#' @param reg_model A function that fits a regression model, e.g., `glm`, `lm`,
#'        `gee::gee`. The function must accept a `formula` and `data` argument
#'        and return an object compatible with the generic `predict()` function.
#' @param formula A formula object specifying the regression model to be fit.
#'        The variable names should correspond to the columns created within the
#'        function: "Y" for the outcome, "F" for the exposures (of alters),
#'         and "X1", "X2", ..., "Xp" for the covariates.
#' @param ... Additional arguments to be passed to the `reg_model` function
#'        (e.g., `family = binomial` for `glm`).
#'
#' @return A list containing three numeric vectors of predictions:
#' \describe{
#'   \item{mu_1_alters}{Predictions for the alters (`X_a`) assuming everyone
#'         exposed to the treatment (F=1).}
#'   \item{mu_0_alters}{Predictions for the alters (`X_a`) assuming everyone
#'         exposed to the control (F=0).}
#'   \item{mu_01_egos}{Predictions for the egos (`X_e`) assuming everyone
#'         exposed to the treatment (F=1), but not treated directly (Z=0).}
#' }
#'
#'
#' @examples
#' # 1. Generate synthetic data
#' set.seed(42)
#' n_a <- 200
#' n_e <- 100
#' p <- 3
#'
#' X_a <- matrix(rnorm(n_a * p), ncol = p)
#' X_e <- matrix(rnorm(n_e * p), ncol = p)
#' F_a <- rbinom(n_a, 1, 0.5)
#'
#' # True relationship: Y depends on F, X1, X2, X3
#' true_beta <- c(1.5, -1, 0.5) # X1, X2, X3
#' linear_pred <- 1 + 0.5 * F_a + X_a %*% true_beta
#' prob <- plogis(linear_pred) # Inverse logit
#' Y_a <- rbinom(n_a, 1, prob)
#'
#' # 2. Define the formula
#' # Note: Covariate names are X1, X2, X3 as created by the function
#' f <- as.formula(Y ~ F + X1 + X2 + X3)
#'
#' # 3. Run the function with glm
#' predictions <- reg_func(
#'   Y_a = Y_a,
#'   X_a = X_a,
#'   F_a = F_a,
#'   X_e = X_e,
#'   reg_model = glm,
#'   formula = f,
#'   family = binomial(link = "logit") # Additional arg for glm
#' )
#'

reg_func_alters <- function(Y_a, X_a, F_a, X_e, reg_model, formula, ...) {
  
  # --- Input Validation ---
  if (!is.vector(Y_a) || !is.numeric(Y_a)) stop("Y_a must be a numeric vector.")
  if (!is.matrix(X_a) || !is.numeric(X_a)) stop("X_a must be a numeric matrix.")
  if (!is.vector(F_a) || !is.numeric(F_a)) stop("F_a must be a numeric vector.")
  if (!all(F_a %in% c(0, 1))) stop("F_a must be a binary vector (0 or 1).")
  if (!is.matrix(X_e) || !is.numeric(X_e)) stop("X_e must be a numeric matrix.")
  if (ncol(X_a) != ncol(X_e)) stop("X_a and X_e must have the same number of columns (p).")
  if (length(Y_a) != nrow(X_a) || length(Y_a) != length(F_a)) {
    stop("Y_a, X_a (rows), and F_a must have the same length (n_a).")
  }
  if (!is.function(reg_model)) stop("reg_model must be a function.")
  if (!inherits(formula, "formula")) stop("formula must be a formula object.")
  
  
  # --- 1. Prepare Data and Fit Model ---
  
  # Convert matrices to data frames for compatibility with formula interface
  X_a_df <- as.data.frame(X_a)
  p <- ncol(X_a_df)
  # Assign generic, predictable names for use in the formula
  colnames(X_a_df) <- paste0("X", 1:p)
  
  # Combine into a single data frame for the regression model
  dat_a <- data.frame(Y = Y_a, F = F_a)
  dat_a <- cbind(dat_a, X_a_df)
  
  # Use do.call to flexibly call the provided regression function
  # with the formula, data, and any additional arguments (...)
  fit <- do.call(reg_model, list(formula = formula, data = dat_a, ...))

  # --- 2. Create Counterfactual Datasets for Prediction ---
  
  # New data for alters where everyone are exposed (F=1)
  newdata_a_f1 <- dat_a
  newdata_a_f1$F <- 1
  
  # New data for alters where none are exposed (F=0)
  newdata_a_f0 <- dat_a
  newdata_a_f0$F <- 0
  
  # Prepare the 'egos' data for prediction. It needs the same column names
  # as the training data.
  X_e_df <- as.data.frame(X_e)
  colnames(X_e_df) <- paste0("X", 1:p)
  
  # New data for egos where everyone is treated (F=1)
  newdata_e_f1 <- X_e_df
  newdata_e_f1$F <- 1
  
  # --- 3. Generate Predictions ---
  
  # Use the generic predict() function with 'type = "response"' for GLM
  mu_1_alters <- predict(fit, newdata = newdata_a_f1, type = "response")
  mu_0_alters <- predict(fit, newdata = newdata_a_f0, type = "response")
  mu_01_egos <- predict(fit, newdata = newdata_e_f1, type = "response")

  # --- 4. Return Results ---
  
  # Return a named list of vectors, ensuring output is clean.
  results <- list(
    mu_1_alters = as.vector(mu_1_alters),
    mu_0_alters = as.vector(mu_0_alters),
    mu_01_egos  = as.vector(mu_01_egos)
  )
  
  return(results)
}

reg_func_egos <- function(Y_e, X_e, Z_e, reg_model, formula, ...) {
  
  # --- Input Validation ---
  if (!is.vector(Y_e) || !is.numeric(Y_e)) stop("Y_e must be a numeric vector.")
  if (!is.matrix(X_e) || !is.numeric(X_e)) stop("X_e must be a numeric matrix.")
  if (!is.vector(Z_e) || !is.numeric(Z_e)) stop("Z_e must be a numeric vector.")
  if (!all(Z_e %in% c(0, 1))) stop("Z_e must be a binary vector (0 or 1).")
  if (length(Y_e) != nrow(X_e) || length(Y_e) != length(Z_e)) {
    stop("Y_e, X_e (rows), and Z_e must have the same length (n_e).")
  }
  if (!is.function(reg_model)) stop("reg_model must be a function.")
  if (!inherits(formula, "formula")) stop("formula must be a formula object.")
  
  
  # --- 1. Prepare Data and Fit Model ---
  
  # Convert matrices to data frames for compatibility with formula interface
  X_e_df <- as.data.frame(X_e)
  p <- ncol(X_e_df)
  # Assign generic, predictable names for use in the formula
  colnames(X_e_df) <- paste0("X", 1:p)
  
  # Combine into a single data frame for the regression model
  dat_e <- data.frame(Y = Y_e, Z = Z_e)
  dat_e <- cbind(dat_e, X_e_df)
  
  # Use do.call to flexibly call the provided regression function
  # with the formula, data, and any additional arguments (...)
  fit <- do.call(reg_model, list(formula = formula, data = dat_e, ...))

  # --- 2. Create Counterfactual Datasets for Prediction ---
  
  # New data for egos where everyone is treated (Z=1)
  newdata_e_z1 <- dat_e
  newdata_e_z1$Z <- 1
  
  # New data for egos where no one is treated (Z=0)
  newdata_e_z0 <- dat_e
  newdata_e_z0$Z <- 0
  
  # --- 3. Generate Predictions ---
  
  # Use the generic predict() function with 'type = "response"' for GLMs
  mu_10_e <- predict(fit, newdata = newdata_e_z1, type = "response")
  mu_00_e <- predict(fit, newdata = newdata_e_z0, type = "response")

  # --- 4. Return Results ---
  
  # Return a named list of vectors, ensuring output is clean.
  results <- list(
    mu_10_e = as.vector(mu_10_e),
    mu_00_e = as.vector(mu_00_e)
  )
  
  return(results)
}
