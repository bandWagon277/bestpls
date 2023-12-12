#' A basic implementation of BVE-PLS (Backward variable elimination)
#' @name BVEPLS
#' @param y description a n*q matrix
#' @param X description a n*p matrix
#' @return a object of class
#' @examples
#' y=matrix(rnorm(100),5,20)
#' X=matrix(rnorm(200),10,20)
#' fit = bve_pls(y,X)
#' @importFrom pls plsr
#' @export

# Load necessary package
library(pls)

#NAMESPACE FILE
#source("VIPPLS.R")

# Wrapper method by backward variable elimination in PLS
bve_pls <- function(y, X, ncomp = 10, test_ratio = 0.25, VIP.threshold = 2) {

  # Combine predictors and response into a data frame
  data_df <- data.frame(y = y, X)

  # Determine the number of samples
  n <- nrow(data_df)
  p <- ncol(X)

  # Set the seed for reproducibility
  set.seed(123)

  # Generate indices for the test set
  test_indices <- sample(1:n, size = round(test_ratio * n))

  # Create the test set
  test_set <- data_df[test_indices, ]
  test_set_X <- X[test_indices, ]
  test_set_y <- y[test_indices, ]

  # Create the training set
  train_set <- data_df[-test_indices, ]
  train_set_X <- X[-test_indices, ]
  train_set_y <- y[-test_indices, ]

  # Create initial variables
  terminated <- FALSE
  rmsep_values <- numeric(0)
  Variable.list <- list()
  is.selected <- rep(TRUE, p)
  variables <- which(is.selected)

  while (!terminated) {
    # Fitting, and finding optimal number of components
    co <- min(ncomp, ncol(train_set_X) - 1)

    # Using plsr from the pls package
    pls.object <- plsr(y ~ ., ncomp = co, data = data.frame(y = I(train_set_y), train_set_X), validation = "LOO")
    opt.comp <- which.min(pls.object$validation$PRESS[1,])

    # Using train set to fit pls model
    pls_fit_model <- plsr(y ~ ., ncomp = opt.comp, data = train_set)

    # Predict on the test set
    predict_pls <- predict(pls_fit_model, ncomp = opt.comp, newdata = test_set)

    # Calculate the prediction error using root-mean-square error of prediction (RMSEP)
    rmsep_values_i <- sqrt(mean((test_set_y - predict_pls)^2))
    rmsep_values <- c(rmsep_values, rmsep_values_i)

    # Calculate VIP scores
    Vip <- VIP(pls_fit_model, opt.comp = opt.comp)
    VIP.index <- which(as.matrix(Vip) < VIP.threshold)

    if (length(VIP.index) <= (ncomp + 1)) {
      VIP.index <- sort(Vip, decreasing = FALSE, index.return = TRUE)$ix[1:ncomp]
    }

    # Update the selected variables
    is.selected[variables[VIP.index]] <- FALSE
    variables <- which(is.selected)
    Variable.list <- c(Variable.list, list(variables))

    # Remove constant variables with zero variance
    variables_remove <- unique(which(apply(train_set_X, 2, var) == 0), which(apply(test_set_X, 2, var) == 0))
    train_set_X <- train_set_X[, -variables_remove]
    test_set_X <- test_set_X[, -variables_remove]

    if (ncol(train_set_X) <= ncomp + 1) {
      # Terminate if fewer than ncomp + 1 variables remain (set termination condition)
      terminated <- TRUE
    }
  }

  # Find the iteration with the minimum prediction error as the optimal iteration
  opt.iter <- which.min(rmsep_values)
  bve_selection <- Variable.list[[opt.iter]]

  # Create a subset of X with the selected variables
  X_sub <- X[bve_selection]

  # Combine response y with those selected variable X_sub as a data frame
  data_df_new <- data.frame(y = y, X_sub)

  # Using new data frame to fit pls model as the final model
  pls_fit_model_final <- plsr(y ~ ., ncomp = opt.comp, data = data_df_new)

  return(pls_fit_model_final)
}
