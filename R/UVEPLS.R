# Wrapped method - UVE
#' Uninformative Variable Elimination in PLS
#'
#' This function performs Multivariate Calibration by Uninformative Variable Elimination on the given datasets.
#' It combines the original dataset with a noise matrix, then iteratively fits a PLSR model, selecting the most informative variables and refining the model.
#' @name UVEPLS
#' @param X The predictor matrix (numeric). If not a matrix or data.frame, it will be converted to a matrix.
#' @param Y The response matrix (numeric).
#' @param ncomp The number of components to use in the PLSR model. Default is 10.
#' @param s The number of iterations for variable elimination. Default is 30.
#' @param split The proportion of data to be used for training in each iteration. Default is 0.8.
#' @param scale Logical; if TRUE, the predictor matrix X will be scaled. Default is TRUE.
#' @return A list containing several components of the final PLSR model: scores, loadings, Yscores, Yloadings, loading.weights, residuals, and fitted.values.
#' @examples
#' # Example usage:
#' # mcuve(X = matrix(rnorm(100), 10, 10), Y = matrix(rnorm(10), 10, 1), ncomp = 15)
#' @export

mcuve <- function(X, Y, ncomp = 10, s = 30, split = 0.8, scale = TRUE){
  if (!is.matrix(X) && !is.data.frame(X)) {
    X <- as.matrix(X)
  }
  Y <- as.matrix(Y)

  # Optional scaling to normalize X
  if (scale) {
    X <- scale(X)
  }

  n <- nrow(X)
  p <- ncol(X)
  # Generate noise matrix N
  N  <- matrix(runif(n*p, min = 0, max = 1), nrow = n, ncol = p)
  # Combine N and X matrix in new matrix Z
  Z  <- cbind(X, N)

  C_j  <- matrix(0, nrow = s, ncol = p*2)

  best_performance <- Inf
  #stop_criterion <- FALSE
  count <- 0

  for(i in 1:s){
    # Split the data set to train data according to preferred proportion
    train <- sample(n,round(n*split))
    Y_train <- Y[train,]
    Z_train <- Z[train, ]
    ncomp <- min(n, ncomp)

    # Fit PLSR model to combined matrix Z_train
    pls_model <- plsr(Y_train ~ Z_train,  ncomp = ncomp, validation = "LOO")
    # Find the optimal number of components
    optimal_number_of_components <- which.min(pls_model$validation$PRESS[1,])
    C_j[i,]  <- coef(pls_model, ncomp = optimal_number_of_components)

    # Compute test statistic for each variable
    vars_statistic <- apply(C_j, 2, function(x) mean(x) / sd(x))
    # Set threshold and eliminate noise variables
    c_max <- max(abs(vars_statistic[-(1:p)]))
    # Update the selected variables indices - current_vars_indices
    informative_vars_indices <- which(abs(vars_statistic[(1:p)]) >= c_max)

    # Test data set
    Z_test <- Z[-train, ]
    Y_test <- Y[-train,]

    # Fit the model and predict with test data
    pls_test <- plsr(Y_train ~ Z_train,  ncomp = optimal_number_of_components)
    predictions <- predict(pls_test, newdata = Z_test)

    # Calculate RMSE for Model Performance
    rmse <- sqrt(mean((Y_test - predictions)^2))

    # Check performance and decide whether stop or not
    if(rmse < best_performance || i <= 20) { # set a least MC simulations to ensure stability
      best_performance <- rmse # better performance
      count <- count + 1
    } else {
      #stop_criterion <- TRUE
      break  # Stop if performance decreases
    }
  }

  # Update X for the final model
  X <- as.matrix(X[, informative_vars_indices])

  # Fit the final model with selected variables
  pls_final_model <- plsr(Y ~ X,  ncomp = ncomp)
  return(list("model"=pls_final_model,"index" = informative_vars_indices))
}
