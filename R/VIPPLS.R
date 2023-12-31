#' Calculation of variable importance in projections (VIP) scores
#' @name VIPPLS
#' @param pls_object a pls model generated by "pls" package
#' @param opt.comp the optimal component in the pls model
#' @param num_variables the number of variables in the pls model; extract from the dimension of coefficients in the pls model
#' @return numeric VIP scores
#' @examples
#' pls_object = plsr(y ~ ., ncomp = opt.comp, data = data)
#' opt.comp = 3
#' num_variables = dim(pls_object$coef)[1]
#' Vip <- VIP(pls_object, opt.comp = opt.comp)
#' @importFrom pls plsr
#' @export

# Filter method using variable importance in projections (VIP)
VIP <- function(pls_object, opt.comp, num_variables = dim(pls_object$coef)[1]) {
  # Calculate Variable Importance in Prediction (VIP) for PLS model
  # Extract loading weights from PLS model
  loading_weights <- pls_object$loading.weights

  # Create a weight matrix with squared loading weights normalized by the sum of squares of each variable
  weight_matrix <- sweep(loading_weights^2, 2, colSums(loading_weights^2), "/")

  # Extract Yloadings, score matrix T, and calculate Q^2
  Q <- pls_object$Yloadings
  T_matrix <- pls_object$scores
  Q2 <- rowSums(t(Q * Q))

  # Calculate SS(q_a t_a) = Q^2 * t_a^T * t_a based on the VIP equation
  Q2TT <- Q2[1:opt.comp] * diag(crossprod(T_matrix))[1:opt.comp]

  # Calculate VIP scores using the formula in research paper
  VIP_scores <- sqrt(num_variables * apply(sweep(weight_matrix[, 1:opt.comp, drop=FALSE],2,Q2TT,"*"), 1, sum) / sum(Q2TT))

  return(VIP_scores)
}
