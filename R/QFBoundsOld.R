#' Quadratic Form Bounds
#'
#' Compute upper and lower bounds on the CDF of a quadratic form in normal random variables.
#'
#' The above was the one line summary which appears at the top of the documentation.  The gory detailed description of what is does goes here.
#'
#' Those gory details can span as many paragraphs as you want.
#'
#' @param obs vector; observed values of the quadratic form for which upper and lower bounds on the CDF
#' @param evals vector of eigen-values of the matrix \eqn{M}.  These need not be all eigen values.
#' @param ncps description here
#' @param E_R description here
#' @param nu description here
#' @param N description here
#' @param resid.op.norm.bd description here
#' @param if.insuff.eigs string indicating what action to take if there are insufficient eigen-values to produce an accurate bound.  If \code{"trivial"} then the bounds \eqn{[0,1]} are returned; if \code{"missing"} then \code{NA} is returned for both bounds.
#' @param lower.tail logical; if \code{TRUE} (default), probability is \eqn{P(y^T M y \le obs)}, otherwise \eqn{P(y^T M y > obs)}
#' @param log logical; if \code{TRUE}, probability \eqn{p} is given as \eqn{log(p)}
#'
#' @return A data frame containing the variables lower and upper which provide the bounds on the CDF of the quadratic form.
#'
#' @examples
#' # Some code here which runs a self-contained example
#'

QFBounds <- function(obs, M, mu, sigma, k = c(20), resid.op.norm.bd = NULL, lower.tail = TRUE, log = FALSE) {

  # All input must be on the original scale.  obs= y^T M y

  if(class(M) != "dsyMatrix") {
    stop("Matrix M must be symmetric of class dsyMatrix (from Matrix package).")
  }
  if(!all(diag(M) == 0)) {
    stop("M must be a hollow matrix (zero on diagonal).")
  }
  if(!is.vector(obs)) {
    stop("obs must be a vector.")
  }
  # obs <- suppressMessages(as.numeric(crossprod(crossprod(M, y), y)))
  N <- nrow(M)
  if(!is.vector(sigma) || length(sigma) != N) {
    stop("sigma must be a vector matching the number of rows of M.")
  }
  if(!is.vector(mu) || length(mu) != N) {
    stop("mu must be a vector matching the number of rows of M.")
  }

  # Eigen-decompose
  mu.tilde <- mu/sigma
  M.tilde <- sweep(M * sigma, 2, sigma, "*")
  e <- eigs(M.tilde, max(k), which = "LM")
  evec.tilde <- e$vectors[, order(abs(e$values), decreasing=TRUE)]
  eval.tilde <- e$values[order(abs(e$values), decreasing=TRUE)]

  # Set resid.op.norm.bd to magnitude of smallest truncated eigenvalue if not specified
  if(is.null(resid.op.norm.bd)) {
    resid.op.norm.bd <- abs(eval.tilde[length(eval.tilde)])
  }

  # Compute ncps
  ncps <- list()

  for(kk in 1:length(k)) {
    ncps[[kk]] <- c(crossprod(evec.tilde[,1:k[kk]], mu.tilde))^2
  }

  # nu and E
  nu <- list()
  E_R <- list()
  for(kk in 1:length(k)) {
    R <- (M.tilde - evec.tilde[,1:k[kk]] %*% (t(evec.tilde[,1:k[kk]]) * eval.tilde[1:k[kk]]))
    nu[[kk]] <- 8 * sum((R %*% mu.tilde)^2) + 4 * sum(R^2)
    E_R[[kk]] <- as.numeric(mu.tilde %*% R %*% mu.tilde + sum(diag(R)))
  }

  # Call complex function
  res <- data.frame(k = c(), obs = c(), lower = c(), upper = c())
  for(kk in 1:length(k)) {
    for(i in 1:length(obs)) {
      res <- rbind(res, c(k = k[kk],
                          obs = obs[i],
                          QFBounds2(obs[i], eval.tilde[1:k[kk]], ncps[[kk]], E_R[[kk]], nu[[kk]], N, resid.op.norm.bd, "missing", lower.tail, log)))
    }
  }
  res
}
