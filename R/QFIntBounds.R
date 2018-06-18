QFIntBounds <- function(obs, M, mu, sigma, k = c(20), resid.op.norm.bd = NULL, lower.tail = TRUE, log = FALSE) {

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

  # Eigen-decompose and rescale M by N in order to keep magnitude of observations in a more maneagable range
  M.tilde <- sweep(M * sigma, 2, sigma, "*")/N
  mu.tilde <- mu/sigma
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
