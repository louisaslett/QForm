library("Rmpfr")
library("glue")

acc <- 128
twiddle.table <- list()
for(pwr in c(5,16)) {
  N <- mpfr(2^pwr, acc)
  twiddle.table[[pwr]] <- matrix(c(mpfr(0:(2^pwr-1), acc),
                                   -2*Const("pi", acc)*mpfr(0:(2^pwr-1), acc)/N),
                                 nrow = as.integer(N), ncol = 2)
  twiddle.table[[pwr]][,1] <- cos(twiddle.table[[pwr]][,2])
  twiddle.table[[pwr]][,2] <- sin(twiddle.table[[pwr]][,2])
}

c.atan2 <- function(y, x) {
  res <- 2*atan(y/(x+sqrt(x^2+y^2)))
  if(any(y==0))
    res[y==0] <- Const("pi", acc)
  res
}

c.log <- function(z) {
  matrix(c(log(sqrt(z[,1]^2+z[,2]^2)),
           c.atan2(z[,2], z[,1])), nrow = nrow(z), ncol = 2)
}

bitrevorder <- function(i) {
  nonzeros <- max(which(intToBits(i)>0))-1
  sapply(0:(i-1), function(x, nonzeros) {
    packBits(c(rev(intToBits(x)[1:nonzeros]), intToBits(0)[1:(32-nonzeros)]), "integer")
  }, nonzeros = nonzeros) + 1
}

c.mul <- function(x, y) {
  matrix(c(
    x[,1]*y[,1] - x[,2]*y[,2],
    x[,2]*y[,1] + x[,1]*y[,2]
  ), nrow = nrow(x), ncol = 2)
}

# Arbitrary precision FFT
fft.ap <- function(z) {
  N <- length(z)
  n.stages <- log2(N)
  if(!(n.stages %in% which(lapply(twiddle.table, length)>0))) {
    stop("z must be 2^16 thru 2^19 long")
  }

  tt <- twiddle.table[[n.stages]]

  new.i <- bitrevorder(N)
  x <- mpfr(matrix(c(Re(z)[new.i], Im(z)[new.i]), ncol = 2), acc)
  x.res <- mpfr(matrix(0, nrow = nrow(x), ncol = ncol(x)), acc)

  for(p in 1:n.stages) {
    alpha <- 2^(p-1)
    bf.start <- 1
    while(bf.start <= (N-alpha)) {
      for(ibf.idx in 0:(alpha-1)) {
        x.res[bf.start,]         <- x[bf.start,, drop = FALSE] + c.mul(x[bf.start+alpha,, drop = FALSE], tt[ibf.idx*2^(n.stages - p) + 1,, drop = FALSE])
        x.res[bf.start + alpha,] <- x[bf.start,, drop = FALSE] - c.mul(x[bf.start+alpha,, drop = FALSE], tt[ibf.idx*2^(n.stages - p) + 1,, drop = FALSE])
        # cat(glue("bf.start = {bf.start}, alpha = {alpha}\n{as.double(x.res[bf.start,])}\n{as.double(x.res[bf.start + alpha,])} \n \n"))
        bf.start <- bf.start + 1
        if(ibf.idx == alpha-1) {
          bf.start <- bf.start + alpha
        }
      }
    }
    x <- x.res
  }
  x
}

fft.ap2 <- function(z) {
  N <- nrow(z)
  n.stages <- log2(N)
  if(!(n.stages %in% which(lapply(twiddle.table, length)>0))) {
    stop("z must be 2^16 thru 2^19 long")
  }

  tt <- twiddle.table[[n.stages]]

  new.i <- bitrevorder(N)
  x <- z
  x.res <- mpfr(matrix(0, nrow = nrow(x), ncol = ncol(x)), acc)

  for(p in 1:n.stages) {
    alpha <- 2^(p-1)
    bf.start <- 1
    while(bf.start <= (N-alpha)) {
      for(ibf.idx in 0:(alpha-1)) {
        x.res[bf.start,]         <- x[bf.start,, drop = FALSE] + c.mul(x[bf.start+alpha,, drop = FALSE], tt[ibf.idx*2^(n.stages - p) + 1,, drop = FALSE])
        x.res[bf.start + alpha,] <- x[bf.start,, drop = FALSE] - c.mul(x[bf.start+alpha,, drop = FALSE], tt[ibf.idx*2^(n.stages - p) + 1,, drop = FALSE])
        # cat(glue("bf.start = {bf.start}, alpha = {alpha}\n{as.double(x.res[bf.start,])}\n{as.double(x.res[bf.start + alpha,])} \n \n"))
        bf.start <- bf.start + 1
        if(ibf.idx == alpha-1) {
          bf.start <- bf.start + alpha
        }
      }
    }
    x <- x.res
  }
  x
}
