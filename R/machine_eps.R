machine_eps <- function(value) {
  .C(C_machine_eps, res = as.double(value))$res
}
