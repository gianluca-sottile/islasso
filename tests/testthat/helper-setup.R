sim_gaus <- function(n = 80, p = 15, seed = 1, rho = 0) {
  simulXy(n = n, p = p, beta = c(1, -1, 0.5, rep(0, p - 3)), seed = seed, rho = rho)
}
# numerical gradient of the penalised objective (standardised scale)
num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) { e <- replace(numeric(length(x)), j, h); (f(x + e) - f(x - e)) / (2 * h) }, 0)
}
