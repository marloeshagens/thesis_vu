### ### ### ### ### RANDOM WALK FUNCTION ### ### ### ### ###

# data generating function for a random walk process Y
# number of periods (T)
# initial value (x0)
# drift (mu) = delta
# variance of innovation

# model: y[t] = delta + y[t-1] + u[t]
randomwalk <- function(T, y0, delta, variance) {
  u <- cumsum(rnorm(n=T, mean=0, sd=sqrt(variance))) # normal randomly distributed  error terms
  t <- 1:T
  y <- y0 + t*delta + u
  return(y)
}