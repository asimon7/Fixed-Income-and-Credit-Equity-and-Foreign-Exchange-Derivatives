GBM_milstein_method = function(y0, mu, sigma, T, n, N) {
  # y0- initial value
  # mu- drift
  # sigma- volatility
  # T- time horizon
  # n- number of intervals
  # N- number of trajectories
  
  dt = T / n
  X = matrix(0, n + 1, N)
  dW = matrix(rnorm(n*N, 0, sqrt(dt)), n, N)
  X[1, ] = y0
  for (i in 1:n) {
    X[i + 1, ] = X[i, ] + mu * X[i, ] * dt + sigma * X[i, ] * dW[i, ] +
      0.5 * sigma^2 * X[i, ] * (dW[i, ]^2 - dt)
  }
  
  X
}