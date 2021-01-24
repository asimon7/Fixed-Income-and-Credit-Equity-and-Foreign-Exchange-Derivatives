library(dplyr)
library(tidyr)
library(ggplot2)

Convergence = function(S0, mu, sigma, t, n, N) {
  dt = T / n
  Brownian_M = matrix(0, n+1, N)
  dW = matrix(rnorm(n*N, 0, sqrt(dt)), n, N)
  Brownian_M[1, ] = S0
  Euler_M = Brownian_M
  Milstein_M = Brownian_M
  dW = matrix(rnorm(n*N, 0, 1), n, N) * sqrt(dt)

  for (i in 1:n) {
    Brownian_M[i + 1,] = Brownian_M[i, ] * exp((mu - sigma^2/2) * dt + dW[i, ] * sigma)
    
    Euler_M[i + 1, ] = Euler_M[i, ] + mu * Euler_M[i, ] * dt + sigma * Euler_M[i,] * dW[i,]
    
    Milstein_M[i + 1, ] = Milstein_M[i, ] + mu * Milstein_M[i, ] * dt + sigma * Milstein_M[i, ] * dW[i, ] +
      0.5 * sigma^2 * Milstein_M[i, ] * (dW[i, ]^2 - dt)
  }
  
  err_em = apply(abs(Brownian_M - Euler_M), 1, sum)
  err_mil = apply(abs(Brownian_M - Milstein_M), 1, sum)
  
  data.frame(
    n = n,
    Euler_strong_error = max(err_em) / N,
    Milstein_strong_error = max(err_mil) / N
  )
}

S0 = 100
mu = 1
sigma = 0.2
t = 1
N = 1000

n_vec = c(2^(10:15))


lapply(n_vec, Convergence, S0 = S0, mu = mu, sigma = sigma, t = t, N = N) %>%
  bind_rows() %>%
  gather(key, value, -n) %>%
  ggplot(aes(x = n, y = value, group = key)) +
  geom_point(aes(color = key))+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")

a = lapply(n_vec, Convergence, S0 = S0, mu = mu, sigma = sigma, t = t, N = N) %>%
  bind_rows()
a

lm(log(Euler_strong_error) ~ log(n), a)

lm(log(Milstein_strong_error) ~ log(n), a)
