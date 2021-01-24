option_pricing_milstein = function(S0, K, barrier, mu, sigma, T, N, n) {
  # N- number of trajectories
  # n- number of intervals (time discretization $dt=T/n$)
  # T- time horizon
  # S0- initial stock price
  # K- strike price
  # barrier- barrier
  # mu- "drift" in the pricing context, it should be equal to risk-free rate     
  # sigma- volatility   

  ###################
  
  X = GBM_milstein_method(y0 = S0, mu = mu, sigma = sigma, T = T, n = n, N = N)
  
  Call = exp(-mu*T) * mean(pmax(X[nrow(X), ] - K, 0))
  
  staterror_Call = exp(-mu*T) * sd(pmax(X[nrow(X), ] - K, 0)) / sqrt(N)
  d1 = (log(S0/K) + (mu + 1/2 * sigma^2) * T) / (sigma*sqrt(T))
  d2 = d1 - sigma * sqrt(T)
  Call_theoretical = S0 * pnorm(d1) - K * exp(-mu*T) * pnorm(d2)
  
  
  Put = exp(-mu*T) * mean(pmax(K - X[nrow(X), ], 0))
  staterror_Put = exp(-mu*T) * sd(pmax(K - X[nrow(X), ], 0)) / sqrt(N)
  Put_theoretical = -S0 * pnorm(-d1) + K * exp(-mu * T) * pnorm(-d2)
  
  CashorNothing_Call = exp(-mu*T) * mean(X[nrow(X), ] > K) # Payoff is 1 or 0
  staterror_CashorNothing_Call = exp(-mu*T) *  sd(X[nrow(X), ] > K) / sqrt(N)
  CashorNothing_Call_theoretical = exp(-mu * T) * pnorm(d2)
  
  CashorNothing_Put = exp(-mu*T) * mean(X[nrow(X),] < K) # Payoff is 1 or 0
  staterror_CashorNothing_Put = exp(-mu*T) * sd(X[nrow(X), ] < K) / sqrt(N)
  CashorNothing_Put_theoretical = exp(-mu * T) * pnorm(-d2)
  
  AssetorNothing_Call = exp(-mu*T) * mean(X[nrow(X), ] * (X[nrow(X), ] > K)) # Payoff is equal to $S_T$
  staterror_AssetorNothing_Call = exp(-mu*T) * sd(X[nrow(X), ] * (X[nrow(X), ] > K)) / sqrt(N)
  AssetorNothing_Call_theoretical = exp(-mu * T) * mean(X[nrow(X),]) * pnorm(d1)
  
  AssetorNothing_Put = exp(-mu*T) * mean(X[nrow(X), ] * (X[nrow(X),] < K)) # Payoff is equal to $S_T$
  staterror_AssetorNothing_Put = exp(-mu*T) * sd(X[nrow(X), ] * (X[nrow(X), ] < K)) / sqrt(N)
  AssetorNothing_Put_theoretical = exp(-mu * T) * mean(X[nrow(X),]) * pnorm(-d1)
  
  d3 = (log(S0/barrier) + (mu + 1/2 * sigma^2)*T)/(sigma * sqrt(T))
  d4 = d3 - sigma^2 * T
  d5 = (log(S0/barrier) - (mu - 1/2 * sigma^2)*T)/(sigma * sqrt(T))
  d6 = d5 - sigma^2 * T
  d7 = (log(S0*K / barrier^2) - (mu - 1/2 * sigma^2) * T)/(sigma * sqrt(T))
  d8 = d7 - sigma^2 * T
  
  a = (barrier / S0)^(-1 + 2 * mu / sigma^2)
  b = (barrier / S0)^(1 + 2 * mu / sigma^2)
  
  UpandIn_Call = exp(-mu*T) * mean(pmax(X[nrow(X), ] - K,0) * (apply(X, 2, max) > barrier))
  staterror_UpandIn_Call = exp(-mu*T) * sd(pmax(X[nrow(X), ] - K, 0) * (apply(X, 2, max) > barrier)) / sqrt(N)
  UpandIn_Call_theoretical = S0 * (pnorm(d3) + b * (pnorm(d6) - pnorm(d8))) -
    K * exp(-mu * T) * (pnorm(d4) + a * (pnorm(d5) - pnorm(d7)))
  
  UpandOut_Call = exp(-mu*T) * mean(pmax(X[nrow(X),] - K, 0) * (apply(X, 2, max) < barrier))
  staterror_UpandOut_Call = exp(-mu*T) * sd(pmax(X[nrow(X), ] - K,0) * (apply(X, 2, max) < barrier)) / sqrt(N)
  UpandOut_Call_theoretical = S0 * (pnorm(d1) - pnorm(d3) - b * (pnorm(d6) - pnorm(d8))) -
                              K * exp(-mu * T) * (pnorm(d2) - pnorm(d4) - a * (pnorm(d5) - pnorm(d7)))
  
  ####
  GeomAsian_Call = exp(-mu*T) * mean(pmax(apply(X, 2, function(i) {exp(mean(log(i)))}) - K, 0))
  staterror_GeomAsian_Call = exp(-mu*T) * sd(pmax(apply(X, 2, function(i) {exp(mean(log(i)))}) - K, 0))/sqrt(N)
  
  #sigma_A = sigma / sqrt(3)
  #b_A = 1/2 * (mu - sigma^2/6)
  
  #d1_A = (log(S0/K) + (b_A + sigma_A^2/2) * T) / (sigma_A * sqrt(T))
  #d2_A = d1 - sigma_A * sqrt(T)
  #GeomAsian_Call_theoretical = S0 * exp((b_A - mu) * T) * pnorm(d1_A) - K * exp(-mu * T) * pnorm(d2_A)
  
  GeomAsian_Put = exp(-mu*T) * mean(pmax(apply(X, 2, function(i) {K - exp(mean(log(i)))}), 0))
  staterror_GeomAsian_Put = exp(-mu*T) * sd(pmax(apply(X, 2, function(i) {K - exp(mean(log(i)))}), 0))/sqrt(N)
  #GeomAsian_Put_theoretical = K * exp(-mu * T) * pnorm(-d2_A) - S0 * exp((b_A - mu) * T) * pnorm(-d1_A)
  
  
  # description = c('Call','Put','CashorNothing_Call', 'CashorNothing_Put', 'AssetorNothing_Call', 'AssetorNothing_Put', 'UpandIn_Call', 'UpandOut_Call', 
  #                 "GeomAsian_Call", "GeomAsian_Put")
  # Prices = c(Call, Put, CashorNothing_Call, CashorNothing_Put, AssetorNothing_Call,
  #            AssetorNothing_Put,UpandIn_Call, UpandOut_Call, GeomAsian_Call, GeomAsian_Put)
  # stat_errors = c(staterror_Call, staterror_Put, staterror_CashorNothing_Call, staterror_CashorNothing_Put, staterror_AssetorNothing_Call,
  #                 staterror_AssetorNothing_Put,staterror_UpandIn_Call, staterror_UpandOut_Call,
  #                 staterror_GeomAsian_Call, staterror_GeomAsian_Put)
  # result = rbind(description, Prices, stat_errors)
  # result
  data.frame(type = c('Call','Put','CashorNothing_Call', 'CashorNothing_Put', 'AssetorNothing_Call', 'AssetorNothing_Put', 'UpandIn_Call', 'UpandOut_Call', 
                      "GeomAsian_Call", "GeomAsian_Put"),
             Price = c(Call, Put, CashorNothing_Call, CashorNothing_Put, AssetorNothing_Call,
                       AssetorNothing_Put,UpandIn_Call, UpandOut_Call, GeomAsian_Call, GeomAsian_Put),
             Theoretical = c(Call_theoretical, Put_theoretical, CashorNothing_Call_theoretical, CashorNothing_Put_theoretical,
                             AssetorNothing_Call_theoretical, AssetorNothing_Put_theoretical,
                             UpandIn_Call_theoretical, UpandOut_Call_theoretical, rep(0, 2)),
             stat_errors = c(staterror_Call, staterror_Put, staterror_CashorNothing_Call, staterror_CashorNothing_Put, staterror_AssetorNothing_Call,
                             staterror_AssetorNothing_Put,staterror_UpandIn_Call, staterror_UpandOut_Call,
                             staterror_GeomAsian_Call, staterror_GeomAsian_Put)
             )
}
