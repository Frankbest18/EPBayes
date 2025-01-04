library(REBayes)
library(Rmosek)
library(asht)

# P value calculator for 1D Parametric MLE

p_tau_j_lambda = function (lambda, par, tau_j, n1, n2) {
  k = par[1]
  d1 = par[2]
  d2 = par[3]
  tau_j = tau_j
  out = 1/beta((n1-1)/2, (n2-1)/2) * ((n1-1)/(n2-1))^((n1-1)/2) * lambda^((n2-1)/2) * (lambda + (n1-1)/(n2-1) *tau_j)^(-(n1+n2-2)/2) * 1/beta(d1/2, d2/2) * lambda^(d1/2 - 1) * (lambda + (d2 * k)/d1)^(-(d1 + d2)/2)
  #out = 1/beta((n1-1)/2, (n2-1)/2) * ((n1-1)/(n2-1))^((n1-1)/2) * tau_j^((n1-3)/2) * lambda^((n2-1)/2) * (lambda + (n1-1)/(n2-1) *tau_j)^(-(n1+n2-2)/2) * 1/beta(d1/2, d2/2) * (d1/(d2 * k))^(-d2/2) * lambda^(d1/2 - 1) * (lambda + (d2 * k)/d1)^(-(d1 + d2)/2)
}

logL = function(par, tau, n1, n2) {
  logliklihood = 0
  k = par[1]
  d1 = par[2]
  d2 = par[3]
  for (tau_j in tau) {
    #p_tau_j = tryCatch(integrate(p_tau_j_lambda, lower = 0, upper = Inf, par = par, tau_j = tau_j, n1 = n1, n2 = n2)$value, error = function(e) {print(par)
      #stop(tau_j)}
    
    p_tau_j = tryCatch(tau_j^((n1-3)/2) * (d1/(d2 * k))^(-d2/2) * integrate(p_tau_j_lambda, lower = 0, upper = Inf, par = par, tau_j = tau_j, n1 = n1, n2 = n2)$value, error = function(e) {print(par)
      stop(tau_j)})
    logliklihood = logliklihood + log(p_tau_j)
  }
  
  -logliklihood
}

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

p_given_lambda = function(n1, n2, Z1, Z2, S1, S2, lambda) {
  
  upper = (Z1 - Z2) * sqrt(n1 + n2 - 2)
  lower = sqrt((n1-1) * (n2 * lambda + n1) / (n1 * n2 * lambda) * S1 + (n2-1) * (n2 * lambda + n1) / (n1 * n2) * S2)
  
  test_stat = upper/lower
  
  p = pt(q = abs(test_stat), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  
  out = p
}

p_joint_parametric_j = function(lambda, n1, n2, Z1, Z2, S1, S2, par_mle, tau_j) {
  k = par_mle[1]
  d1 = par_mle[2]
  d2 = par_mle[3]
  
  p_value_given_lambda = p_given_lambda(n1, n2, Z1, Z2, S1, S2, lambda)
  
  g_lambda_mle = 1 / beta(d1/2, d2/2) * (d1/(d2*k))^(-d2/2) * lambda^(d1/2 - 1) * (lambda + (k * d2)/d1)^(-(d1+d2)/2)
  
  f_tau_j_given_lambda = p_tau_j_given_lambda(n1, n2, tau_j, lambda)
  
  #f_tau_j = integrate(p_tau_j_lambda, lower = 0, upper = Inf, par = par_mle, tau_j = tau_j, n1 = n1, n2 = n2)$value
  
  f_tau_j = tau_j^((n1-3)/2) * (d1/(d2 * k))^(-d2/2) * integrate(p_tau_j_lambda, lower = 0, upper = Inf, par = par_mle, tau_j = tau_j, n1 = n1, n2 = n2)$value
  
  out = p_value_given_lambda * f_tau_j_given_lambda * g_lambda_mle / f_tau_j
  
}

P_value_1D_PMLE = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list) {
  
  tau = S1_list/S2_list
  
  par_mle = optim(c(5, 5, 5), logL, tau = tau, n1 = n1, n2 = n2, lower = c(0.01, 0.01, 0.01), upper = c(100, 100, 100), method = 'L-BFGS-B')$par
  
  #print(par_mle)
  #print(logL(par_mle, tau, n1, n2))
  
  P_value_list_pmle = rep(0, m)
  
  for (i in c(1:m)) {
    P_value_list_pmle[i] = integrate(p_joint_parametric_j, lower = 0, upper = Inf, n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], par_mle = par_mle, tau_j = tau[i])$value
  }
  
  return (P_value_list_pmle)
}

# P value calculator for 1D Non-Parametric MLE

mass_given_tau_j = function(n1, n2, grid, mass, tau_j) {
  f_tau_j_given_grid = p_tau_j_given_lambda(n1, n2, tau_j, grid)
  f_tau_j_grid = f_tau_j_given_grid * mass
  post_mass = f_tau_j_grid / sum(f_tau_j_grid)
  return(post_mass)
}

p_value_npmle_j = function(n1, n2, Z1, Z2, S1, S2, grid, mass, tau_j) {
  P_value_joint = p_given_lambda(n1, n2, Z1, Z2, S1, S2, grid)
  post_mass_j = mass_given_tau_j(n1, n2, grid, mass, tau_j)
  P_value_j = sum(P_value_joint * post_mass_j)
  return (P_value_j)
}

NPMLE_1D = function(tau, B, m, n1, n2, lower_quantile, upper_quantile) {
  lower = quantile(tau, lower_quantile)
  upper = quantile(tau, upper_quantile)
  log_u = seq(log(lower), log(upper), length = B)
  u = exp(log_u)
  d = rep(1,B)
  w = rep(1, m)/m
  A = outer(tau, u, FUN = p_tau_j_given_lambda, n1 = n1, n2 = n2)
  result = KWPrimal(A, d, w)
  mass = result$f/sum(result$f)
  
  output = list('grid' = u, 'mass' = mass)
  
}

P_value_1D_NPMLE = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, NPMLE_par) {
  
  tau = S1_list/S2_list
  
  NPMLE_result = NPMLE_1D(tau = tau, m = m, n1 = n1, n2 = n2, B = NPMLE_par[1], lower_quantile = NPMLE_par[2], upper_quantile = NPMLE_par[3])
  
  P_value_list_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle[i] = p_value_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], NPMLE_result$grid, NPMLE_result$mass, tau[i])
  }
  
  return(P_value_list_npmle)
}

# P value calculator for 2D Non-Parametric MLE

p_s_j_given_sigma2 = function(n1, n2, s1_j, s2_j, var_pair) {
  var1 = var_pair[1]
  var2 = var_pair[2]
  out = ((n1-1)/var1) * 1 / (2^((n1-1)/2) * gamma((n1-1)/2)) * ((n1-1) * s1_j/var1)^((n1-3)/2) * exp((-1/2) * (n1-1) * s1_j/var1) * ((n2-1)/var2) * 1 / (2^((n2-1)/2) * gamma((n2-1)/2)) * ((n2-1) * s2_j/var2)^((n2-3)/2) * exp((-1/2) * (n2-1) * s2_j/var2)
}

p_given_var1_var2 = function(n1, n2, Z1, Z2, var_pair) {
  
  upper = (Z1 - Z2)
  lower = sqrt(var_pair[1]/n1 + var_pair[2]/n2)
  
  test_stat = as.numeric(unlist(upper/lower))
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  
  out = p
}

mass_given_tau_j_2D = function(var_pairs, mass, s1_j, s2_j, n1, n2) {
  f_s_j_given_grid = p_s_j_given_sigma2(n1, n2, s1_j, s2_j,var_pairs)
  f_s_j_grid = f_s_j_given_grid * mass
  post_mass = f_s_j_grid / sum(f_s_j_grid)
  return(post_mass)
}

p_value_npmle_2D_j = function(n1, n2, Z1, Z2, var_pairs, mass, s1_j, s2_j) {
  P_value_joint_list = p_given_var1_var2(n1, n2, Z1, Z2, var_pairs)
  post_mass_j = mass_given_tau_j_2D(var_pairs, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list * post_mass_j)
  return (P_value_j)
}

NPMLE_2D = function(S1_list, S2_list, B1, B2, m, n1, n2, lower_quantile, upper_quantile) {
  lower1 = quantile(S1_list, lower_quantile)
  upper1 = quantile(S1_list, upper_quantile)
  log_u1 = seq(log(lower1), log(upper1), length = B1)
  u1 = exp(log_u1)
  
  lower2 = quantile(S2_list, lower_quantile)
  upper2 = quantile(S2_list, upper_quantile)
  log_u2 = seq(log(lower2), log(upper2), length = B2)
  u2 = exp(log_u2)
  
  var_df = data.frame('var1' = rep(u1, each = B2), 'var2' = rep(u2, B1))
  
  d = rep(1, B1 * B2)
  w = rep(1, m)/m
  
  A = matrix(0, nrow = m, ncol = B1 * B2)
  for (i in 1:(B1*B2)) {
    A[, i] = p_s_j_given_sigma2(n1, n2, S1_list, S2_list, c(var_df[i, 1], var_df[i, 2]))
  }
  
  result = KWPrimal(A, d, w)
  mass =  result$f/sum(result$f)
  
  output = list('grid' = var_df, 'mass' = mass)
  
  return (output)
  
}

P_value_2D_NPMLE = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, grid, mass) {
  
  P_value_list_npmle_2D = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle_2D[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], grid, mass, S1_list[i], S2_list[i])
  }
  
  return (P_value_list_npmle_2D)
}

# P value calculation for 1D projection from 2D NPMLE

P_value_1D_projection = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, grid, mass) {
  tau = S1_list/S2_list
  P_value_list_npmle_projection = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle_projection[i] = p_value_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], grid, mass, tau[i])
  }
  
  return (P_value_list_npmle_projection)
}

# P value calculation for Welch Approximation

Welch_approximation_p = function(n1, n2, Z1, Z2, S1, S2) {
  dfw = (S1/n1 + S2/n2)^2/(1 / (n1-1) * (S1/n1)^2 + 1 / (n2-1) * (S2/n2)^2)
  se2 = S1/n1 + S2/n2
  tw = (Z1 - Z2)/sqrt(se2)
  p = pt(q = abs(tw), df = dfw, lower.tail = FALSE) * 2
  return (p)
}

P_value_Welch = function (n1, n2, m, Z1_list, Z2_list, S1_list, S2_list) {
  
  P_value_list_Welch = rep(0, m)
  for (i in 1:m) {
    P_value_list_Welch[i] = Welch_approximation_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i])
  }
  
  return (P_value_list_Welch)
}

# P value calculation for Beherens-Fisher Test

P_value_BF = function(m, X1, X2) {
  
  P_list = rep(0, m)
  
  for (i in c(1:m)) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_i = bfTest(X1_i, X2_i)$p.value
    P_list[i] = P_i
  }
  
  return (P_list)
}

# P value calculation for Equal-variance t-test

pooled_t_test_p = function(n1, n2, Z1, Z2, S1, S2) {
  Spool = sqrt(((n1-1) * S1 + (n2-1) * S2) / (n1 + n2 - 2))
  se = Spool * sqrt(1/n1 + 1/n2)
  t = (Z1 - Z2)/se
  p = pt(q = abs(t), df = n1 + n2 -2, lower.tail = FALSE) * 2
  return (p)
}

P_value_pooled_t_test = function (n1, n2, m, Z1_list, Z2_list, S1_list, S2_list) {
  
  P_value_list_pooled_t_test = rep(0, m)
  for (i in 1:m) {
    P_value_list_pooled_t_test[i] = pooled_t_test_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i])
  }
  
  return (P_value_list_pooled_t_test)
}

# P value calculation for Equal-variance-distr NPMLE

p_s_j_given_sigma2_EV = function(n, s_j, var) {
  out = ((n-1)/var) * 1 / (2^((n-1)/2) * gamma((n-1)/2)) * ((n-1) * s_j/var)^((n-3)/2) * exp((-1/2) * (n-1) * s_j/var)
}

EVD_NPMLE_1D = function(S1_list, S2_list, B, m, n1, n2, lower_quantile, upper_quantile) {
  
  S_list = c(S1_list, S2_list)
  
  lower = quantile(S_list, lower_quantile)
  upper = quantile(S_list, upper_quantile)
  log_u = seq(log(lower), log(upper), length = B)
  u = exp(log_u)
  d = rep(1,B)
  w = rep(1, m * 2) / (m * 2) 
  A1 = outer(S1_list, u, FUN = p_s_j_given_sigma2_EV, n = n1)
  A2 = outer(S2_list, u, FUN = p_s_j_given_sigma2_EV, n = n2)
  A = rbind(A1, A2)
  result = KWPrimal(A, d, w)
  mass = result$f/sum(result$f)
  
  var_df = data.frame('var1' = rep(u, each = B), 'var2' = rep(u, B))
  pair_mass = as.vector(outer(mass, mass, FUN = "*"))
  
  output = list('grid' = var_df, 'mass' = pair_mass)
}

P_value_EVD_NPMLE = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, EVD_NPMLE_1D_parameter) {
  EVD_NPMLE_result = EVD_NPMLE_1D(S1_list, S2_list, B = EVD_NPMLE_1D_parameter[1], m, n1, n2, lower_quantile = EVD_NPMLE_1D_parameter[2], upper_quantile = EVD_NPMLE_1D_parameter[3])
  
  P_value_list_EVD_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_EVD_npmle[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], EVD_NPMLE_result$grid, EVD_NPMLE_result$mass, S1_list[i], S2_list[i])
  }
  
  return(P_value_list_EVD_npmle)
}

# P value calculation for Equal-variance NPMLE

p_s_j_given_sigma2 = function(n, s_j, var) {
  out = ((n-1)/var) * 1 / (2^((n-1)/2) * gamma((n-1)/2)) * ((n-1) * s_j/var)^((n-3)/2) * exp((-1/2) * (n-1) * s_j/var)
}

EV_NPMLE_1D = function(S1_list, S2_list, B, m, n1, n2, lower_quantile, upper_quantile) {
  
  S_list = c(S1_list, S2_list)
  
  lower = quantile(S_list, lower_quantile)
  upper = quantile(S_list, upper_quantile)
  log_u = seq(log(lower), log(upper), length = B)
  u = exp(log_u)
  d = rep(1,B)
  w = rep(1, m * 2) / (m * 2) 
  A1 = outer(S1_list, u, FUN = p_s_j_given_sigma2_EV, n = n1)
  A2 = outer(S2_list, u, FUN = p_s_j_given_sigma2_EV, n = n2)
  A = rbind(A1, A2)
  result = KWPrimal(A, d, w)
  mass = result$f/sum(result$f)
  
  output = list('grid' = u, 'mass' = mass)
}

p_given_sigma2 = function(n1, n2, Z1, Z2, var) {
  
  upper = (Z1 - Z2)
  lower = sqrt(var * (1/n1 + 1/n2))
  
  test_stat = upper/lower
  
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  
  out = p
}

p_s1_s2_given_sigma2 = function(n1, n2, s1_j, s2_j, var) {
  p_s1_given_sigma2 = p_s_j_given_sigma2_EV(n1, s1_j, var)
  p_s2_given_sigma2 = p_s_j_given_sigma2_EV(n2, s2_j, var)
  return (p_s1_given_sigma2 * p_s2_given_sigma2)
}

mass_given_s1_s2 = function(grid, mass, s1_j, s2_j, n1, n2) {
  p_s1_s2_given_grid = p_s1_s2_given_sigma2(n1, n2, s1_j, s2_j, grid)
  p_s1_s2_grid = p_s1_s2_given_grid * mass
  post_mass = p_s1_s2_grid / sum(p_s1_s2_grid)
  return(post_mass)
}

p_value_EV_npmle_j = function(n1, n2, Z1, Z2, grid, mass, s1_j, s2_j) {
  P_value_joint_list_j = p_given_sigma2(n1, n2, Z1, Z2, grid)
  post_mass_j = mass_given_s1_s2(grid, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list_j * post_mass_j)
  return (P_value_j)
}

P_value_EV_NPMLE = function(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, EV_NPMLE_1D_parameter) {
  EV_NPMLE_result = EV_NPMLE_1D(S1_list, S2_list, B = EV_NPMLE_1D_parameter[1], m, n1, n2, lower_quantile = EV_NPMLE_1D_parameter[2], upper_quantile = EV_NPMLE_1D_parameter[3])
  
  P_value_list_EV_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_EV_npmle[i] = p_value_EV_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], EV_NPMLE_result$grid, EV_NPMLE_result$mass, S1_list[i], S2_list[i])
  }
  
  return(P_value_list_EV_npmle)
}

# BH, Power, FDR

my_BH = function(P, alpha) {
  n = length(P)
  sorted_index = order(P, decreasing = FALSE)
  k_list = which(P[sorted_index][1:n] <= (1:n) * alpha / n)
  if (length(k_list) == 0) {
    return (c())
  }
  else {
    k = max(k_list)
    discovery = sorted_index[1:k]
    return (discovery)
  }
}

Power = function(discovery, flag_list) {
  dis_count = sum(flag_list)
  true_dis_count = sum(flag_list[discovery])
  power = true_dis_count/dis_count
  return (power)
}

FDP = function(discovery, flag_list) {
  #null_count = length(flag_list) - sum(flag_list)
  false_dis = length(discovery) - sum(flag_list[discovery])
  fdp = false_dis/max(length(discovery) ,1)
  return (fdp)
}

# Data Generator

data_generator = function(n1, n2, data_generation_parameter, equ_var) {
  
  n1=n1
  n2=n2
  k=data_generation_parameter$k
  d1=data_generation_parameter$d1
  d2=data_generation_parameter$d2
  m=data_generation_parameter$m
  mu1=data_generation_parameter$mu1
  mu2=data_generation_parameter$mu2
  mean_var2 = data_generation_parameter$mean_var2
  var_var2 = data_generation_parameter$var_var2
  pi0 = data_generation_parameter$pi0
  mu0 = data_generation_parameter$mu0
  
  null_count = as.integer(pi0 * m)
  dis_count = m - null_count
  flag_list = c(rep(0, null_count), rep(1, dis_count)) # 0 for null, 1 for truth
  
  lambda = k * rf(m, d1, d2)
  var2 = abs(rnorm(m, mean_var2, sqrt(var_var2)))
  var1 = lambda * var2
  if (equ_var) {
    var1 = var2
  }
  
  X1 = matrix(0, m, n1)
  X2 = matrix(0, m, n2)
  
  for (j in 1:m) {
    if (j <= null_count) {
      X1[j, ] = rnorm(n1, mu0, sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mu0, sqrt(var2[j]))
    }
    else {
      X1[j, ] = rnorm(n1, rnorm(1, mean = 0, sd = sqrt(mu1 * var1[j])), sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mu2, sqrt(var2[j]))
    }
  }
  
  output = list('X1' = X1, 'X2' = X2, 'flag_list' = flag_list)
  
  return (output)
}

information_extractor = function(X1, X2) {
  
  n1 = ncol(X1)
  n2 = ncol(X2)
  m = nrow(X1)
  
  Z1_list = rowMeans(X1)
  Z2_list = rowMeans(X2)
  S1_list = apply(X1,1,var)
  S2_list = apply(X2,1,var)
  
  information = list('n1' = n1, 'n2' = n2, 'm' = m, 'Z1_list' = Z1_list, 'Z2_list' = Z2_list, 'S1_list' = S1_list, 'S2_list' = S2_list)
  
  return (information)
}

dir_name = function(equ_var, n1, n2, data_generation_parameter, NPMLE_1D_parameter, NPMLE_2D_parameter, alpha) {
  
  equ = 'unequal'
  if (equ_var) {
    equ = 'equal'
  }
  
  n1=n1
  n2=n2
  k=data_generation_parameter$k
  d1=data_generation_parameter$d1
  d2=data_generation_parameter$d2
  m=data_generation_parameter$m
  mu1=data_generation_parameter$mu1
  mu2=data_generation_parameter$mu2
  mean_var2 = data_generation_parameter$mean_var2
  var_var2 = data_generation_parameter$var_var2
  pi0 = data_generation_parameter$pi0
  mu0 = data_generation_parameter$mu0
  
  B = NPMLE_1D_parameter[1]
  l1 = NPMLE_1D_parameter[2]
  u1 = NPMLE_1D_parameter[3]
  
  B1 = NPMLE_2D_parameter[1]
  B2 = NPMLE_2D_parameter[2]
  l2 = NPMLE_2D_parameter[3]
  u2 = NPMLE_2D_parameter[4]
  
  base_dir = paste('Simulation_result_local/', equ, sep = '')
  
  if (!dir.exists(base_dir)) {
    print('Create Base Directory')
    dir.create(base_dir)
  }
  
  head = paste(base_dir, '/(', paste(n1, n2, k, d1, d2,  m,  mu1,  mu2,  mean_var2,  var_var2,  pi0, mu0, B, l1, u1, B1, B2, l2, u2, alpha, rounds, sep = ','), ')', sep = '')
  
  return(head)
}

file_name = function(rounds, algorithm_list) {
  time = Sys.time()
  file = paste(time, ': R=', rounds, '', sep = '')
  for (code in algorithm_list) {
    file = paste(file, code, sep = ',')
  }
  return (paste('(', file, ')', sep = ''))
}

simulator = function(seed, data_generation_parameter, NPMLE_1D_parameter, NPMLE_2D_parameter, EVD_NPMLE_1D_parameter, EV_NPMLE_1D_parameter, alpha, rounds, algorithm_list, equ_var) {
  
  if (!dir.exists('Simulation_result')) {
    print('Create Data Directory')
    dir.create('Simulation_result')
  }
  
  set.seed(seed)
  
  args = commandArgs(TRUE)
  
  n1 = as.integer(args[1])
  n2 = as.integer(args[2])
  
  dir = dir_name(equ_var, n1, n2, data_generation_parameter, NPMLE_1D_parameter, NPMLE_2D_parameter, alpha)
  file = file_name(rounds, algorithm_list)
  
  if (!dir.exists(dir)) {
    print('Create Directory')
    dir.create(dir)
  }
  
  algorithm_name = c('1D_MLE', '1D_NPMLE', '2D_NPMLE', '1D_Proj', 'Welch', 'B_F', 't_test', 'EV_NPMLE', 'EV2_NPMLE')
  
  FDP_of_algorithms = matrix(0, 9, rounds)
  Power_of_algorithms = matrix(0, 9, rounds)
  
  print('Simulation Start')
  print(paste('Parameters:', file.path(dir, file)))
  
  #P_list = NA
  
  r = 1
  
  while(r <= rounds) {
    
    print(paste('Start of round', r))
    
    Rerun = FALSE
    
    output_r = data_generator(n1, n2, data_generation_parameter, equ_var)
    X1 = output_r$X1
    X2 = output_r$X2
    flag_list = output_r$flag_list
    #flag_list = rep(1, m)
    information = information_extractor(X1, X2)
    n1 = information$n1
    n2 = information$n2
    m = information$m
    Z1_list = information$Z1_list
    Z2_list = information$Z2_list
    S1_list = information$S1_list
    S2_list = information$S2_list
    
    NPMLE_2D_result = NA
    
    for (code in algorithm_list) {
      
      if (code == 1) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_1D_PMLE(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list)
        #hist(P_list, main = '1D_PMLE', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      else if (code == 2) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_1D_NPMLE(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, NPMLE_1D_parameter)
        #hist(P_list, main = '1D_NPMLE', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 3) {
        print(paste('start of', algorithm_name[code]))
        if (sum(is.na(NPMLE_2D_result)) > 0) {
          NPMLE_2D_result = NPMLE_2D(S1_list = S1_list, S2_list = S2_list, m = m, n1 = n1, n2 = n2, B1 = NPMLE_2D_parameter[1], B2 = NPMLE_2D_parameter[2], lower_quantile = NPMLE_2D_parameter[3], upper_quantile = NPMLE_2D_parameter[4])
        }
        P_list = P_value_2D_NPMLE(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, NPMLE_2D_result$grid, NPMLE_2D_result$mass)
        #hist(P_list, main = '2D_NPMLE', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        
        if (power == 0) {
          print(paste('Error detected and rerun round', r))
          Rerun = TRUE
          break
        }
        
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 4) {
        print(paste('start of', algorithm_name[code]))
        if (sum(is.na(NPMLE_2D_result)) > 0) {
          NPMLE_2D_result = NPMLE_2D(S1_list = S1_list, S2_list = S2_list, m = m, n1 = n1, n2 = n2, B1 = NPMLE_2D_parameter[1], B2 = NPMLE_2D_parameter[2], lower_quantile = NPMLE_2D_parameter[3], upper_quantile =NPMLE_2D_parameter[4])
        }
        lambda_projection_list = NPMLE_2D_result$grid[, 1] / NPMLE_2D_result$grid[, 2]
        P_list = P_value_1D_projection(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, lambda_projection_list, NPMLE_2D_result$mass)
        #hist(P_list, main = '1D_Projection', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 5) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_Welch(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list)
        #hist(P_list, main = 'Welch', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 6) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_BF(m, X1, X2)
        #hist(P_list, main = 'Welch', breaks = seq(0, 1, 0.025))
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 7) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_pooled_t_test(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list)
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 8) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_EVD_NPMLE(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, EVD_NPMLE_1D_parameter)
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        
        if (power == 0) {
          print(paste('Error detected and rerun round', r))
          Rerun = TRUE
          break
        }
        
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
      
      if (code == 9) {
        print(paste('start of', algorithm_name[code]))
        P_list = P_value_EV_NPMLE(n1, n2, m, Z1_list, Z2_list, S1_list, S2_list, EV_NPMLE_1D_parameter)
        discovery = my_BH(P_list, alpha)
        power = Power(discovery, flag_list)
        fdp = FDP(discovery, flag_list)
        
        if (power == 0) {
          print(paste('Error detected and rerun round', r))
          Rerun = TRUE
          break
        }
        
        Power_of_algorithms[code, r] = power
        FDP_of_algorithms[code, r] = fdp
        print(c('round' = r, 'algorithm' = algorithm_name[code], 'power' = power, 'fdp' = fdp))
      }
    }
    
    if (Rerun) {
      next
    }
    
    print(paste('End of round', r))
    print('')
    
    print('Updating Round Data')
    
    power_df_r = data.frame('Round (Power)' = c(1:r), '1D_MLE' = Power_of_algorithms[1, 1:r], '1D_NPMLE' = Power_of_algorithms[2, 1:r], '2D_NPMLE' = Power_of_algorithms[3, 1:r], '1D_Proj' = Power_of_algorithms[4, 1:r], 'Welch' = Power_of_algorithms[5, 1:r], 'B_F' = Power_of_algorithms[6, 1:r], 't_test' = Power_of_algorithms[7, 1:r], 'EV_NPMLE' = Power_of_algorithms[8, 1:r], 'EV2_NPMLE' = Power_of_algorithms[9, 1:r])
    write.csv(power_df_r, file = paste(file.path(dir, file), '_power.csv', sep = ''))
    fdp_df_r = data.frame('Round (FDP)' = c(1:r), '1D_MLE' = FDP_of_algorithms[1, 1:r], '1D_NPMLE' = FDP_of_algorithms[2, 1:r], '2D_NPMLE' = FDP_of_algorithms[3, 1:r], '1D_Proj' = FDP_of_algorithms[4, 1:r], 'Welch' = FDP_of_algorithms[5, 1:r], 'B_F' = FDP_of_algorithms[6, 1:r], 't_test' = FDP_of_algorithms[7, 1:r], 'EV_NPMLE' = FDP_of_algorithms[8, 1:r], 'EV2_NPMLE' = FDP_of_algorithms[9, 1:r])
    write.csv(fdp_df_r, file = paste(file.path(dir, file), '_fdp.csv', sep = ''))
    
    r = r + 1
    
  }
  
  print('Simulation Over')
  print('Saving Simulation Result')
  
  Power_list = rowMeans(Power_of_algorithms)
  FDR_list = rowMeans(FDP_of_algorithms)
  simulation_result = data.frame('Algorithm' = algorithm_name[algorithm_list], 'Power' = Power_list[algorithm_list], 'FDR' = FDR_list[algorithm_list])
  write.csv(simulation_result, file = paste(file.path(dir, file), '_summary.csv', sep = ''))
  
  return (simulation_result)
}

alpha = 0.1
rounds = 2
NPMLE_1D_parameter = c(1000, 0.01, 1.0)
NPMLE_2D_parameter = c(80, 80, 0.01, 1.0)
EVD_NPMLE_1D_parameter = c(80, 0.01, 1.0)
EV_NPMLE_1D_parameter = c(1000, 0.01, 1.0)
algorithm_list = c(9)
equ_var = FALSE
seed = Sys.time()
data_generation_parameter = data.frame('k' = 2, 'd1' = 8, 'd2' = 12, 'm' = 5000, 'mu1' = 12, 'mu2' = 0, 'mean_var2' = 6, 'var_var2' = 4, 'pi0' = 0.9, 'mu0' = 0)

result = simulator(seed, data_generation_parameter, NPMLE_1D_parameter, NPMLE_2D_parameter, EVD_NPMLE_1D_parameter, EV_NPMLE_1D_parameter, alpha, rounds, algorithm_list, equ_var)

print(result)
