library(REBayes)
library(Rmosek)
library(asht)

# P value calculator for 1D Non-Parametric MLE

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

p_given_lambda = function(n1, n2, Z1, Z2, S1, S2, W1, W2, lambda) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  tau = S1/S2
  
  c = tau / (tau + nA/nB)
  gamma = lambda / (lambda + nA/nB)
  
  phi = sqrt((((1-c) / (1-gamma)) * (n2-1) + (c / gamma) * (n1-1)) / (n1 + n2 - 2))
  
  left = (Z1 - Z2) / sqrt(S1/nA + S2/nB)
  
  test_stat = left/phi
  
  p = pt(q = abs(test_stat), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  
  return(p)
}

mass_given_tau_j = function(n1, n2, grid, mass, tau_j) {
  f_tau_j_given_grid = p_tau_j_given_lambda(n1, n2, tau_j, grid)
  f_tau_j_grid = f_tau_j_given_grid * mass
  post_mass = f_tau_j_grid / sum(f_tau_j_grid)
  return(post_mass)
}

p_value_npmle_j = function(n1, n2, Z1, Z2, S1, S2, W1, W2, grid, mass, tau_j) {
  P_value_joint = p_given_lambda(n1, n2, Z1, Z2, S1, S2, W1, W2, grid)
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

P_value_1D_NPMLE = function(info, NPMLE_par) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  tau = S1_list/S2_list
  
  print('NPMLE 1D starts')
  
  NPMLE_result = NPMLE_1D(tau = tau, m = m, n1 = n1, n2 = n2, B = NPMLE_par[1], lower_quantile = NPMLE_par[2], upper_quantile = NPMLE_par[3])
  
  print('NPMLE 1D ends')
  
  P_value_list_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle[i] = p_value_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ], NPMLE_result$grid, NPMLE_result$mass, tau[i])
  }
  
  return(P_value_list_npmle)
}

# P value calculator for 2D Non-Parametric MLE

p_s_j_given_sigma2 = function(n1, n2, s1_j, s2_j, var_pair) {
  var1 = var_pair[1]
  var2 = var_pair[2]
  out = ((n1-1)/var1) * 1 / (2^((n1-1)/2) * gamma((n1-1)/2)) * ((n1-1) * s1_j/var1)^((n1-3)/2) * exp((-1/2) * (n1-1) * s1_j/var1) * ((n2-1)/var2) * 1 / (2^((n2-1)/2) * gamma((n2-1)/2)) * ((n2-1) * s2_j/var2)^((n2-3)/2) * exp((-1/2) * (n2-1) * s2_j/var2)
}

p_given_var1_var2 = function(n1, n2, Z1, Z2, var_pair, W1, W2) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  upper = (Z1 - Z2)
  lower = sqrt(var_pair[1]/nA + var_pair[2]/nB)
  
  test_stat = as.numeric(unlist(upper/lower))
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  
  return(p)
}

mass_given_tau_j_2D = function(var_pairs, mass, s1_j, s2_j, n1, n2) {
  f_s_j_given_grid = p_s_j_given_sigma2(n1, n2, s1_j, s2_j,var_pairs)
  f_s_j_grid = f_s_j_given_grid * mass
  post_mass = f_s_j_grid / sum(f_s_j_grid)
  return(post_mass)
}

p_value_npmle_2D_j = function(n1, n2, Z1, Z2, var_pairs, mass, s1_j, s2_j, W1, W2) {
  P_value_joint_list = p_given_var1_var2(n1, n2, Z1, Z2, var_pairs, W1, W2)
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

P_value_2D_NPMLE = function(info, grid, mass) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_npmle_2D = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle_2D[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], grid, mass, S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  }
  
  return (P_value_list_npmle_2D)
}

# P value calculation for 1D projection from 2D NPMLE

P_value_1D_projection = function(info, grid, mass) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  tau = S1_list/S2_list
  P_value_list_npmle_projection = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle_projection[i] = p_value_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ], grid, mass, tau[i])
  }
  
  return (P_value_list_npmle_projection)
}

# P value calculation for Welch Approximation

Welch_approximation_p = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  
  dfw = (S1/nA + S2/nB)^2/(1 / (n1-1) * (S1/nA)^2 + 1 / (n2-1) * (S2/nB)^2)
  se2 = S1/nA + S2/nB
  tw = (Z1 - Z2)/sqrt(se2)
  p = pt(q = abs(tw), df = dfw, lower.tail = FALSE) * 2
  return (p)
}

P_value_Welch = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_Welch = rep(0, m)
  for (i in 1:m) {
    P_value_list_Welch[i] = Welch_approximation_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ])
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

# P value calculation for Pooled t-test

pooled_t_test_p = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  
  Spool = sqrt(((n1-1) * S1 + (n2-1) * S2) / (n1 + n2 - 2))
  se = Spool * sqrt(1/nA + 1/nB)
  t = (Z1 - Z2)/se
  p = pt(q = abs(t), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  return (p)
}

P_value_pooled_t_test = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_pooled_t_test = rep(0, m)
  for (i in 1:m) {
    P_value_list_pooled_t_test[i] = pooled_t_test_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ])
  }
  
  return (P_value_list_pooled_t_test)
}

# P value calculation for EVD-NPMLE

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

P_value_EVD_NPMLE = function(info, EVD_NPMLE_1D_parameter) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  EVD_NPMLE_result = EVD_NPMLE_1D(S1_list, S2_list, B = EVD_NPMLE_1D_parameter[1], m, n1, n2, lower_quantile = EVD_NPMLE_1D_parameter[2], upper_quantile = EVD_NPMLE_1D_parameter[3])
  
  P_value_list_EVD_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_EVD_npmle[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], EVD_NPMLE_result$grid, EVD_NPMLE_result$mass, S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  }
  return(P_value_list_EVD_npmle)
}

# P value calculation for EV-NPMLE

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

p_given_sigma2 = function(n1, n2, Z1, Z2, var, W1, W2) {
  
  nA = sum(W1)
  nB = sum(W2)
  upper = (Z1 - Z2)
  lower = sqrt(var * (1/nA + 1/nB))
  
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

p_value_EV_npmle_j = function(n1, n2, Z1, Z2, grid, mass, s1_j, s2_j, W1, W2) {
  P_value_joint_list_j = p_given_sigma2(n1, n2, Z1, Z2, grid, W1, W2)
  post_mass_j = mass_given_s1_s2(grid, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list_j * post_mass_j)
  return (P_value_j)
}

P_value_EV_NPMLE = function(info, EV_NPMLE_1D_parameter) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  EV_NPMLE_result = EV_NPMLE_1D(S1_list, S2_list, B = EV_NPMLE_1D_parameter[1], m, n1, n2, lower_quantile = EV_NPMLE_1D_parameter[2], upper_quantile = EV_NPMLE_1D_parameter[3])
  
  P_value_list_EV_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_EV_npmle[i] = p_value_EV_npmle_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], EV_NPMLE_result$grid, EV_NPMLE_result$mass, S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  }
  
  return(P_value_list_EV_npmle)
}

# Extract sufficient statistics

information_extractor = function(X1, X2, W1, W2) {
  
  n1 = ncol(X1)
  n2 = ncol(X2)
  m = nrow(X1)
  
  Z1_list = rowSums(X1 * W1) / rowSums(W1)
  Z2_list = rowSums(X2 * W2) / rowSums(W2)
  S1_list = 1 / (n1-1) * rowSums(W1 * (X1 - Z1_list)^2)
  S2_list = 1 / (n2-1) * rowSums(W2 * (X2 - Z2_list)^2)
  
  information = list('n1' = n1, 'n2' = n2, 'm' = m, 'Z1_list' = Z1_list, 'Z2_list' = Z2_list, 'S1_list' = S1_list, 'S2_list' = S2_list, 'W1_matrix' = W1, 'W2_matrix' = W2)
  
  return (information)
}

# BH

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
