library(REBayes)
library(Rmosek)
library(asht)

# P value calculator for 1D Non-Parametric MLE

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

P_value_1D_NPMLE = function(info, NPMLE_par) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
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

NPMLE_2D = function(info, B1, B2, lower_quantile, upper_quantile) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
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
  
  P_value_list_npmle_2D = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_npmle_2D[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], grid, mass, S1_list[i], S2_list[i])
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

P_value_Welch = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
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

# P value calculation for Pooled t-test

P_value_pooled_t_test = function (m, X1, X2) {
  
  P_value_list_pooled_t_test = rep(0, m)
  for (i in 1:m) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_value_list_pooled_t_test[i] = t.test(X1_i, X2_i, var.equal = TRUE)$p.value
  }
  
  return (P_value_list_pooled_t_test)
}

# P value calculation for EV-NPMLE

p_s_j_given_sigma2_EV = function(n, s_j, var) {
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
  
  var_df = data.frame('var1' = rep(u, each = B), 'var2' = rep(u, B))
  pair_mass = as.vector(outer(mass, mass, FUN = "*"))
  
  output = list('grid' = var_df, 'mass' = pair_mass)
}

P_value_EV_NPMLE = function(info, EV_NPMLE_1D_parameter) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  EV_NPMLE_result = EV_NPMLE_1D(S1_list, S2_list, B = EV_NPMLE_1D_parameter[1], m, n1, n2, lower_quantile = EV_NPMLE_1D_parameter[2], upper_quantile = EV_NPMLE_1D_parameter[3])
  
  P_value_list_EV_npmle = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_EV_npmle[i] = p_value_npmle_2D_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], EV_NPMLE_result$grid, EV_NPMLE_result$mass, S1_list[i], S2_list[i])
  }
  return(P_value_list_EV_npmle)
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

# extract sufficient statistics

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

# solver

solver = function(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, EV_NPMLE_1D_parameter, alpha, algorithm_list = c(1,2,3,4,5,6)) {
  
  algorithm_name = c('1D_NPMLE', '2D_NPMLE', 'Welch', 'Pooled_t', 'B_F', 'EV_NPMLE')
  
  print('Solver Starts')
  
  information = information_extractor(X1, X2)
  
  NPMLE_2D_result = NA
  
  P_list_1 = NA
  P_list_2 = NA
  P_list_3 = NA
  P_list_4 = NA
  P_list_5 = NA
  P_list_6 = NA
  
  for (code in algorithm_list) {
    
    if (code == 1) {
      print(paste('start of', algorithm_name[code]))
      P_list_1 = P_value_1D_NPMLE(information, NPMLE_1D_parameter)
    }
    
    if (code == 2) {
      print(paste('start of', algorithm_name[code]))
      if (sum(is.na(NPMLE_2D_result)) > 0) {
        NPMLE_2D_result = NPMLE_2D(information, B1 = NPMLE_2D_parameter[1], B2 = NPMLE_2D_parameter[2], lower_quantile = NPMLE_2D_parameter[3], upper_quantile = NPMLE_2D_parameter[4])        
      }
      
      P_list_2 = P_value_2D_NPMLE(information, NPMLE_2D_result$grid, NPMLE_2D_result$mass)
      #discovery = my_BH(P_list_3, alpha)
    }
    
    if (code == 3) {
      print(paste('start of', algorithm_name[code]))
      P_list_3 = P_value_Welch(information)
      #discovery = my_BH(P_list_5, alpha)
    }
    
    if (code == 4) {
      print(paste('start of', algorithm_name[code]))
      P_list_4 = P_value_pooled_t_test(information$m, X1, X2)
    }
    
    if (code == 5) {
      print(paste('start of', algorithm_name[code]))
      P_list_5 = P_value_BF(information$m, X1, X2)
    }
    
    if (code == 6) {
      print(paste('start of', algorithm_name[code]))
      P_list_6 = P_value_EV_NPMLE(information, EV_NPMLE_1D_parameter)
    }
    
  }
  
  print('Solver Ends')
  print('')
  
  if (sum(is.na(NPMLE_2D_result)) > 0) {
    return (list('1D_NPMLE' = P_list_1, '2D_NPMLE' = P_list_2, 'Welch' = P_list_3, 'Pooled_t' = P_list_4, 'B_F' = P_list_5, 'EV_NPMLE' = P_list_6))
  }
  
  return (list('1D_NPMLE' = P_list_1, '2D_NPMLE' = P_list_2, 'Welch' = P_list_3, 'Pooled_t' = P_list_4, 'B_F' = P_list_5, 'EV_NPMLE' = P_list_6, 'grid' = NPMLE_2D_result$grid, 'mass' = NPMLE_2D_result$mass))
}
