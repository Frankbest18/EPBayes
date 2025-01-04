cutoff_pvalue = function(pvalue_list, alpha) {
  adjusted_pvalue = p.adjust(pvalue_list, method = 'BH')
  cutoff_p_value = max(pvalue_list[adjusted_pvalue <= alpha])
  return (cutoff_p_value)
}

# Helper function of BF

ifunc =  function(x, t_BF, tau, n1, n2, W1, W2){
  nA = sum(W1)
  nB = sum(W2)
  
  return (pt(t_BF * sqrt(1 + nB / nA * tau) + sqrt(nB / nA * tau)*x, n2-1) * dt(x, n1-1))
}

BF_pvalue_BF = function(t_BF, tau, n1, n2, W1, W2) {
  
  p.L = integrate(function(x) ifunc(x, t_BF, tau, n1, n2, W1, W2), -Inf, Inf)$value
  
  return (min(2*p.L, 2*(1-p.L)))
}

BF_pvalue_BF_target_function = function(t_BF, tau, n1, n2, W1, W2, target_p_value) {
  return(BF_pvalue_BF(t_BF, tau, n1, n2, W1, W2) - target_p_value)
}

BF_pvalue_BF_solver_positive = function(tau, n1, n2, W1, W2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_BF_target_function(t_BF, tau, n1, n2, W1, W2, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of 1DNPMLE

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

BF_pvalue_given_lambda_1DNPMLE = function(t_BF, n1, n2, tau, W1, W2, lambda) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  c = tau / (tau + nA/nB)
  gamma = lambda / (lambda + nA/nB)
  
  upper = t_BF * sqrt(n1 + n2 -2)
  lower = sqrt((n1-1) * c / gamma + (n2-1) * (1-c) / (1-gamma))
  
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

BF_pvalue_1DNPMLE = function(t_BF, n1, n2, W1, W2, grid, mass, tau) {
  P_value_joint = BF_pvalue_given_lambda_1DNPMLE(t_BF, n1, n2, tau, W1, W2, grid)
  post_mass = mass_given_tau_j(n1, n2, grid, mass, tau)
  P_value = sum(P_value_joint * post_mass)
  return (P_value)
}

BF_pvalue_1DNPMLE_target_function = function(t_BF, n1, n2, W1, W2, grid, mass, tau, target_p_value) {
  return(BF_pvalue_1DNPMLE(t_BF, n1, n2, W1, W2, grid, mass, tau) - target_p_value)
}

BF_pvalue_1DNPMLE_solver_positive = function(n1, n2, W1, W2, grid, mass, tau, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_1DNPMLE_target_function(t_BF, n1, n2, W1, W2, grid, mass, tau, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of Welch

BF_pvalue_Welch = function(t_BF, tau, n1, n2, W1, W2) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  df_estimate = (tau/nA + 1/nB)^2/(1 / (n1-1) * (tau/nA)^2 + 1 / (n2-1) * (1/nB)^2)
  
  p = pt(abs(t_BF), df = df_estimate, lower.tail = FALSE) * 2
  
  return (p)
}

BF_pvalue_Welch_target_function = function(t_BF, tau, n1, n2, W1, W2, target_p_value) {
  return(BF_pvalue_Welch(t_BF, tau, n1, n2, W1, W2) - target_p_value)
}

BF_pvalue_Welch_solver_positive = function(tau, n1, n2, W1, W2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_Welch_target_function(t_BF, tau, n1, n2, W1, W2, target_p_value), interval = c(0, 100))$root
  
}

# Helper function of Pooled t-test 

BF_pvalue_t = function(t_BF, tau, n1, n2, W1, W2) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  test_statistics = t_BF * sqrt((tau / nA + 1 / nB)/((1/nA + 1/nB) * (((n1-1) * tau + (n2-1))/(n1+n2-2))))
  p = pt(abs(test_statistics), df = n1+n2-2, lower.tail = FALSE) * 2
  
  return (p)
}

BF_pvalue_t_target_function = function(t_BF, tau, n1, n2, W1, W2, target_p_value) {
  return(BF_pvalue_t(t_BF, tau, n1, n2, W1, W2) - target_p_value)
}

BF_pvalue_t_solver_positive = function(tau, n1, n2, W1, W2, target_p_value) {
  
  uniroot(function(t_BF) BF_pvalue_t_target_function(t_BF, tau, n1, n2, W1, W2, target_p_value), interval = c(0, 100))$root
  
}










