# Generate data

check_stability = function(coef_mat) {
  
  coef_mat = t(coef_mat[-nrow(coef_mat), ])
  
  lower_part = matrix(0, K*(L - 1), K*L)
  lower_part[1:K*(L - 1), 1:K * (L - 1)] = diag(K*(L - 1))
  
  coef_mat = rbind(coef_mat, lower_part)
  
  max_eig = max(abs(eigen(coef_mat)$values))
    
    # cat('\n', max_eig)
    
    if (max_eig > 0.85) {
      return(1)
    }
  
  
  return(0)
}

runif_gen = function(n, a, c) {
  res_vec1 = vector('numeric', n)
  b = c(-a[2], -a[1])
  for (i in 1:n) {
    a_s = runif(1, a[1], a[2])
    b_s = runif(1, b[1], b[2])
    
    ch = (1 + (runif(1) > 0.5))
    
    res_vec1[i] = c(a_s, b_s)[ch]
  }
  
  res_vec2 = vector('numeric', n)
  d = c(-c[2], -c[1])
  for (i in 1:n) {
    a_s = runif(1, c[1], c[2])
    b_s = runif(1, d[1], d[2])
    
    ch = (1 + (runif(1) > 0.5))
    
    res_vec2[i] = c(a_s, b_s)[ch]
  }
  
  return(sample(c(res_vec1, res_vec2), size = n))
}

NN = 10
K = 3
L = 2
T_obs = 215

stab = 1

# Generate VAR Dynamics

# 1. Generate innovations

atj_sh_sd = matrix(runif((K * K - K) / 2, 0, 0.055), (K * K - K) / 2, 1)
str_sh_var_sh_sd = runif(K, 0, 0.01)

a_ij_sh_mat = matrix(0, T_obs, nrow(atj_sh_sd))

for (i in 1:nrow(atj_sh_sd)) {
  a_ij_sh_mat[, i] = rnorm(T_obs, sd = atj_sh_sd[i, 1])
}

a_ij_mat = matrix(0, T_obs, nrow(atj_sh_sd))

for (i in 1:nrow(a_ij_mat)) {
  
  if (i == 1) {
    mean_mat = sample(c(.15, .3, .4, 0, .25), nrow(atj_sh_sd))
    a_ij_mat[1, ] = a_ij_sh_mat[1, ] + mean_mat
  } else {
    
    for (j in 1:ncol(a_ij_mat)) {
      a_ij_mat[i, j] = a_ij_mat[i - 1, j] * 0.92 + 0.08 * mean_mat[j] + 
        a_ij_sh_mat[i, j]
      
    }
  }
}

# Generate stoch vol of str sh

str_stoch_vol_ln = matrix(0, T_obs, K)

for (i in 1:K) {
  
  str_stoch_vol_ln[, i] = rnorm(T_obs, sd = str_sh_var_sh_sd[i])
  
}

start = matrix(c(-0.7, -1.4, -1.8), 1, 3)
# start = matrix(sample(c(-0.7, -2.4, -2.8), replace = F, size = 3), 1, 3)
str_stoch_vol_ln = apply(str_stoch_vol_ln, 2, cumsum) + 
  repmat_cpp(start, 215, 1)

str_stoch_vol = exp(str_stoch_vol_ln)
R_list = list()
A_list = list()
innovations = matrix(0, T_obs, K)
for (i in 1:T_obs) {
  A = diag(K)
  col = 1
  row = 2
  for (z in 1:ncol(a_ij_mat)) {
    A[row, col] = a_ij_mat[i, z]
    
    col = col + 1
    if (col == row) {
      row = row + 1
      col = 1
    }
  }
  H = diag(str_stoch_vol[i, ])
  A_list[[i]] = A
  R = solve(A) %*% H %*% t(solve(A))
  # R = matrix(c(1, .4, .3,
  #              .4, 1, .45,
  #              .3, .45, 1), 3, 3, byrow = T)
  R_list[[i]] = R
  innovations[i, ] = matrix(rnorm(K), 1, K) %*% chol(R)
}

# 2. Generate coefficients

rho_sh_sd = matrix(runif(K * K * L + K, 0, 0.055), K * K * L + K, 1)

rho_list = list()
for (i in 1:T_obs) {
  rho_mat = matrix(0, K * L + 1, K)
  
  if (i == 1) {
    
    stab = 1
    
    while(stab) {
      rho_start1 = matrix(sample(size = K * K,
                                 seq(-0.3, 0.5, by = 0.1), 
                                 replace = T),  
                          K, K)
      
      diag(rho_start1) = sample(size = K,
                                seq(0.6, 0.8, by = 0.1), 
                                replace = T)
      
      rho_start2 = matrix(sample(size = (K * (L - 1) + 1) * K,
                                 seq(-0.2, 0.2, by = 0.05), 
                                 replace = T),  
                          K * (L - 1) + 1, K)
      
      rho_start = rbind(rho_start1, rho_start2)
      rho_mat = matrix(rnorm(K * K * L + K, sd = rho_sh_sd), K * L + 1, K) + 
        rho_start
      stab = check_stability(rho_mat)
    }
    
    rho_list[[i]] = rho_mat
    
  } else {
    
    stab = 1
    while(stab) {
      rho_mat = matrix(rnorm(K * K * L + K, sd = rho_sh_sd), K * L + 1, K) + 
        0.9 * rho_list[[i - 1]] + 0.1 * rho_start
      stab = check_stability(rho_mat)
    }
    
    rho_list[[i]] = rho_mat
  }
  
  cat('\n', i)
}

# Write down start values for B
mu = rho_start[nrow(rho_start), ]
rho1_mat = as.vector(rho_start[1:K, 1:K])
rho2_mat = as.vector(rho_start[(K + 1):(2 * K), 1:K])

B_start = matrix(c(mu, rho1_mat, rho2_mat))

# Generate fact

fact = matrix(0, T_obs, K)

for (i in 1:L) {
  fact[i,] = innovations[i, ]
}

for (i in (L + 1):T_obs) {
  X = as.matrix(c(fact[i - 1,], fact[i - 2,], 1))
  fact[i,] = t(rho_list[[i]]) %*% X + innovations[i, ]
}


# Generate AR and observables
obs_er_ar_coef_sh = runif(NN, 0, 0.04)
shocks_to_var_of_true_error_in_obs_eq = runif(NN, 0, 0.01)

# Now we want to generate errors in observation equation
# I. Generate var(epsilon), where e_t = rho_t x e_{t-1} + epsilon_t
start = sample(c(2.2, 2.3, 2.4, 2.5), NN, replace = T)

obs_stoch_vol_ln = matrix(0, T_obs, NN)

for (i in 1:NN) {
  obs_stoch_vol_ln[,i] = rnorm(T_obs, sd = shocks_to_var_of_true_error_in_obs_eq[i])
}

obs_stoch_vol_ln = apply(obs_stoch_vol_ln, 2, cumsum) - 
  repmat_cpp(t(start), T_obs, 1)

obs_stoch_vol = exp(obs_stoch_vol_ln)

obs_shocks = matrix(0, T_obs, NN)

for (i in 1:nrow(obs_shocks)) {
  for (j in 1:ncol(obs_shocks)) {
    obs_shocks[i, j] = rnorm(1, 0, sqrt(obs_stoch_vol[i, j]))
  }
}

# II. Rho_e

Rho_e = matrix(0, T_obs, NN)
Rho_e_fixed = matrix(0.5, 1, NN)
Rho_e[1, ] = Rho_e_fixed
for (i in 2:nrow(Rho_e)) {
  for (j in 1:ncol(Rho_e)) {
    Rho_e[i, j] = Rho_e[i - 1, j] * 0.95 + 0.05 * 0.5 + 
      rnorm(1, 0, obs_er_ar_coef_sh[j])
  }
}

# III. Generate observable errors

obs_error = matrix(0, T_obs, NN)
obs_error[1, ] = obs_shocks[1, ]
for (i in 2:T_obs) {
  obs_error[i, ] = obs_error[i-1, ] * Rho_e[i, ] + obs_shocks[i, ]
}

# Generate F_load
Fload = matrix(runif_gen(NN*K, c(1,2), c(4,5)), NN, K)
Fload[1:K, 1:K] = diag(K)

Y = Fload %*% t(fact) + t(obs_error)
Y = t(Y)




fwrite(Y, 'syntetic_Y.csv')
# fwrite(factors, 'syntetic_fact.csv')
# fwrite(a_ij_mat, 'syntetic_a_ij.csv')
# fwrite(innovations, 'syntetic_innovations.csv')
# fwrite(str_stoch_vol, 'syntetic_stoch_vol.csv')
# fwrite(R_mat, 'R_mat.csv')