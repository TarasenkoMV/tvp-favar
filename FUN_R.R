# FUN for DFM

standartise = function(x) {
  if (!is.matrix(x)) stop('x is not a matrix')
  
  for (i in 1:ncol(x)) {
    
    x[, i] = (x[, i] - mean(x[, i])) / sqrt(var(x[, i]))
    
  }
  
  return(x)
}

PC_estimate = function(x, k, do_rot = F) {
  prc = prcomp(x, center = F, scale = F)
  
  new_x = prc$x[,1:k]
  
  if (do_rot) {
    new_x = t(prc$rotation[1:k, 1:k] %*% t(prc$x[, 1:k]))
  }
  
  return(new_x)
}

do_rot = function(x, k) {
  prc = prcomp(x, center = F, scale = F)
  rot = prc$rotation[1:k, 1:k]
  return(rot)
}

calc_prior_mean_and_eps = function(data, factors, fast_slow, ordering) {
  
  loadings_pr = matrix(0, ncol(data), ncol(factors))
  res_mat = matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    y = data[, i, drop = F]
    
    if (fast_slow[i]) {
      x = factors
    } else {
      x = factors[, !as.logical(ordering)]
    }
    
    if (fast_slow[i]) {
      loadings_pr[i,] = t(solve(t(x) %*% x) %*% t(x) %*% y)
      res_mat[,i] = y - t(loadings_pr[i,,drop = F] %*% t(x))
    } else {
      loadings_pr[i, !as.logical(ordering)] = t(solve(t(x) %*% x) %*% t(x) %*% y)
      res_mat[,i] = y - t(loadings_pr[i,!as.logical(ordering),drop = F] %*% t(x))
    }
  }
  
  return(list(loadings_pr, res_mat))
}

do_lag = function(x, L) {
  if (!is.matrix(x)) stop('x is not a matrix')
  if (L == 1) {
    
    x = lag(x, 1)
    x = x[-1, ]
    return(x)
  } else {
    
    lag_list = list()
    
    for (i in 1:L) {
      lag_list[[i]] = lag(x, i)
    }
    
    x = Reduce(cbind, lag_list)[-c(1:L),]
    return(x)
  }
}

get_ar_prior = function(y, x) {
  
  T_K = nrow(x) - ncol(x) - 1
  xTx = solve(t(x) %*% x)
  beta = xTx %*% t(x) %*% y
  e = y - x %*% beta
  
  sigma_hat = (t(e) %*% e) / T_K
  
  omega = diag(sigma_hat) %>% 
    lapply(function(z) z * xTx)
  
  return(list(beta, sigma_hat, omega, e))
}

tsprior = function (priordat, nlag, ndraws = 4000) 
{
  M <- dim(priordat)[2]
  K <- M + nlag * (M^2)
  numa <- 0.5 * M * (M - 1)
  yt <- t(priordat[(nlag + 1):(dim(priordat)[1]), ])
  tau <- dim(yt)[2]
  Zt <- makeregs(priordat, 0, nlag)
  vbar <- matrix(0, K, K)
  xhy <- matrix(0, K, 1)
  for (i in 1:tau) {
    zhat1 <- Zt[((i - 1) * M + 1):(i * M), ]
    vbar <- vbar + t(zhat1) %*% zhat1
    xhy <- xhy + t(zhat1) %*% matrix(yt[, i], M, 1)
  }
  aols <- solve(vbar) %*% xhy
  sse2 <- matrix(0, M, M)
  for (i in 1:tau) {
    e <- matrix(yt[, i] - Zt[((i - 1) * M + 1):(i * M), ] %*% 
                  aols, M, 1)
    sse2 <- sse2 + e %*% t(e)
  }
  hbar <- sse2/tau
  hbarinv <- solve(hbar)
  vbar <- matrix(0, K, K)
  for (i in 1:tau) {
    zhat1 <- Zt[((i - 1) * M + 1):(i * M), ]
    vbar <- vbar + t(zhat1) %*% hbarinv %*% zhat1
  }
  vbar <- solve(vbar)
  achol <- t(chol(hbar))
  ssig <- diag(achol)
  for (i in 1:M) {
    for (j in 1:M) {
      achol[j, i] <- achol[j, i]/ssig[i]
    }
  }
  achol <- solve(achol)
  a0 <- matrix(0, numa, 1)
  ic <- 1
  for (i in 2:M) {
    for (j in 1:(i - 1)) {
      a0[ic, 1] <- achol[i, j]
      ic <- ic + 1
    }
  }
  ssig1 <- log(ssig^2)
  hbar1 <- solve(tau * hbar)
  hdraw <- matrix(0, M, M)
  a02mo <- matrix(0, numa, numa)
  a0mean <- matrix(0, numa, 1)
  for (irep in 1:ndraws) {
    hdraw <- solve(rWishart(1, tau, hbar1)[, , 1])
    achol <- t(chol(hdraw))
    ssig <- diag(achol)
    for (i in 1:M) {
      for (j in 1:M) {
        achol[j, i] <- achol[j, i]/ssig[i]
      }
    }
    achol <- solve(achol)
    a0draw <- matrix(0, numa, 1)
    ic <- 1
    for (i in 2:M) {
      for (j in 1:(i - 1)) {
        a0draw[ic, 1] <- achol[i, j]
        ic <- ic + 1
      }
    }
    a02mo <- a02mo + a0draw %*% t(a0draw)
    a0mean <- a0mean + a0draw
  }
  a02mo <- a02mo/ndraws
  a0mean <- a0mean/ndraws
  a02mo <- a02mo - a0mean %*% t(a0mean)
  return(list(B_OLS = aols, VB_OLS = vbar, A_OLS = a0, sigma_OLS = ssig1, 
              VA_OLS = a02mo))
}


makeregs = function(dat, cut, p) {
  t <- dim(dat)[1]
  M <- dim(dat)[2]
  K <- M + p * (M^2)
  Z <- matrix(NA, t * M, K)
  for (i in (cut + p + 1):t) {
    ztemp <- diag(M)
    for (j in 1:p) {
      xtemp <- t(diag(M) %x% dat[i - j, ])
      ztemp <- cbind(ztemp, xtemp)
    }
    Z[((i - 1) * M + 1):(i * M), ] <- ztemp
  }
  return(Z[((cut + p) * M + 1):(t * M), ])
}

getmix = function() {
  q <- c(0.0073, 0.10556, 2e-05, 0.04395, 0.34001, 0.24566, 
         0.2575)
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, 
         -1.08819)
  u2 <- c(5.79596, 2.61369, 5.1795, 0.16735, 0.64009, 0.34023, 
          1.26261)
  return(list(q = q, m = m, u2 = u2))
}

# Compare things


# Try1 - Настройки для VAR части перенесены из прошлого скрипта
# Try2. k_Q_ увеличена до 0.1. Данные не менялись. Снижено число burn in до 30к

save_res = function(model, name) {
  
  beta_est_all = model$beta_est_all
  fact_est_all = model$fact_est_all
  Ht_est_all = model$Ht_est_all
  rho_obs_all = model$rho_obs_all
  hlaste_all = model$hlaste_all
  
  # name = 'try2.xlsx'
  beta_ar_true = matrix(0, T_obs, 3)
  beta_ar_est = matrix(0, T_obs, 3)
  
  for (i in 1:T_obs) {
    beta_ar_true[i,1] = rho_list[[i]][1,1]
    beta_ar_true[i,2] = rho_list[[i]][2,2]
    beta_ar_true[i,3] = rho_list[[i]][3,3]
    
    beta_ar_est[i,1] = beta_est_all[5, i]
    beta_ar_est[i,2] = beta_est_all[10, i]
    beta_ar_est[i,3] = beta_est_all[15, i]
  }
  
  write.xlsx(cbind(beta_ar_true[,1],
                   beta_ar_est[,1]),
             name,
             sheetName = 'ar1')
  
  write.xlsx(cbind(beta_ar_true[,2],
                   beta_ar_est[,2]),
             name,
             sheetName = 'ar2',
             append = T)
  
  write.xlsx(cbind(beta_ar_true[,3],
                   beta_ar_est[,3]),
             name,
             sheetName = 'ar3',
             append = T)
  
  c21_est = matrix(0, T_obs, 1)
  c31_est = matrix(0, T_obs, 1)
  c32_est = matrix(0, T_obs, 1)
  
  c21_tr = matrix(0, T_obs, 1)
  c31_tr = matrix(0, T_obs, 1)
  c32_tr = matrix(0, T_obs, 1)
  
  for (i in 1:T_obs) {
    c21_tr[i, 1] = R_list[[i]][2, 1]
    c31_tr[i, 1] = R_list[[i]][3, 1]
    c32_tr[i, 1] = R_list[[i]][3, 2]
    
    row1 = (i - 1) * 4 + 1
    row2 = i * 4
    c21_est[i, 1] = Ht_est_all[row1:row2, ][2,1]
    c31_est[i, 1] = Ht_est_all[row1:row2, ][3,1]
    c32_est[i, 1] = Ht_est_all[row1:row2, ][3,2]
  }
  
  write.xlsx(cbind(c21_tr[,1],
                   c21_est[,1]),
             name,
             sheetName = 'c21',
             append = T)
  
  write.xlsx(cbind(c31_tr[,1],
                   c31_est[,1]),
             name,
             sheetName = 'c31',
             append = T)
  
  write.xlsx(cbind(c32_tr[,1],
                   c32_est[,1]),
             name,
             sheetName = 'c32',
             append = T)
  
  rho_obs_tr1 = matrix(0, T_obs, 1)
  rho_obs_tr2 = matrix(0, T_obs, 1)
  rho_obs_tr3 = matrix(0, T_obs, 1)
  rho_obs_tr4 = matrix(0, T_obs, 1)
  rho_obs_tr5 = matrix(0, T_obs, 1)
  
  rho_obs_est1 = matrix(0, T_obs, 1)
  rho_obs_est2 = matrix(0, T_obs, 1)
  rho_obs_est3 = matrix(0, T_obs, 1)
  rho_obs_est4 = matrix(0, T_obs, 1)
  rho_obs_est5 = matrix(0, T_obs, 1)
  
  for (i in 1:T_obs) {
    rho_obs_tr1[,1] = Rho_e[,1]
    rho_obs_tr2[,1] = Rho_e[,2]
    rho_obs_tr3[,1] = Rho_e[,3]
    rho_obs_tr4[,1] = Rho_e[,4]
    rho_obs_tr5[,1] = Rho_e[,5]
    
    rho_obs_est1[,1] = rho_obs_all[,1]
    rho_obs_est2[,1] = rho_obs_all[,2]
    rho_obs_est3[,1] = rho_obs_all[,3]
    rho_obs_est4[,1] = rho_obs_all[,4]
    rho_obs_est5[,1] = rho_obs_all[,5]
  }
  
  write.xlsx(cbind(rho_obs_tr1[,1],
                   rho_obs_est1[,1]),
             name,
             sheetName = 'rho_obs_1',
             append = T)
  
  write.xlsx(cbind(rho_obs_tr2[,1],
                   rho_obs_est2[,1]),
             name,
             sheetName = 'rho_obs_2',
             append = T)
  
  write.xlsx(cbind(rho_obs_tr3[,1],
                   rho_obs_est3[,1]),
             name,
             sheetName = 'rho_obs_3',
             append = T)
  
  write.xlsx(cbind(rho_obs_tr4[,1],
                   rho_obs_est4[,1]),
             name,
             sheetName = 'rho_obs_4',
             append = T)
  
  write.xlsx(cbind(rho_obs_tr5[,1],
                   rho_obs_est5[,1]),
             name,
             sheetName = 'rho_obs_5',
             append = T)
  
  var_eps_tr = matrix(0, T_obs, 5)
  var_eps_est = matrix(0, T_obs, 5)
  
  for (i in 1:T_obs) {
    for (j in 1:5) {
      var_eps_tr[i, j] = obs_stoch_vol[i, j]
      var_eps_est[i, j] = hlaste_all[i, j]
    }
  }
  
  write.xlsx(cbind(var_eps_tr[,1],
                   var_eps_est[,1]),
             name,
             sheetName = 'obs_var_1',
             append = T)
  
  write.xlsx(cbind(var_eps_tr[,2],
                   var_eps_est[,2]),
             name,
             sheetName = 'obs_var_2',
             append = T)
  
  write.xlsx(cbind(var_eps_tr[,3],
                   var_eps_est[,3]),
             name,
             sheetName = 'obs_var_3',
             append = T)
  
  write.xlsx(cbind(var_eps_tr[,4],
                   var_eps_est[,4]),
             name,
             sheetName = 'obs_var_4',
             append = T)
  
  write.xlsx(cbind(var_eps_tr[,5],
                   var_eps_est[,5]),
             name,
             sheetName = 'obs_var_5',
             append = T)
  
  fact_tr = Fact
  fact_est = fact_est_all
  
  write.xlsx(cbind(fact_tr[,1],
                   fact_est[,1]),
             name,
             sheetName = 'fact1',
             append = T)
  
  write.xlsx(cbind(fact_tr[,2],
                   fact_est[,2]),
             name,
             sheetName = 'fact2',
             append = T)
  
  write.xlsx(cbind(fact_tr[,3],
                   fact_est[,3]),
             name,
             sheetName = 'fact3',
             append = T)
  
}

insert_col = function(x, col, ordering) {
  new_x = matrix(0, nrow(x), ncol(x) + 1)
  col_pos = which(ordering == 1)
  
  if (col_pos == 1) {
    new_x[,1] = col
    new_x[,-1] = x
  } else if (col_pos == (ncol(x) + 1)) {
    new_x[,1:(col_pos - 1)] = x[,1:(col_pos - 1)]
    new_x[,col_pos] = col
  } else {
    new_x[,1:(col_pos - 1)] = x[,1:(col_pos - 1)]
    new_x[,col_pos] = col
    new_x[,-c(1:col_pos)] = x[,-c(1:(col_pos-1))]
  }
  
  return(new_x)
}



save_model = function(model, i, calc_irf) {
  dir.create(file.path(paste0('favar_res', i)), showWarnings = F)
  path = paste0('~/favar_res', i, '/')
  
  fwrite(model$beta_est_all, paste0(path, 'beta.csv'), sep = ';')
  fwrite(model$fact_est_all, paste0(path, 'fact.csv'), sep = ';')
  fwrite(model$Ht_est_all, paste0(path, 'Ht.csv'), sep = ';')
  fwrite(model$rho_obs_all, paste0(path, 'rho.csv'), sep = ';')
  fwrite(model$hlaste_all, paste0(path, 'hlast.csv'), sep = ';')
  fwrite(model$fload_all, paste0(path, 'fload.csv'), sep = ';')
  
  
  if (calc_irf) {
    for (j in 1:dim(model$irf_mean)[3]) {
      dir.create(file.path(paste0('favar_res', i, '/', j,'/')), showWarnings = F)
      path1 = paste0(path, j, '/')
      fwrite(model$irf_mean[,,j], paste0(path1, j, 'mean.csv'), sep = ';')
      fwrite(model$irf_left_ci[,,j], paste0(path1, j, 'low.csv'), sep = ';')
      fwrite(model$irf_right_ci[,,j], paste0(path1, j, 'high.csv'), sep = ';')
    }
  }
}




estimate_tvp_favar = function(data, Z, k_B_ = 4, k_A_ = 4, k_sig_ = 1, 
                              k_Q_ = 0.1, k_S_ = 0.1, k_W_ = 0.01,
                              Reps = 35000, Burn = 30000,
                              print_num = 250,
                              Z_incl = T,
                              L = 4, tau = 160, K = 5,
                              ordering = c(0, 0, 0, 1),
                              calc_irf = T,
                              horizon = 36,
                              vars_to_calc_irf = NULL,
                              fast_slow = NULL,
                              rotate_prior_fact = T
                              ) {
  
    # data = data; Z = Z; k_B_ = grid$k_B[i]; k_A_ = grid$k_A_[i];
    # k_sig_ = grid$k_sig_[i];
    # k_Q_ = grid$k_Q_[i]; k_S_ = grid$k_S_[i]; k_W_ = grid$k_W_[i];
    # Reps = 16000; Burn = 8000; print_num = 1000; 
    # Z_incl = T; L = 5; tau = 160; K = 5;
    # ordering = c(0, 0, 0, 0, 0, 1);
    # calc_irf = T; horizon = 36;
    # vars_to_calc_irf = vars_to_calc_irf;
    # fast_slow = as.logical(fast_slow)
  
    
    
    # data = data
    # Z = Z
    # k_B_ = grid$k_B[1]
    # k_A_ = grid$k_A_[1]
    # k_sig_ = grid$k_sig_[1]
    # k_Q_ = grid$k_Q_[1]
    # k_S_ = grid$k_S_[1]
    # k_W_ = grid$k_W_[1]
    # Reps = 10000; Burn = 5000; print_num = 50
    # Z_incl = T; L = 2; tau = 160; K = 5;
    # ordering = c(0, 0, 0, 0, 0, 1)
    # calc_irf = T; horizon = 36;
    # vars_to_calc_irf = vars_to_calc_irf;
    # fast_slow = as.logical(fast_slow)
    # rotate_prior_fact = T

    # k_B_ = 4; k_A_ = 4; k_sig_ = 1; k_Q_ = 0.1; k_S_ = 0.1; k_W_ = 0.01;
    # fast_slow = as.logical(fast_slow)
    # Z_incl = T
    # Reps = 10000; Burn = 5000;
    # print_num = 50;
    # Z_incl = T;
    # L = 2; tau = 160; K = 3;
    # ordering = c(0, 0, 0, 0, 0, 1);
    # calc_irf = F;
    # horizon = 36
    
    data = standartise(data)
    
    s = 0
    if (L == 1) {
      s = 1
    }
    
    Lx = 1
    maxdraws = 100
    check = 1
    T_obs = nrow(data)
    NN = ncol(data)
    # K = 3
    
    fact = PC_estimate(data, K, do_rot = rotate_prior_fact)
    our_rot = do_rot(data, K)
    
    if (Z_incl) {
      fact = insert_col(fact, Z, ordering = ordering)
    }
    
    # get prior mean for factor loading
    
    mean_res_list = calc_prior_mean_and_eps(data, fact, fast_slow, ordering)
    
    fload0 = mean_res_list[[1]]
    # mean_res_list - лист, состоящий из 2-х элементов
    # 1-й это не меняющиеся во времени оценки для матрицы,
    # через которую ненаблюдаемые факторы действуют на 
    # наблюдаемые переменные
    # 2-й элемент это остатки соответствующих уравнений
    
    scale = 3.5e-02
    
    # Выборка
    # Хотим приор для ненаблюдаемых переменных
    xy_list = preparex(fact, L, 1)
    y0w = xy_list[[1]]
    x0w = xy_list[[2]]
    
    ols_list = get_ar_prior(y0w, x0w)
    # Остатки по столбцам:
    e_ar_fact = ols_list[[4]]
    
    # Residual in observation eq
    # Хотим приор для epsilon ar(Lx)
    
    res_obs = mean_res_list[[2]]
    
    b00e = list()
    s00e = list()
    p00e = list()
    Q0e = list()
    Qe_list = list()
    hlaste = matrix(0, T_obs+1, NN)
    
    for (i in 1:NN) {
      
      xy_list = preparex(res_obs[, i] %>% as.matrix, Lx, 0)
      y0e = xy_list[[1]]
      x0e = xy_list[[2]]
      
      ols_list = get_ar_prior(y0e, cbind(x0e, 1))
      b00e[[i]] = ols_list[[1]][1, 1]
      
      # b00e[[i]] = 0.5
      s00e[[i]] = ols_list[[2]]
      p00e[[i]] = ols_list[[3]][[1]][1, 1]
      Q0e[[i]] = p00e[[i]] * scale * nrow(y0e)
      Qe_list[[i]] = p00e[[i]] * scale * nrow(y0e)
      hlaste[,i] = ols_list[[2]]
    }
    
    
    # Остатки из obs eq
    
    SS0 = 10  
    g0 = 0.1^2    
    Tg0 = 1 
    ge = matrix(1, NN, 1) * g0 
    
    beta_eps_ar = matrix(0, nrow(data), NN)
    
    for (i in 1:ncol(data)) {
      
      beta_eps_ar[, i] = repmat_cpp(t(b00e[[i]]), T_obs, 1)
      
    }
    
    pmat00 = as.matrix(rep(fact[1, ], L + s)) %>% t
    
    if (Z_incl) {
      vmat00 = matrix(0, (K + 1)*(L+s), (K + 1)*(L+s))
      diag(vmat00[1:(K + 1), 1:(K + 1)]) = 1
    } else {
      vmat00 = matrix(0, K*(L+s), K*(L+s))
      diag(vmat00[1:K, 1:K]) = 1
    }
    
    eps_AR_error = list()
    
    # Альтернативные начальные значения для части с At
    

    
    Qprior_ = tau
    Wprior_ = dim(fact)[2] + 1
    Sprior_ = (1:(dim(fact)[2])) + 1
    
    t_ = dim(fact)[1]
    M_ = dim(fact)[2]
    numa_ = 0.5 * M_ * (M_ - 1)
    K_ = M_ + L * (M_^2)
    y_ = t(fact[(L + 1):t_, ])
    t_ = dim(y_)[2]
    consQ_ = 1e-04
    consS_ = 1e-04
    consH_ = 0.01
    consW_ = 1e-04
    prior_ = tsprior(fact[1:(tau + L), ], L, 2000)
    
    B_0_prmean_ = prior_$B_OLS
    B_0_prvar_ = k_B_ * prior_$VB_OLS
    A_0_prmean_ = prior_$A_OLS
    A_0_prvar_ = k_A_ * prior_$VA_OLS
    sigma_prmean_ = prior_$sigma_OLS
    sigma_prvar_ = k_sig_ * diag(M_)
    Q_prmean_ = ((k_Q_)^2) * Qprior_ * prior_$VB_OLS
    Q_prvar_ = Qprior_
    W_prmean_ = ((k_W_)^2) * Wprior_ * diag(M_)
    W_prvar_ = Wprior_
    S_prmean_ = vector("list", M_ - 1)
    S_prvar_ = matrix(0, M_ - 1, 1)
    
    ind <- 1
    for (ii in 2:M_) {
      
      S_prmean_[[ii - 1]] <- ((k_S_)^2) * Sprior_[ii - 1] * prior_$VA_OLS[((ii - 
                                                                              1) + (ii - 3) * (ii - 2)/2):ind, ((ii - 1) + (ii - 
                                                                                                                              3) * (ii - 2)/2):ind]
      S_prvar_[ii - 1, 1] <- Sprior_[ii - 1]
      ind <- ind + ii
      
    }
    
    Ht_ = matrix(1, t_, 1) %x% (consH_ * diag(M_))
    Htchol_ = matrix(1, t_, 1) %x% (sqrt(consH_) * diag(M_))
    Qdraw_ = consQ_ * diag(K_)
    Sdraw_ = consS_ * diag(numa_)
    Sblockdraw_ <- vector("list", M_ - 1)
    ijc <- 1
    for (jj in 2:M_) {
      Sblockdraw_[[jj - 1]] <- Sdraw_[((jj - 1) + (jj - 3) * 
                                         (jj - 2)/2):ijc, 
                                      ((jj - 1) + (jj - 3) * (jj - 2)/2):ijc]
      ijc = ijc + jj
    }
    Wdraw_ <- consW_ * diag(M_)
    Btdraw_ <- matrix(0, K_, t_)
    Atdraw_ <- matrix(0, numa_, t_)
    Sigtdraw_ <- matrix(0, M_, t_)
    sigt_ <- matrix(1, t_, 1) %x% (0.01 * diag(M_))
    statedraw_ <- 5 * matrix(1, t_, M_)
    Zs_ <- matrix(1, t_, 1) %x% diag(M_)
    tmp_ <- getmix()
    q_ <- tmp_$q
    m_ <- tmp_$m
    u2_ <- tmp_$u2
    
    # Предцикл
    jgibbs = 1
    igibbs = 1
    
    kept_factors = list()
    kept_aij = list()
    kept_stv = list()
    kept_Ht = list()
    kept_aij = list()
    kept_beta = list()
    kept_obs_rho = list()
    kept_hlaste = list()
    kept_fload = list()
    
    if (Z_incl) {
      data = cbind(data, Z)
    }
    
    if (calc_irf) {
      kept_irf = list()
      sh_pos = which(ordering == 1)
    }
  
  
  
  # Начало цикла
  while(jgibbs <= Reps) {
    
    if (jgibbs == 1) {
      cat('Now: ', Sys.time() %>% as.character %>% substr(1, 19), ' starting \n', sep = '')
    }
    
    if (jgibbs %% print_num == 0) {
      cat('Now: ', Sys.time() %>% as.character %>% substr(1, 19), " ", jgibbs, '-th iteration \n', sep = '')
    }
    
    #1. Sample TVP-Var
    
    Z_ = do_regs(fact, 0, L)
    y_ = t(fact[(L+1):T_obs,])
    
    Btdraw_11 = carterkohn(y_, Z_, Ht_, Qdraw_,
                           K_, M_, t_,
                           B_0_prmean_, B_0_prvar_)
    
    Btdraw_ = Btdraw_11$bdraws
    
    Btdraw_1 = cbind(repmat_cpp(Btdraw_[, 1, drop = F], 1, L), Btdraw_)
    
    sse_2_ = Btdraw_[, 2:t_] - Btdraw_[, 1:(t_ - 1)]
    sse_2_ = sse_2_ %*% t(sse_2_)
    Qdraw_ = solve(rWishart(1, t_ + Q_prvar_, solve(sse_2_ + Q_prmean_))[, , 1])
    
    
    # 2. Sample sthoch covar matrix
    yhat_ = alphahelper(y_, Z_, Btdraw_)
    Zc_ = -t(yhat_)
    
    sigma2temp_ <- t(exp(Sigtdraw_))
    
    Atdraw_ = {}
    ind_ = 1
    
    for (ii in 2:M_) {
      Atblockdraw_ <- carterkohn(yhat_[ii, , drop = FALSE],
                                 Zc_[, 1:(ii - 1), drop = FALSE],
                                 sigma2temp_[, ii, drop = FALSE],
                                 matrix(Sblockdraw_[[ii - 1]], ii - 1, ii - 1),
                                 ii - 1, 1, t_,
                                 A_0_prmean_[((ii - 1) + (ii - 3) * (ii - 2)/2):ind_, ],
                                 A_0_prvar_[((ii - 1) + (ii - 3) * (ii - 2) / 2):ind_,
                                            ((ii - 1) + (ii - 3) * (ii - 2)/2):ind_, drop = FALSE])$bdraws
      
      Atdraw_ <- rbind(Atdraw_, Atblockdraw_)
      ind_ <- ind_ + ii
    }
    
    sse_2_ <- Atdraw_[, 2:t_] - Atdraw_[, 1:(t_ - 1)]
    sse_2_ <- sse_2_ %*% t(sse_2_)
    ijc_ <- 1
    
    for (jj in 2:M_) {
      Sinv_ <- solve(sse_2_[((jj - 1) + (jj - 3) * (jj -
                                                      2) / 2):ijc_, ((jj - 1) + (jj - 3) * (jj - 2)/2):ijc_] +
                       S_prmean_[[jj - 1]])
      Sblockdraw_[[jj - 1]] <- solve(rWishart(1, t_ - 1 +
                                                S_prvar_[jj - 1], Sinv_)[, , 1])
      ijc_ <- ijc_ + jj
    }
    
    capAt_ <- sigmahelper1(Atdraw_, M_)
    aux_ <- sigmahelper2(capAt_, yhat_, q_, m_, u2_, Sigtdraw_,
                         Zs_, Wdraw_, sigma_prmean_, sigma_prvar_)
    
    Sigtdraw_ = aux_$Sigtdraw
    sigt_ = aux_$sigt
    sse_2 = Sigtdraw_[, 2:t_] - Sigtdraw_[, 1:(t_ - 1)]
    sse_2 = sse_2 %*% t(sse_2)
    Wdraw_ = solve(rWishart(1, t_ + W_prvar_, solve(sse_2 +
                                                      W_prmean_))[, , 1])
    aux2_ = sigmahelper3(capAt_, sigt_)
    Ht_ = aux2_$Ht
    
    if (Z_incl) {
      Ht_1 = rbind(repmat_cpp(Ht_[1:(K + 1), 1:(K + 1)], (M_ * T_obs - nrow(Ht_)) / (K+1), 1), Ht_)
    }
    
    #3. Sample TVP-Ar for errors in obs eq
    
    # cat('3')
    
    for (i in 1:NN) {
      res = res_obs[,i,drop = F] - mean(res_obs[,i,drop = F])
      res = rbind(res[1:Lx,], res)
      yx_list = preparex(res, Lx, 0)
      
      ye = yx_list[[1]]
      xe = yx_list[[2]]
      Qe = Qe_list[[i]]
      hlasteps = hlaste[,i] %>% as.matrix
      b00eps = b00e[[i]]
      p00eps = p00e[[i]]
      
      sampled_AR_list = AR_KFS(ye, xe, Qe %>% as.matrix, hlasteps, t(b00eps),
                               p00eps %>% as.matrix, Lx, check = 0, maxdraws, 0)
      
      if (sampled_AR_list[[4]] == 0) {
        beta_eps_ar[, i] = sampled_AR_list[[1]]
        eps_AR_error[[i]] = sampled_AR_list[[2]]
      } else {
        warning("AR coefficients are unstable in obs residuals")
        if (length(eps_AR_error) == 0) stop('error sry ((')
      }
      #4. Draw Qe
      resbetae = diff(beta_eps_ar[[i]])
      scaleQ = t(resbetae) %*% resbetae + Q0e[[i]]
      Qe_list[[i]] = iwpq(T_obs, inv_c(scaleQ))
    }
    
    #4. Sample SVOL for observation eq
    
    for (i in 1:NN) {
      
      hlasteps = hlaste[,i] %>% as.matrix
      hnew = do_svol(hlasteps, ge[i,1], log(s00e[[i]]),
                     SS0, eps_AR_error[[i]] %>% as.matrix)
      hlaste[,i] = hnew
      
      gerrors = diff(log(hlaste[,i])) %>% as.matrix
      ge[i,1] = IG(Tg0, g0, gerrors)
      
    }
    
    #5. Sample Loading
    
    if (Z_incl) {
      fload = matrix(0, NN + 1, K + 1)
    } else {
      fload = matrix(0, NN, K)
    }
    
    res_obs = matrix(0, T_obs, NN)
    
    # fload[1:K, !as.logical(ordering)] = our_rot
    fload[1:K, !as.logical(ordering)] = diag(1, K)
    xx = fact[,!as.logical(ordering)]
    for (i in 1:K) {
      yy = data[,i] %>% as.matrix
      res_obs[,i] = yy - xx %*% as.matrix(fload[i,!as.logical(ordering)])
    }
    
    for (i in (K+1):NN) {
      if (Z_incl) {
        if (fast_slow[i]) {
          xx = fact
        } else {
          xx = fact[,!as.logical(ordering)]
        }
      } else {
        xx = fact
      }
      
      yy = data[,i] %>% as.matrix
      
      # Автокор
      yys = remSC(yy, beta_eps_ar[,i] %>% as.matrix)
      xxs = remSC(xx, beta_eps_ar[,i] %>% as.matrix)
      
      # Гетероскедастичность
      yyss = yys / sqrt(hlaste[-1,i])
      xxss = xxs / repmat_cpp(sqrt(hlaste[-1,i] %>% as.matrix), 1, K + as.numeric(fast_slow[i]))
      
      yyss = yyss[-1,] %>% as.matrix
      xxss = xxss[-1,] %>% as.matrix
      
      coef_list = fast_reg(yyss, xxss)
      
      Fl = coef_list[[1]]
      
      if (Z_incl) {
        if (fast_slow[i]) {
          fload[i,] = t(Fl)
        } else {
          fload[i,!as.logical(ordering)] = t(Fl)
        }
      } else {
        fload[i,] = t(Fl)
      }
      
      res_obs[,i] = yy - xx %*% Fl
    }
    
    fload[nrow(fload), as.logical(ordering)] = 1
    
    #6. Sample Factors
    
    if (Z_incl) {
      dataF = prepare_for_factors_CK(data, cbind(beta_eps_ar, 0), Lx)
    } else {
      dataF = prepare_for_factors_CK(data, beta_eps_ar, Lx)
    }
    
    if (Z_incl) {
      beta2 = Factors_KFS1(dataF = dataF, fload = fload,
                           beta2e = cbind(beta_eps_ar, 0), Lx = Lx,
                           L = L, pmat00 = pmat00,
                           vmat00 = vmat00, Ar_trans_mat = Btdraw_1,
                           K = K + 1, hlaste = cbind(hlaste, 0), Ht = Ht_1,
                           s = s)
    } else {
      beta2 = Factors_KFS1(dataF = dataF, fload = fload,
                           beta2e = beta_eps_ar, Lx = Lx,
                           L = L, pmat00 = pmat00,
                           vmat00 = vmat00, Ar_trans_mat = Btdraw_1,
                           K = K, hlaste = hlaste, Ht = Ht_1,
                           s = s)
    }
    
    # cat('10')
    if (beta2[[1]][1] == 1) {
      stop('Error in Factors sampler. Forecast covariance is not PD')
    } else {
      beta2 = beta2[[1]]
    }
    
    fact = beta2[, 1:(K + as.numeric(Z_incl))]
    
    if (jgibbs > Burn) {
      kept_fload[[jgibbs - Burn]] = fload
      kept_beta[[jgibbs - Burn]] = Btdraw_1
      kept_factors[[jgibbs - Burn]] = fact
      kept_Ht[[jgibbs - Burn]] = Ht_1
      kept_obs_rho[[jgibbs - Burn]] = beta_eps_ar
      kept_hlaste[[jgibbs - Burn]] = hlaste[-1, ]
      
      if (calc_irf) {
        kept_irf[[jgibbs - Burn]] = calc_irf_cpp(Btdraw_1, Ht = Ht_1, K = K + 1, 
                                                 Fload = fload[as.logical(vars_to_calc_irf),], 
                                                 sh_pos = sh_pos, 
                                                 horizon = horizon
                                                 )
      }
    }
    
    jgibbs = jgibbs + 1
    
  }
  
  cat('Now: ', Sys.time() %>% as.character %>% substr(1, 19), 
      ' calculating Irf \n', sep = '')
  
  beta_est_all = matrix(0, dim(Btdraw_1)[1], dim(Btdraw_1)[2])
  fact_est_all = matrix(0, dim(fact)[1], dim(fact)[2])
  Ht_est_all = matrix(0, dim(Ht_1)[1], dim(Ht_1)[2])
  rho_obs_all = matrix(0, dim(beta_eps_ar)[1], dim(beta_eps_ar)[2])
  hlaste_all = matrix(0, dim(hlaste[-1, ])[1], dim(hlaste[-1, ])[2])
  fload_all = matrix(0, dim(fload)[1], dim(fload)[2])
  
  if (calc_irf) {
    irf_mean = calc_irf_mean_and_quantile(irf_list = kept_irf, 
                                          t_ = dim(Btdraw_1)[2], 
                                          horizon = horizon, 
                                          quantile_vec = c(0.05, 0.95),
                                          NN = sum(vars_to_calc_irf)
    )
    rm(kept_irf)
    gc()
  }
  
  for (i in 1:length(kept_beta)) {
    beta_est_all = beta_est_all + kept_beta[[i]]
    fact_est_all = fact_est_all + kept_factors[[i]]
    Ht_est_all = Ht_est_all + kept_Ht[[i]]
    rho_obs_all = rho_obs_all + kept_obs_rho[[i]]
    hlaste_all = hlaste_all + kept_hlaste[[i]]
    fload_all = fload_all + kept_fload[[i]]
    
    if (i == length(kept_beta)) {
      beta_est_all = beta_est_all / length(kept_beta)
      fact_est_all = fact_est_all / length(kept_beta)
      Ht_est_all = Ht_est_all / length(kept_beta)
      rho_obs_all = rho_obs_all / length(kept_beta)
      hlaste_all = hlaste_all / length(kept_beta)
      fload_all = fload_all / length(kept_beta)
    }
  }
  
  cat('\nNow: ', Sys.time() %>% as.character %>% 
        substr(1, 19), ' Done ', sep = '')
  
  if (calc_irf) {
    return(list(beta_est_all = beta_est_all,
                fact_est_all = fact_est_all,
                Ht_est_all = Ht_est_all,
                rho_obs_all = rho_obs_all,
                hlaste_all = hlaste_all,
                fload_all = fload_all,
                irf_mean = irf_mean[[1]],
                irf_left_ci = irf_mean[[2]],
                irf_right_ci = irf_mean[[3]]
    )
    )
  }
  
  return(list(beta_est_all = beta_est_all,
              fact_est_all = fact_est_all,
              Ht_est_all = Ht_est_all,
              rho_obs_all = rho_obs_all,
              hlaste_all = hlaste_all,
              fload_all = fload_all))
}
