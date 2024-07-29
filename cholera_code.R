# Simplified code to make plots
# Cormac LaPrete, Lindsay Keegan
# 07/29/2024

# load libraries
library(ggplot2)
library(dplyr)

# Functions ####################################################################
# beta as a funciton of R0
beta_func <- function(R0) {
  q <- 0
  
  R0 / ((p_A*nu_A / gamma_A) +
          (p_MU*nu_M/alpha_M) + (p_MT*nu_M / ((1-q)*alpha_M + q*delta*tau)) +
          (p_MU + (p_MT*(1-q)*alpha_M / ((1-q)*alpha_M + q*delta*tau))) * (nu_sh*nu_M / (gamma_M + mu_M)) +
          p_MT * (q*delta*tau / ((1-q)*alpha_M + q*delta*tau)) * (nu_abx*nu_M / gamma_Mabx) +
          (p_SU / alpha_S) + (p_ST / tau) +
          p_SU * (nu_sh / (gamma_S + mu_S)) + p_ST * (nu_abx / gamma_Sabx))
}

# beta as a funciton of all parameters
beta_func_all <- function(R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                          mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                          p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N) {
  q <- 0
  
  R0 / ((p_A*nu_A / gamma_A) +
          (p_MU*nu_M/alpha_M) + (p_MT*nu_M / ((1-q)*alpha_M + q*delta*tau)) +
          (p_MU + (p_MT*(1-q)*alpha_M / ((1-q)*alpha_M + q*delta*tau))) * (nu_sh*nu_M / (gamma_M + mu_M)) +
          p_MT * (q*delta*tau / ((1-q)*alpha_M + q*delta*tau)) * (nu_abx*nu_M / gamma_Mabx) +
          (p_SU / alpha_S) + (p_ST / tau) +
          p_SU * (nu_sh / (gamma_S + mu_S)) + p_ST * (nu_abx / gamma_Sabx))
}

# R (effective) as a function of q
R_eff_func <- function(q) {
  beta * (
    (p_A*nu_A / gamma_A) +
      (p_MU*nu_M/alpha_M) + (p_MT*nu_M / ((1-q)*alpha_M + q*delta*tau)) +
      (p_MU + (p_MT*(1-q)*alpha_M / ((1-q)*alpha_M + q*delta*tau))) * (nu_sh*nu_M / (gamma_M + mu_M)) +
      p_MT * (q*delta*tau / ((1-q)*alpha_M + q*delta*tau)) * (nu_abx*nu_M / gamma_Mabx) +
      (p_SU / alpha_S) + (p_ST / tau) +
      p_SU * (nu_sh / (gamma_S + mu_S)) + p_ST * (nu_abx / gamma_Sabx)
  )
}

# R (effective) as a function of q and beta
R_eff_func_qbeta <- function(q, beta) {
  beta * (
    (p_A*nu_A / gamma_A) +
      (p_MU*nu_M/alpha_M) + (p_MT*nu_M / ((1-q)*alpha_M + q*delta*tau)) +
      (p_MU + (p_MT*(1-q)*alpha_M / ((1-q)*alpha_M + q*delta*tau))) * (nu_sh*nu_M / (gamma_M + mu_M)) +
      p_MT * (q*delta*tau / ((1-q)*alpha_M + q*delta*tau)) * (nu_abx*nu_M / gamma_Mabx) +
      (p_SU / alpha_S) + (p_ST / tau) +
      p_SU * (nu_sh / (gamma_S + mu_S)) + p_ST * (nu_abx / gamma_Sabx)
  )
}

# R (effective) as a function of all parameters
R_eff_func_all <- function(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                           mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                           p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta) {
  beta * (
    (p_A*nu_A / gamma_A) +
      (p_MU*nu_M/alpha_M) + (p_MT*nu_M / ((1-q)*alpha_M + q*delta*tau)) +
      (p_MU + (p_MT*(1-q)*alpha_M / ((1-q)*alpha_M + q*delta*tau))) * (nu_sh*nu_M / (gamma_M + mu_M)) +
      p_MT * (q*delta*tau / ((1-q)*alpha_M + q*delta*tau)) * (nu_abx*nu_M / gamma_Mabx) +
      (p_SU / alpha_S) + (p_ST / tau) +
      p_SU * (nu_sh / (gamma_S + mu_S)) + p_ST * (nu_abx / gamma_Sabx)
  )
}

#### Inf Functions (all variables) #############################################
# R_inf as a function of all variables (via R_eff)
R_inf_func <- function(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                       mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                       p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta) {
  R_eff <- R_eff_func_all(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                          mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                          p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      ((alpha_M*gamma_M / (gamma_M + mu_M)) * (sigma*p_MU / alpha_M +
                                                 sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_A) +
         (sigma*p_SU*gamma_S / (gamma_S + mu_S))
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

# R_abx_inf as a function of all variables (via R_eff)
R_abx_inf_func <- function(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                           mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                           p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta) {
  R_eff <- R_eff_func_all(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                          mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                          p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      (q*delta*tau*(sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_ST)
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

# D_inf as a function of all variables (via R_eff)
D_inf_func <- function(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                       mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                       p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta) {
  R_eff <- R_eff_func_all(q, R0, sigma, gamma_A, gamma_M, gamma_S, gamma_Mabx, gamma_Sabx,
                          mu_M, mu_S, alpha_M, alpha_S, eps_A, p_A, eps_S, eps_MT, eps_ST,
                          p_MU, p_MT, p_SU, p_ST, delta, tau, nu_sh, nu_A, nu_M, nu_abx, N, beta)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      ((alpha_M*mu_M / (gamma_M + mu_M)) * (sigma*p_MU / alpha_M +
                                              sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_SU*mu_S / (gamma_S + mu_S))
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

#### Inf Functions (fucntion of q) #############################################
# R_inf as a function of q (via R_eff)
R_inf_func_q <- function(q) {
  R_eff <- R_eff_func(q)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      ((alpha_M*gamma_M / (gamma_M + mu_M)) * (sigma*p_MU / alpha_M +
                                                 sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_A) +
         (sigma*p_SU*gamma_S / (gamma_S + mu_S))
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

# R_abx_inf as a function of q (via R_eff)
R_abx_inf_func_q <- function(q) {
  R_eff <- R_eff_func(q)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      (q*delta*tau*(sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_ST)
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

# D_inf as a function of q (via R_eff)
D_inf_func_q <- function(q) {
  R_eff <- R_eff_func(q)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      ((alpha_M*mu_M / (gamma_M + mu_M)) * (sigma*p_MU / alpha_M +
                                              sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_SU*mu_S / (gamma_S + mu_S))
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}

# R_abx_inf as a function of q (via R_eff) and beta
R_abx_inf_func_qbeta <- function(q, beta) {
  R_eff <- R_eff_func_qbeta(q, beta)
  
  if(R_eff < 1) {
    return(0)
  } else {
    S_inf_func <- function(S_inf) { S_inf - exp(R_eff*(S_inf - 1)) }
    S_inf <- uniroot(S_inf_func, c(0.00001,0.9999))$root
    
    return(
      (q*delta*tau*(sigma*p_MT / ((1-q)*alpha_M + q*delta*tau)) +
         (sigma*p_ST)
      ) * (1/R_eff) * (log(N) - log(S_inf))
    )
  }
}


# Build data set ###############################################################
# load Sharia's parameters
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# load("~/Downloads/used.param.Rdata")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("used.param.Rdata")
params <- used.param
params_df <- params$`prop.m.abx=0.5`

# Set up data frame with all parameters for each run
df_out <- data.frame(run_num=NA, q=NA, R0=NA, sigma=NA,
                     gamma_A=NA, gamma_M=NA, gamma_S=NA, gamma_Mabx=NA, gamma_Sabx=NA,
                     mu_M=NA, mu_S=NA, alpha_M=NA, alpha_S=NA,
                     eps_A=NA, p_A=NA, eps_S=NA, eps_MT=NA, eps_ST=NA,
                     p_MU=NA, p_MT=NA, p_SU=NA, p_ST=NA,
                     delta=NA, tau=NA,
                     nu_sh=NA, nu_A=NA, nu_M=NA, nu_abx=NA,
                     N=NA, beta=NA,
                     R_inf=NA, R_abx_inf=NA, D_inf=NA, S_inf=NA,
                     stringsAsFactors=FALSE)
df_out <- df_out[1:10100,]

# fill out data frame with Sharia's parameters
js <- 1:100
for(j in js){
  p_list <- params_df[j][[1]]
  
  j_seq <- ((j-1)*100+j):(((j-1)*100+j)+100)
  df_out$run_num[j_seq] <- j
  
  df_out$q[j_seq] <- seq(0, 1, 0.01)
  
  df_out$R0[j_seq] <- p_list["R0"][[1]]
  df_out$sigma[j_seq] <- p_list["sigma"][[1]]
  df_out$gamma_A[j_seq] <- p_list["gamma.a"][[1]]
  df_out$gamma_M[j_seq] <- (1-p_list["theta.m"][[1]]) * p_list["gamma.m"][[1]]
  df_out$gamma_S[j_seq] <- (1-p_list["theta.s"][[1]]) * p_list["gamma.s"][[1]]
  df_out$gamma_Mabx[j_seq] <- p_list["gamma.m.abx"][[1]]
  df_out$gamma_Sabx[j_seq] <- p_list["gamma.s.abx"][[1]]
  df_out$mu_M[j_seq] <- p_list["theta.m"][[1]] * p_list["mu.m"][[1]]
  df_out$mu_S[j_seq] <- p_list["theta.s"][[1]] * p_list["mu.s"][[1]]
  df_out$alpha_M[j_seq] <- p_list["alpha.m"][[1]]
  df_out$alpha_S[j_seq] <- p_list["alpha.s"][[1]]
  df_out$eps_A[j_seq] <- p_list["epsilon.a"][[1]]
  df_out$p_A[j_seq] <- p_list["epsilon.a"][[1]]
  df_out$eps_S[j_seq] <- p_list["epsilon.s"][[1]]
  df_out$eps_MT[j_seq] <- p_list["epsilon.m.T"][[1]]
  df_out$eps_ST[j_seq] <- p_list["epsilon.s.T"][[1]]
  df_out$p_MU[j_seq] <- (1-p_list["epsilon.a"][[1]])*(1-p_list["epsilon.s"][[1]])*(1-p_list["epsilon.m.T"][[1]])
  df_out$p_MT[j_seq] <- (1-p_list["epsilon.a"][[1]])*(1-p_list["epsilon.s"][[1]])*p_list["epsilon.m.T"][[1]]
  df_out$p_SU[j_seq] <- (1-p_list["epsilon.a"][[1]])*p_list["epsilon.s"][[1]]*(1- p_list["epsilon.s.T"][[1]])
  df_out$p_ST[j_seq] <- (1-p_list["epsilon.a"][[1]])*p_list["epsilon.s"][[1]]* p_list["epsilon.s.T"][[1]]
  df_out$delta[j_seq] <- 0.5
  df_out$tau[j_seq] <- p_list["tau"][[1]]
  df_out$nu_sh[j_seq] <- p_list["v.sh"][[1]]
  df_out$nu_A[j_seq] <- p_list["v.a"][[1]]
  df_out$nu_M[j_seq] <- p_list["v.m"][[1]]
  df_out$nu_abx[j_seq] <- p_list["v.abx"][[1]]
  df_out$N[j_seq] <- 1
  df_out$beta[j_seq] <- beta_func_all(p_list["R0"][[1]], p_list["sigma"][[1]], p_list["gamma.a"][[1]],
                                      (1-p_list["theta.m"][[1]]) * p_list["gamma.m"][[1]],
                                      (1-p_list["theta.s"][[1]]) * p_list["gamma.s"][[1]],
                                      p_list["gamma.m.abx"][[1]], p_list["gamma.s.abx"][[1]],
                                      p_list["theta.m"][[1]] * p_list["mu.m"][[1]], p_list["theta.s"][[1]] * p_list["mu.s"][[1]],
                                      p_list["alpha.m"][[1]], p_list["alpha.s"][[1]],
                                      p_list["epsilon.a"][[1]], p_list["epsilon.a"][[1]],
                                      p_list["epsilon.s"][[1]], p_list["epsilon.m.T"][[1]], p_list["epsilon.s.T"][[1]],
                                      (1-p_list["epsilon.a"][[1]])*(1-p_list["epsilon.s"][[1]])*(1-p_list["epsilon.m.T"][[1]]),
                                      (1-p_list["epsilon.a"][[1]])*(1-p_list["epsilon.s"][[1]])*p_list["epsilon.m.T"][[1]],
                                      (1-p_list["epsilon.a"][[1]])*p_list["epsilon.s"][[1]]*(1- p_list["epsilon.s.T"][[1]]),
                                      (1-p_list["epsilon.a"][[1]])*p_list["epsilon.s"][[1]]* p_list["epsilon.s.T"][[1]],
                                      0.5, p_list["tau"][[1]], p_list["v.sh"][[1]], p_list["v.a"][[1]],
                                      p_list["v.m"][[1]], p_list["v.abx"][[1]], 1)
}

# Calculate final sizes
df_out$R_inf <- sapply(1:10100, function(i) R_inf_func(df_out$q[i], df_out$R0[i], df_out$sigma[i],
                                                       df_out$gamma_A[i], df_out$gamma_M[i], df_out$gamma_S[i], df_out$gamma_Mabx[i], df_out$gamma_Sabx[i],
                                                       df_out$mu_M[i], df_out$mu_S[i],
                                                       df_out$alpha_M[i], df_out$alpha_S[i],
                                                       df_out$eps_A[i], df_out$p_A[i], df_out$eps_S[i], df_out$eps_MT[i], df_out$eps_ST[i],
                                                       df_out$p_MU[i], df_out$p_MT[i], df_out$p_SU[i], df_out$p_ST[i],
                                                       df_out$delta[i], df_out$tau[i],
                                                       df_out$nu_sh[i], df_out$nu_A[i], df_out$nu_M[i], df_out$nu_abx[i],
                                                       df_out$N[i], df_out$beta[i]))
df_out$R_abx_inf <- sapply(1:10100, function(i) R_abx_inf_func(df_out$q[i], df_out$R0[i], df_out$sigma[i],
                                                               df_out$gamma_A[i], df_out$gamma_M[i], df_out$gamma_S[i], df_out$gamma_Mabx[i], df_out$gamma_Sabx[i],
                                                               df_out$mu_M[i], df_out$mu_S[i],
                                                               df_out$alpha_M[i], df_out$alpha_S[i],
                                                               df_out$eps_A[i], df_out$p_A[i], df_out$eps_S[i], df_out$eps_MT[i], df_out$eps_ST[i],
                                                               df_out$p_MU[i], df_out$p_MT[i], df_out$p_SU[i], df_out$p_ST[i],
                                                               df_out$delta[i], df_out$tau[i],
                                                               df_out$nu_sh[i], df_out$nu_A[i], df_out$nu_M[i], df_out$nu_abx[i],
                                                               df_out$N[i], df_out$beta[i]))
df_out$D_inf <- sapply(1:10100, function(i) D_inf_func(df_out$q[i], df_out$R0[i], df_out$sigma[i],
                                                       df_out$gamma_A[i], df_out$gamma_M[i], df_out$gamma_S[i], df_out$gamma_Mabx[i], df_out$gamma_Sabx[i],
                                                       df_out$mu_M[i], df_out$mu_S[i],
                                                       df_out$alpha_M[i], df_out$alpha_S[i],
                                                       df_out$eps_A[i], df_out$p_A[i], df_out$eps_S[i], df_out$eps_MT[i], df_out$eps_ST[i],
                                                       df_out$p_MU[i], df_out$p_MT[i], df_out$p_SU[i], df_out$p_ST[i],
                                                       df_out$delta[i], df_out$tau[i],
                                                       df_out$nu_sh[i], df_out$nu_A[i], df_out$nu_M[i], df_out$nu_abx[i],
                                                       df_out$N[i], df_out$beta[i]))

# Add outputs to data frame
R_inf_vals <- aggregate(R_inf ~ q, data=df_out, FUN=mean)
R_abx_inf_vals <- aggregate(R_abx_inf ~ q, data=df_out, FUN=mean)
D_inf_vals <- aggregate(D_inf ~ q, data=df_out, FUN=mean)
df <- data.frame(qs=seq(0,1,0.01), R_inf_vals$R_inf, R_abx_inf_vals$R_abx_inf, D_inf_vals$D_inf)

# Setting parameters to mean values
R0 <- mean(df_out$R0)
sigma <- mean(df_out$sigma)
gamma_A <- mean(df_out$gamma_A)
gamma_M <- mean(df_out$gamma_M)
gamma_S <- mean(df_out$gamma_S)
gamma_Mabx <- mean(df_out$gamma_Mabx)
gamma_Sabx <- mean(df_out$gamma_Sabx)
mu_M <- mean(df_out$mu_M)
mu_S <- mean(df_out$mu_S)
alpha_M <- mean(df_out$alpha_M)
alpha_S <- mean(df_out$alpha_S)
eps_A <- mean(df_out$eps_A)
p_A <- mean(df_out$p_A)
eps_S <- mean(df_out$eps_S)
eps_MT <- mean(df_out$eps_MT)
eps_ST <- mean(df_out$eps_ST)
p_MU <- (1-p_A)*(1-eps_S)*(1-eps_MT)
p_MT <- (1-p_A)*(1-eps_S)*eps_MT
p_SU <- (1-p_A)*eps_S*(1-eps_ST)
p_ST <- (1-p_A)*eps_S*eps_ST
delta <- 0.5
tau <- mean(df_out$tau)
nu_sh <- mean(df_out$nu_sh)
nu_A <- mean(df_out$nu_A)
nu_M <- mean(df_out$nu_M)
nu_abx <- mean(df_out$nu_abx)
N <- 1
beta <- beta_func(R0)


# Plotting #####################################################################
# Plot over multiple R0s
# old one when threshold was 1.378:
# R0s <- c(1.05, 1.1, 1.2, 1.3, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.5, 1.75, 2.0)
# new one when threshold is 1.315:
R0s <- c(1.05, 1.1, 1.2, 1.29, 1.3, 1.31, 1.32, 1.35, 1.38, 1.39, 1.4, 1.41, 1.5, 1.75, 2.0)
R0s <- rev(R0s)
betas <- sapply(R0s, beta_func)

qs <- seq(0, 1, 0.01)
df <- NULL
for(i in 1:length(betas)){
  beta <- betas[i]
  dat <- data.frame(q=qs,
                    prop=qs*delta*tau/((1-qs)*alpha_M + qs*delta*tau),
                    R_inf_vals=sapply(qs, R_inf_func_q), 
                    R_abx_inf_vals=sapply(qs, R_abx_inf_func_q),
                    D_inf_vals=sapply(qs, D_inf_func_q),
                    R0=rep(R0s[i],length(qs)))
  df <- rbind(df, dat)
}

# saving list of blues for color
# https://eyetracking.upol.cz/color/
blues <- c("#0B1829FF","#232C43FF","#2F415FFF","#435E7FFF","#4C6C94FF","#A4ABB0FF","#D6D8D0FF","#FBE3C2FF","#ECB27DFF","#E69C6BFF","#D37750FF","#B9563FFF","#92351EFF")

blues <- c("#005C55FF","#1F6E5EFF","#3F7F68FF","#5E9171FF","#7DA37BFF","#9CB484FF","#BCC68DFF","#FBE3C2FF","#ECB27DFF","#E69C6BFF","#D37750FF","#B9563FFF","#92351EFF")

blues <- c("#081D58FF","#253494FF","#225EA8FF","#1D91C0FF","#41B6C4FF","#7FCDBBFF","#C7E9B4FF", "#FBE3C2FF","#F2C88FFF","#ECB27DFF","#E69C6BFF","#D37750FF","#B9563FFF","#92351EFF", "#662506FF")


R0_list <- R0s[R0s != 1.05]

thresh_dat <- NULL

for(i in 1:length(R0_list)){
  R0_thresh <- df %>% 
    filter(R0 == R0_list[i]) %>% 
    filter(R_abx_inf_vals < R_abx_inf_vals[1]) %>% 
    slice(1)
  thresh_dat <- rbind(thresh_dat, R0_thresh)
}


# all infected (summed) vs. proportion treated
library(cowplot)
#pdf("infected-prop.pdf", width=6, height=6)
plot.A <- ggplot(df) +
  geom_line(aes(x = prop, y = 100 * (R_abx_inf_vals + R_inf_vals + D_inf_vals),
                group=R0, color=as.factor(R0)), size = 1) +
  scale_color_manual(values=blues, guide = guide_legend(reverse = TRUE)) +
  theme_cowplot() +
  background_grid(minor = c("xy")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") + 
  labs(x="Proportion of Moderates Treated", y="Percent of Population Infected (%)", color=expression(R[0]))
#dev.off()
#labs(title=paste("Final Percent of Population Infected\nvs. Proportion of Moderates Treated\nfor various R0 values"),
#     x="Proportion of Moderates Treated", y="Percent of Population Infected (%)", color="R0")

# R_abx_inf vs. proportion treated (percent)
#pdf("treated-prop.pdf", width=6, height=6)



plot.B <- ggplot() +
  geom_line(data = df, aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size = 1) +
  scale_color_manual(values=blues, guide = guide_legend(reverse = TRUE)) +
  theme_cowplot() +
  geom_point(data = thresh_dat, aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size =2.5) + 
  background_grid(minor = c("xy")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")+
  labs(x="Proportion of Moderates Treated", y="Percent of Population Treated (%)", color=expression(R[0]))
#dev.off()
#labs(title="Final Percent of Population Treated\nvs. Proportion of Moderates Treated\nfor various R0 values",
#     x="Proportion of Moderates Treated", y="Percent of Population Treated (%)", color=expression(R[0]))

plot.legend <- ggplot() +
  geom_line(data = df, aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size = 1) +
  scale_color_manual(values=blues, guide = guide_legend(reverse = TRUE)) +
  theme_cowplot() +
  geom_point(data = thresh_dat, aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size =2) + 
  background_grid(minor = c("xy")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Proportion of Moderates Treated", y="Percent of Population Treated (%)", color=expression(R[0]))

legend<- cowplot::get_plot_component(plot.legend, "guide-box", return_all = TRUE)[[1]] 

plot_grid(plot.A, plot.B, legend, nrow = 1, rel_widths = c(1,1,0.15), labels = c("A", "B", ""))

# Attempts with other color palettes (Lindsay and Mac experimenting)

# Solve R_opt to check #########################################################
R_opt <- (p_A*nu_A/gamma_A +
              (p_MU+p_MT)*nu_M*(1/alpha_M + nu_sh/(gamma_M+mu_M)) +
              p_SU*(1/alpha_S + nu_sh/(gamma_S+mu_S)) +
              p_ST*(1/tau + nu_abx/gamma_Sabx)
) / 
  (p_A*nu_A/gamma_A +
     p_MU*nu_M*(1/alpha_M + nu_sh/(gamma_M + mu_M)) +
     p_MT*nu_M*(1/(delta*tau) + nu_abx/gamma_Mabx) +
     p_SU*(1/alpha_S + nu_sh/(gamma_S + mu_S)) +
     p_ST*(1/tau + nu_abx/gamma_Sabx)
  )
R_opt # check that this is 1.31538


# Adding thresholds to plots ###################################################
# converts q to proportion treated
proportion_treated <- function(q) { q*delta*tau/((1-q)*alpha_M + q*delta*tau) }

q_thresh_func <- function(q) {
  R_eff <- R_eff_func(q)
  
  if(R_eff >= 1) {
    S_inf_func_q <- function(S_inf_q) { S_inf_q - exp(R_eff*(S_inf_q - 1)) }
    S_inf_q <- uniroot(S_inf_func_q, c(0.00001,0.9999))$root
  } else {
    S_inf_q <- 0.999
  }
  
  return(R_abx_inf_func_q(0) - R_abx_inf_func_q(q))
}

## Plotting threshold lines
# more R0 values
thresholds_df <- data.frame(R0 = seq(1.001, 2, 0.001),
                            prop_less_doses = NA,
                            prop_stop = NA,
                            r_abx_inf = NA,
                            prop = NA)
# Don't worry about the error messages below!
for(R_0 in thresholds_df$R0){
  # find proportion for less doses
  beta <- beta_func(R_0)
  
  S_inf_func_0 <- function(S_inf_0) { S_inf_0 - exp(R_0*(S_inf_0 - 1)) }
  S_inf_0 <- uniroot(S_inf_func_0, c(0.00001,0.9999))$root
  
  q_thresh <- NA
  try(q_thresh <- uniroot(q_thresh_func, c(0.00001,0.9999))$root)
  thresholds_df$prop_less_doses[thresholds_df$R0 == R_0] <- proportion_treated(q_thresh)
  
  # save R_abx_inf for plots
  try(thresholds_df$r_abx_inf[thresholds_df$R0 == R_0] <- R_abx_inf_func_q(q_thresh))
  
  # find proportion for stopping
  solve_stop <- function(q) { R_eff_func(q) - 1 }
  stop_q <- NA
  try(stop_q <- uniroot(solve_stop, c(0.00001,0.9999))$root)
  thresholds_df$prop_stop[thresholds_df$R0 == R_0] <- proportion_treated(stop_q)
}

# Solve R_dut to check #########################################################
rev_q_thresh_func <- function(R_0) {
  R_eff <- R_eff_func_qbeta(1, beta_func(R_0))
  
  if(R_eff >= 1) {
    S_inf_func_q <- function(S_inf_q) { S_inf_q - exp(R_eff*(S_inf_q - 1)) }
    S_inf_q <- uniroot(S_inf_func_q, c(0.00001,0.9999))$root
  } else {
    S_inf_q <- 0.999
  }
  
  return(R_abx_inf_func_qbeta(0, beta_func(R_0)) - R_abx_inf_func_qbeta(1, beta_func(R_0)))
}
R_dut <- uniroot(rev_q_thresh_func, c(1.3,1.4))$root
R_dut # check that this is 1.399602


# plot thresholds by R0 vs. proportion treated
# fix data set the silly way for ggplot to work
thresholds_df <- rbind(thresholds_df, thresholds_df)
thresholds_df$prop <- thresholds_df$prop_stop
thresholds_df$prop[1:nrow(thresholds_df)/2] <- thresholds_df$prop_less_doses[1:nrow(thresholds_df)/2]
thresholds_df$prop_type <- rep(c("Reduced Doses", "Transmission Stopped"), each = nrow(thresholds_df)/2)



# flipping it makes more sense for the discussion
#pdf("thresh_lines.pdf", width=6, height=6)
ggplot(thresholds_df) +
  geom_ribbon(aes(y = prop_less_doses, xmin = 1, xmax = R0), fill="#004F63FF", alpha=0.5) +
  geom_ribbon(aes(y = prop_stop, xmin = 1, xmax = R0), fill="white", alpha=1) +
  geom_ribbon(aes(y = prop_stop, xmin = 1, xmax = R0), fill="#EE8577FF", alpha=0.5) +
  geom_line(aes(x = R0, y = prop, color = prop_type), size = 2) +
  scale_color_manual(values = c("#004F63FF", "#EE8577FF"), labels = c("Dose Utilization\nThreshold\n", "Outbreak Prevention\nThreshold")) +
  theme_cowplot() +
  geom_vline( aes(xintercept = R_opt), color = "#EE8577FF", lty = "dashed", size = 1.1) +
  geom_vline( aes(xintercept = R_dut), color = "#004F63FF", lty = "dashed", size = 1.1) +
  xlim(1,1.45) +
  #background_grid(minor = c("xy")) +
  labs(y="Proportion of Moderates Treated", x=expression(R[0]), color="Threshold Type:") +
  theme(legend.position = c(0.45, 0.25)) +
  ylim(0,1)
# add in region labels
#dev.off()


## Plotting threshold Xs on top of the R0 vs. proportion treated plot ##########
# set R0 values
thresholds_df <- data.frame(R0 = R0s,
                            prop_less_doses = NA,
                            prop_stop = NA,
                            r_abx_inf = NA,
                            prop = NA)
for(R_0 in thresholds_df$R0){
  # find proportion for less doses
  beta <- beta_func(R_0)
  
  S_inf_func_0 <- function(S_inf_0) { S_inf_0 - exp(R_0*(S_inf_0 - 1)) }
  S_inf_0 <- uniroot(S_inf_func_0, c(0.00001,0.9999))$root
  
  q_thresh <- NA
  try(q_thresh <- uniroot(q_thresh_func, c(0.00001,0.9999))$root)
  thresholds_df$prop_less_doses[thresholds_df$R0 == R_0] <- proportion_treated(q_thresh)
  
  # save R_abx_inf for plots
  try(thresholds_df$r_abx_inf[thresholds_df$R0 == R_0] <- R_abx_inf_func_q(q_thresh))
  
  # find proportion for stopping
  solve_stop <- function(q) { R_eff_func(q) - 1 }
  stop_q <- NA
  try(stop_q <- uniroot(solve_stop, c(0.00001,0.9999))$root)
  thresholds_df$prop_stop[thresholds_df$R0 == R_0] <- proportion_treated(stop_q)
}

# plot thresholds by R0 vs. proportion treated
# fix data set the silly way for ggplot to work
thresholds_df <- rbind(thresholds_df, thresholds_df)
thresholds_df$prop <- thresholds_df$prop_stop
thresholds_df$prop[1:nrow(thresholds_df)/2] <- thresholds_df$prop_less_doses[1:nrow(thresholds_df)/2]
thresholds_df$prop_type <- rep(c("Reduced Doses", "Transmission Stopped"), each = nrow(thresholds_df)/2)

# plot thresholds on top of doses given vs. proportion treated
df$trunc_prop <- df$prop
for(r0 in R0s){
  if(r0 > 1.06){
    x <- thresholds_df$prop_less_doses[thresholds_df$R0 == r0]
    df$trunc_prop[df$R0 == r0 & df$prop < x] <- NA
    if(is.na(x)){
      x <- thresholds_df$prop_stop[thresholds_df$R0 == r0]
      df$trunc_prop[df$R0 == r0] <- NA
    }
  }
}
blues_plus <- c("#03254c","#0c2950","#20365c","#314368","#3e4e72","#4c5a7d","#576586",
                "#6e7a99","#6e7a99","#8691ac","#a9b3c9","#c8d2e2","#d9e2f0","#4CBB17")

# plot with Xs
#pdf("with_X_thresh.pdf", width=6, height=6)
ggplot(df) +
  geom_line(aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size = 1) +
  geom_line(aes(x = trunc_prop, y = 100 * R_abx_inf_vals,group=R0, color=as.factor(R0)),
            size = 2, linetype = "dotted", lineend = "round") +
  scale_color_manual(values=blues, guide = guide_legend(reverse = TRUE)) +
  geom_point(data = thresholds_df, aes(x = prop_less_doses, y = 100 * r_abx_inf), size = 1, color = "#3B9212", shape = 4, stroke = 2) +
  geom_point(data = thresholds_df, aes(x = prop_stop, y = 0), size = 1, color = "red", shape = 4, stroke = 2) +
  theme_cowplot() +
  background_grid(minor = c("xy")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Proportion of Moderates Treated", y="Percent of Population Treated (%)", color=expression(R[0]))
#dev.off()

# Xs and green lines
ggplot(df) +
  geom_line(aes(x = prop, y = 100 * R_abx_inf_vals, group=R0, color=as.factor(R0)), size = 1) +
  geom_line(aes(x = trunc_prop, y = 100 * R_abx_inf_vals,group=R0, color="#4CBB17"), size = 1) +
  scale_color_manual(values=blues_plus, labels = c(rev(R0s),"Lower\nDoses"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = thresholds_df, aes(x = prop_less_doses, y = 100 * r_abx_inf), size = 1, color = "#3B9212", shape = 4, stroke = 2) +
  geom_point(data = thresholds_df, aes(x = prop_stop, y = 0), size = 1, color = "red", shape = 4, stroke = 2) +
  theme_cowplot() +
  background_grid(minor = c("xy")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Proportion of Moderates Treated", y="Percent of Population Treated (%)", color=expression(R[0]))
