# Include independent censoring
# Censoring time is a mixture of exponential and truncated normal

# longitudinal measure are taken at the 0, 1/3, 2/3 of
# the minimum OS for all subjetcs

# alpha_X is shrinked 20 times, X > 0

library(splines)
library(MASS)
library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(cowplot)
library(truncnorm)
library(tidyr)
# return TTM_data as the generated dataset

# set.seed(123)

n_patient = 1000 # Number of individuals
n_covariate = 2 # Number of covariates
n_randomeffect = 1 # Number of random effects
length_time = 36 # total length in months

########################## Global variables

# typical situation W(t*) = c(1, t*) \beta: Quran P9.

# prespecified W(t)
W_time = function(t){
  return(c(1, t))
}

# death para
alpha0_X = c(1, 2, 3) / 4
eta_01 = 0.01
eta_02 = 0.05
eta_0A = 0.05
eta_0b = 0

# dropout para
alpha1_X = c(1, 2, 3)/3
eta_11 = 0.01
eta_12 = 0.04
eta_1A = 0.04
eta_1b = 0

# Global coefficients
psi_X = c(10, 20)
beta_f = 0
tau2 = 100

# Simplified Sigma matrix
Sigma = matrix(100, ncol = 1)
# Pre-compute constants for hazard functions
hazard = function(t, t_split = 15,
                  alpha_X, X, A, b,
                  eta_01, eta_02,
                  eta_A, eta_b) {
  # Vectorized ifelse and avoid unnecessary list conversion
  hazard_tmp = (ifelse(t > t_split, eta_02, eta_01)) *
    exp(sum(alpha_X * unlist(X)) + A * eta_A + b * eta_b)
  return(hazard_tmp)
}

cumulative_hazard = function(t, t_split = 15,
                             alpha_X, X, A, b,
                             eta_01, eta_02,
                             eta_A, eta_b) {
  # Pre-compute time component
  time_component = ifelse(t <= t_split,
                          eta_01 * t,
                          eta_01 * t_split + eta_02 * (t - t_split))

  cumu_hazard_tmp = time_component *
    exp(sum(alpha_X * unlist(X)) + A * eta_A + b * eta_b)

  return(cumu_hazard_tmp)
}

# # Optimized beta functions
# beta_0_function = function(t){
#   20 - 4 * log(0.5 * t + 1)
# }
#
# beta_A_function = function(t){
#   10 + 6 * log(2 * t + 2)
# }


# beta_0_function = function(t){
#   0.5*t
# }
#
# beta_A_function = function(t){
#   0.2*t
# }


beta_0_function = function(t){
  5 - 2 * exp(0.05 * t)
}

beta_A_function = function(t){
  10 + 3 * exp(0.03 * t)
}

Y_death = function(D){
  exp(0.1 * D)
}


#################### Local variables

# Newton-Raphson iteration with tolerance
newton_raphson = function(init_val, cum_haz_func, haz_func,
                          u, max_iter = 100, tol = 0.01) {
  time_val = init_val
  for (i in 1:max_iter) {
    cum_haz = cum_haz_func(time_val)
    haz = haz_func(time_val)

    # Avoid division by zero
    if (abs(haz) < 1e-10) break

    time_new = time_val - (cum_haz + log(u)) / haz

    if (abs(time_new - time_val) < tol) {
      return(time_new)
    }
    time_val = time_new
  }
  return(time_val)
}

# generate one patient - optimized
generate_one = function(id, n_covariate, Sigma,
                        alpha0_X, alpha1_X, psi_X,
                        p = 0.5, max_iter = 100) {

  # treatment assignment
  A = rbinom(1, 1, p)

  # covariate generation
  X = c(rbinom(1, 1, 0.5), abs(rnorm(1, sd = 2)))
  b = rnorm(1) # random intercept for coxme

  # death time generation
  death_u = runif(1)
  death_time = newton_raphson(
    abs(rnorm(1, sd = 10)),
    function(t) cumulative_hazard(t, t_split = 15, alpha0_X, c(1, X), A, b,
                                  eta_01, eta_02, eta_0A, eta_0b),
    function(t) hazard(t, t_split = 15, alpha0_X, c(1, X), A, b,
                       eta_01, eta_02, eta_0A, eta_0b),
    death_u
  )

  # drop-out time generation
  drop_u = runif(1)
  dropout_time = newton_raphson(
    abs(rnorm(1, sd = 10)),
    function(t) cumulative_hazard(t, t_split = 15, alpha1_X, c(1, X), A, b,
                                  eta_11, eta_12, eta_1A, eta_1b),
    function(t) hazard(t, t_split = 15, alpha1_X, c(1, X), A, b,
                       eta_11, eta_12, eta_1A, eta_1b),
    drop_u
  )

  # censoring
  censoring_time = ifelse(rbinom(1, 1, 0.5),
                          rexp(1, rate = 0.05),
                          rtruncnorm(1, a = 0, b = 50, mean = 40, sd = 5))

  # Overall Survival
  OS_time = min(death_time, dropout_time, censoring_time)
  death = (death_time < censoring_time) & (death_time < dropout_time)
  dropout = (dropout_time < censoring_time) & (dropout_time < death_time)
  event = ifelse(death == 1, 1, ifelse(dropout == 1, 2, 0))

  # random effects (including random intercept)
  # b_longi = mvrnorm(n = 1, mu = 0, Sigma = Sigma)
  b_longi = b
  randomslope = mvrnorm(n = 1, mu = 0, Sigma = Sigma)
  # Return as a named vector for better readability
  patient_record = c(id = id, A = A,
                     X1 = X[1], X2 = X[2], u = b, b_longi = b_longi, randomslope = randomslope,
                     OS_time = OS_time, event = event, death = death, dropout = dropout,
                     death_time = death_time, dropout_time = dropout_time,
                     censoring_time = censoring_time)

  return(patient_record)
}

# Vectorized data generation
data_list = vector("list", n_patient)
for (i in 1:n_patient) {
  data_list[[i]] = generate_one(id = i, n_covariate, Sigma,
                                alpha0_X, alpha1_X, psi_X)
}
TTM_data = do.call(rbind, data_list)

# Convert to data.frame for easier manipulation
TTM_data = as.data.frame(TTM_data)

# # longi measures are taken at 0, 1/3, 2/3 of minimum OS_time
# tpros = min(TTM_data$OS_time) * (0:2) / 3
# tpros = expand.grid(id = 1:n_patient, t_pros = tpros)
# TTM_data = merge(TTM_data, tpros, by = "id")

# Merge back to create longitudinal dataset
TTM_data = TTM_data %>%
  rowwise() %>%
  mutate(
    t_pros = list(
      # if (OS_time <= 10) {
      seq(0, OS_time, by = 2)
      # } else {
      #   unique(c(0:10, seq(10, OS_time, by = 2), OS_time))
      # }
    )
  ) %>%
  ungroup() %>%
  unnest(t_pros) %>%
  as.data.frame()


# Vectorized longitudinal outcome calculation
t_retro = TTM_data$death_time - TTM_data$t_pros
Y = beta_0_function(t_retro) +
  TTM_data$A * beta_A_function(t_retro) +
  (psi_X[1] * TTM_data$X1 + psi_X[2] * TTM_data$X2) +
  # 5* TTM_data$b_longi +
  # 5* TTM_data$A *TTM_data$randomslope +
  rnorm(nrow(TTM_data), sd = sqrt(tau2))

TTM_data = cbind(TTM_data, t_retro, Y)

data_unique_event = unique(TTM_data[, c("id", "death", "dropout")])
mean(data_unique_event$death)
mean(data_unique_event$dropout)
