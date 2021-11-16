#### SIMULATE HYPOTHETICAL SPILLOVERS IN AN EXPERIMENT ON A NETWORK ###
#### WARNING: THIS TAKES ABOUT 10 HOURS TO SOURCE IN A DEDICATED COMPUTE CLUSTER ####
#### YOU DO NOT NEED TO SOURCE THIS FILE TO REPRODUCE THE MAIN FINDINGS ####
#### Creates sim_df.rda
## Simulate a network either by distance
## Create hypothetical spillover effect within an upper bound
## Use supervised learning model selection protocol to estimate upper bound
## Evaluation metrics should be
## Bias
## Consistency
## Power
## Estimated upper bound

#### SETUP ####
## Packages
library(tidyverse) # utilities
library(sna) # social network data analysis
library(DeclareDesign) # design-based estimators
# devtools::install_github('szonszein/interference')
library(interference)
library(caret) # supervised learning pipeline
library(data.table) # data manipulation utilities
library(future) # parallel processing
library(future.apply) # parallel apply
library(here) # relative path management

## Load road adjacency network of electoral areas in Ghana
## From Bowers et al (2018 Social Networks)
net = as.matrix(read.csv(here("data/GHAadjmatrix.csv"),row.names=1))

# Remove row and col names so matrix is recognized as symmetric
rownames(net) = NULL
colnames(net) = NULL

# Set diagonal to zero so observations are not considered neighbors to themselves
diag(net) = 0

# calculate geodesic
net_dist = geodist(net)$gdist

# node degree
deg = degree(net)

#### PROTOCOL ####
## We want to evaluate the properties of an algorithm that guesses the optimal spillover upper bound
## Then estimates the average spillover effect
## Repeat this many times with many simulations and parameter values

#### COMPUTE EXPOSURE PROBABILITIES ####
make_prob_exposure = function(dist, k, alpha, reps = 1000){
  # observations
  n = nrow(dist)
  
  # transform distance matrix into adjacency matrix
  adj = ifelse(dist <= k & dist != 0, 1, 0)
  
  # compute treatment vectors
  pot_trt_vector = make_tr_vec_bernoulli(adj, p = alpha, R = reps, cluster = "no")
  
  # compute probability of exposure mapping
  prob_exp_mat = make_exposure_prob(pot_trt_vector, adj, make_exposure_map_AS, list(hop = 1))
  
  # isolate indirect exposure probability
  prob_exp = make_prob_exposure_cond(prob_exp_mat)["ind1", ]
  
  out = data.frame(cbind(a = rep(alpha, n), k = rep(k, n), prob_exp))
  
  return(out)
}

# test: ok
# make_prob_exposure(dist = net_dist, k = 2, alpha = 0.4)

# grid of potential alpha-k combinations
designs = expand.grid(
  alpha = seq(.1, .4, by =.1),
  k = 1:4
)

# compute corresponding exposure probs and make a data frame
exp_probs = NULL

set.seed(20201110)

for(i in 1:nrow(designs)){
  temp = make_prob_exposure(
    dist = net_dist,
    k = designs[i, "k"],
    alpha = designs[i, "alpha"]
  )
  
  exp_probs = rbind(exp_probs, temp)
}

# check
# with(exp_probs, table(a, k))

#### FUNCTION TO ESTIMATE AVERAGE SPILLOVER EFFECT UP TO K ####
## Compute Horvitz-Thompson estimate based on k:
## Use lm so we have a predict class in caret
## formula: outcome ~ treatment
## k: spillover upper bound in network distance units
## exposure: data frame with k-specific binary exposure indicators as columns
## weights: data frame with k-specific inverse probability of exposure weights
## data: data frame containing outcome and treatment
lmk = function(formula, k, exposure, weights, data) {
  # extract z
  z = model.frame(formula, data)[, 2]

  # extract k-specific exposure
  data$exposure = exposure[, k]
  
  # extract k-specific weights from data vector
  data$ipw =  weights[, k]
  
  # exclude treated units
  new_data = data[which(z == 0), ] 
  
  # estimate with lm (could have also just been difference in means or other method)
  y = model.frame(formula, new_data)[, 1]
  
  exposure = new_data[, "exposure"]
  
  ipw = new_data[, "ipw"]
  
  model = lm_robust(y ~ exposure, weights = ipw)
  
  return(model)
}

# test: ok
# test_lm = lmk(y_obs ~ z, k = 2, exposure = exposure, data = dat, weights = weights)

#### CREATE CARET FUNCTION ####
## Based on a given value of k
## compute RMSE of lm using caret
## parameters similar to lmk
caret_lmk = function(formula, k, exposure, weights, nfolds = 5, data){
  # extract z
  z = model.frame(formula, data)[, 2]

  # extract k-specific exposure
  data$exposure = exposure[, k]
  
  # extract k-specific weights from data vector
  data$ipw =  weights[, k]
  
  # exclude treated units
  new_data = data[which(z == 0), ]
  
  # isolate y, exposure, weights
  y = model.frame(formula, new_data)[,1]
  
  exposure = new_data[, "exposure"]
  
  ipw = new_data[, "ipw"]
  
  tune_data = cbind(y, exposure, ipw)
  
  # caret model
  mod = train(
    form = y ~ exposure,
    data = tune_data,
    method = "lm",
    trControl = trainControl(method = "cv", number = nfolds),
    weights = ipw
  )
  
  # return results
  out = mod$results %>% mutate(k = k)
  
  return(out)
}

# test: ok
# caret_lmk(y_obs ~ z, k = 2, exposure = exposure, weights = weights, data = dat)

#### SIMULATION ####
## FUNCTION
# Tune lmk_caret and then estimate lmk with best_k
# Compare with oracle that knows true k
sim_spill = function(alpha, gamma, tau, k, sim = 1, nfolds = 5, reps = 1000){
  ## SETUP
  # net, net_dist, and deg defined outside of function
  # observations
  n = nrow(net_dist)

  # baseline outcome
  y00 = runif(n)
  
  # outcome under treatment
  y10 = y00*tau
  
  # outcome under infection
  y01 = y00*tau
  
  # assign treatment
  z = rbinom(n, 1, alpha)
  
  # indicate whether in neighborhood of treated unit
  # This is also our "true" exposure mapping
  exp_true = ifelse(rowSums(net_dist[, which(z == 1)] <= k & net_dist[, which(z == 1)] != 0) > 0, 1, 0)
  
  # Depending on k and infection probability gamma, assign infection
  inf = ifelse(exp_true == 1, rbinom(1, 1, gamma), 0)
  
  # Realize observed outcome
  y_obs = ifelse(z == 1, y10, ifelse(inf == 1, y01, y00))
  
  # merge
  dat = data.frame(y00, y10, y01, y_obs, z, inf, exp_true, deg)
  
  # compute average treated node degree
  avg_trt_deg = dat %>% filter(z == 1) %>% dplyr::summarize(mean_deg = mean(deg)) %>% .$mean_deg
  
  # and degree/treatment correlation
  dz_cor = with(dat, cor(deg, z))
  
  # Estimand is average exposure effect
  estimand = dat %>% filter(z == 0) %>%  dplyr::summarize(est = mean(y01) - mean(y00)) %>% .$est
  
  ## Compute 4 exposure mappings
  # transform distance matrix into adjacency matrix
  adj_1 = ifelse(net_dist <= 1 & net_dist != 0, 1, 0)
  adj_2 = ifelse(net_dist <= 2 & net_dist != 0, 1, 0)
  adj_3 = ifelse(net_dist <= 3 & net_dist != 0, 1, 0)
  adj_4 = ifelse(net_dist <= 4 & net_dist != 0, 1, 0)
  
  # observed indirect exposure
  exp_1 = make_exposure_map_AS(adj_1, z, hop = 1)[, "ind1"]
  exp_2 = make_exposure_map_AS(adj_2, z, hop = 1)[, "ind1"]
  exp_3 = make_exposure_map_AS(adj_3, z, hop = 1)[, "ind1"]
  exp_4 = make_exposure_map_AS(adj_4, z, hop = 1)[, "ind1"]
  
  exposure = cbind(exp_1, exp_2, exp_3, exp_4)
  
  # fetch exposure probs
  prob_exp_1 = exp_probs %>% filter(a == alpha & k == 1) %>% .$prob_exp
  prob_exp_2 = exp_probs %>% filter(a == alpha & k == 2) %>% .$prob_exp
  prob_exp_3 = exp_probs %>% filter(a == alpha & k == 3) %>% .$prob_exp
  prob_exp_4 = exp_probs %>% filter(a == alpha & k == 4) %>% .$prob_exp
  
  # Compute inverse weights
  wt_1 = exp_1/prob_exp_1 + (1 - exp_1)/(1 - prob_exp_1)
  wt_2 = exp_1/prob_exp_2 + (1 - exp_2)/(1 - prob_exp_2)
  wt_3 = exp_1/prob_exp_3 + (1 - exp_3)/(1 - prob_exp_3)
  wt_4 = exp_1/prob_exp_4 + (1 - exp_4)/(1 - prob_exp_4)
  
  weights = cbind(wt_1, wt_2, wt_3, wt_4)
  
  ## MODEL SELECTION
  # tune caret lmk
  kgrid = 1:4
  
  tune = rbindlist(lapply(kgrid, caret_lmk, 
                          formula = y_obs ~ z, 
                          exposure = exposure,
                          weights = weights,
                          data = dat,
                          nfolds = nfolds))
  
  # select k with largest RMSE within 1SE of smallest RMSE to avoid overfitting
  min_plus_1se = sum(tune[which(tune$RMSE == min(tune$RMSE)), c("RMSE", "RMSESD")])
  
  candidates = tune[which(tune$RMSE <= min_plus_1se), ]
  
  best_k = as.numeric(candidates[which(candidates$RMSE == max(candidates$RMSE)), "k"])
  
  # estimate
  algorithm = tidy(lmk(y_obs ~ z, k = best_k, exposure = exposure, weights = weights, data = dat)) %>% 
    filter(term != "(Intercept)") %>% 
    mutate(model = "algorithm", k_hat = best_k)
  
  oracle = tidy(lmk(y_obs ~ z, k = k, exposure = exposure, weights = weights, data = dat)) %>% 
    filter(term != "(Intercept)") %>% 
    mutate(model = "oracle", k_hat = k)
  
  ## OUTPUT
  out = rbind(algorithm, oracle)
  
  out$avg_deg = avg_trt_deg
  out$dz_cor = dz_cor
  out$alpha = alpha
  out$gamma = gamma
  out$tau = tau
  out$k = k
  out$estimand = estimand
  out$sim = sim
  
  # reorder
  out = out %>% 
    select(model, estimand, estimate, std.error, statistic, p.value, conf.low, conf.high, df,
           k, k_hat, alpha, gamma, tau, dz_cor, avg_deg, sim)
    
  return(out)
}

# test: ok
# sim_spill(alpha = 0.1, gamma = 1, tau = 0.26, k = 2)

## PARAMETER GRID

# Original parameters in Bowers et al 2018
# alphas <- seq(.05,.5,by=.05)
# taus <- (1-.37*c(1,2))
# Maximum time-to-infection
# k <- 1

# New parameters
# alpha = seq(.1, .4, by = .1)
# infection probability if exposed gamma = seq(.5, 1, by = .1)
# k = 1:2
params = expand.grid(
  alpha = seq(.1, .4, by = .1),
  gamma = seq(.5, 1, by = .1),
  tau = c(.26, .63), # roughly 1 and 2 sds in uniform distribution 
  k = 1:3,
  sim = 1:1000
)

## SIMULATE AND EXPORT RESULT
plan(multicore)

set.seed(20201108)

sim_list = future_Map(sim_spill, 
                      params$alpha,
                      params$gamma,
                      params$tau,
                      params$k,
                      params$sim,
                      future.seed = TRUE)

plan(sequential)

sim_df = data.frame(rbindlist(sim_list))

save(sim_df, file = here("data/sim_df.rda"))
