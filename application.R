#### REANALYZE EXPERIMENT BY ICHINO AND SCHUNDELN (2012) ####
#### EFFECT OF ELECTION OBSERVERS ON VOTER REGISTRATION IRREGULARITIES ####
#### https://doi.org/10.1017/S0022381611001368 ####

### SETUP ####
library(here) # relative path management
library(tidyverse) # Useful stuff
library(data.table) # data manipulation tools
library(broom) #tidy glmnet output
library(DeclareDesign) # Design-based features
library(glmnet) # lasso
library(sparsereg) # LASSOplus
library(MESS) # compute adaptive lasso weights
library(reshape2) # melt function


# ggplot global options
theme_set(theme_bw(base_size = 20))

# Tables
library(kableExtra)

#### LOAD DATA ####
# Experiment data
load(here("IchinoSchuendeln_spilloversJOP.RData"))

ghana = x

# Distances between ELAs
load(here("IchinoSchuendeln_spilloversJOP_ELAdistances.RData"))

distances = x

rm(x)

#### PREPROCESSING ####
# Subset to control only, this gives a sample of 791 observations
ghana_ctl = ghana %>% filter(Tela == 0)


#### RESEARCH DESIGN TABLE ####
tab = ghana %>% 
  group_by(Tcon, Tela) %>% 
  tally() %>% 
  ungroup()

tab = tab %>% mutate(Tcon = ifelse(Tcon == 1, "Treatment", "Control"),
                     Tela = ifelse(Tela == 1, "Treatment", "Control"))

colnames(tab) = c("Constituencies (clusters)", "Electoral areas (ELAs)", "Observations")

tab %>% 
  kbl(caption = "Research design in Ichino and Schundeln (2012)", booktabs = TRUE) %>% 
  save_kable(here("tab1.pdf"))
  

#### REPRODUCE ORIGINAL FINDINGS ####
lm_0 = lm_robust(percchangeregELA0804 ~ assignedTin5C + assignedTin0510C + Tcon, 
                 data = ghana_ctl, 
                 clusters = block, 
                 fixed_effects = block,
                 se_type = "stata")

### CREATE DISTANCE MATRIX ###
dist = matrix(nrow = length(ghana$ELA_GIS), ncol = length(ghana$ELA_GIS))

rownames(dist) = unique(distances$ELA_GIS)

colnames(dist) = unique(distances$ELA_GIS)

rows = distances$ELA_GIS

cols = distances$ELA_GIS_other

dist[which(rownames(dist) %in% rows[1]), which(colnames(dist) %in% cols[1])]

find_dist = function(i){
  d = distances[distances$ELA_GIS == rows[i] & distances$ELA_GIS_other == cols[i], "distance"]
  
  return(d)
}

for(i in 1:length(rows)){
  dist[which(rownames(dist) %in% rows[i]), which(colnames(dist) %in% cols[i])] = find_dist(i)
}

diag(dist) = 0

dist = ifelse(is.na(dist), Inf, dist)

# Subset matrix to control vs treated
trt_ids = ghana[ghana$Tela == 1,]$ELA_GIS

trt_dist = dist[which(!(rownames(dist) %in% trt_ids)), which(colnames(dist) %in% trt_ids)]

#### CREATE EXPOSURE PREDICTORS AT 5KM INTERVALS ####
make_exposure = function(net, k){
  out = rowSums(ifelse(net > k & net <= k + 5, 1, 0))
  return(out)
}

exposure = NULL 

k = seq(0, 45, by = 5)

for(i in 1:length(k)){
  temp = make_exposure(trt_dist, k = k[i])
  exposure = cbind(exposure, temp)
}

# check
dim(exposure)

colnames(exposure) = c("exp_00_05", "exp_05_10", "exp_10_15", "exp_15_20", "exp_20_25",
                       "exp_25_30", "exp_30_35", "exp_35_40", "exp_40_45", "exp_45_50")

exposure = data.frame(exposure)

# Merge into data
ghana_exp = cbind(ghana_ctl, exposure)

#### FIT AND TUNE LASSO ####
# Make subset of outcome and exposure variables
ghana_exp = ghana_exp %>% dplyr::select(percchangeregELA0804, exp_00_05, exp_05_10, exp_10_15,
                                        exp_15_20, exp_20_25, exp_25_30, exp_30_35, exp_35_40,
                                        exp_40_45, exp_45_50)

g_y = ghana_exp$percchangeregELA0804

g_x = as.matrix(ghana_exp[, -1])

# Lasso
# alpha = 1 is lasso, alpha = 0 is ridge. Anything in between is elastic net.
# nfolds = 10 stands for 10-fold cross-validation
set.seed(20210301)
lasso_cv = cv.glmnet(y = g_y, x = g_x, alpha = 1, nfolds = 10)

lasso_fit = glmnet(y = g_y, x = g_x, alpha = 1)

coef(lasso_cv, s = "lambda.min") # 0-5, 10-15, 20-25, 30-35, 40-45

coef(lasso_cv, s = "lambda.1se") # nothing, use range

# Visualize
# Create fit tidy output
lasso_df = tidy(lasso_fit) %>% 
  filter(term != "(Intercept)")

# Plot
lasso_plot = ggplot(lasso_df) +
  aes(x = lambda, y = estimate, color = term) +
  labs(x = expression(lambda), y = "Coefficient", color = "Predictor") +
  geom_hline(yintercept = 0) +
  geom_line(size = 1) +
  annotate("rect", xmin = lasso_cv$lambda.min, xmax = lasso_cv$lambda.1se,
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#0072B2") +
  geom_vline(xintercept = c(lasso_cv$lambda.min, lasso_cv$lambda.1se), size = 1,
             colour = "#0072B2", linetype = "dashed") +
  scale_color_viridis_d(labels = c("0-5 km", "5-10 km", "10-15 km", "15-20 km","20-25 km",
                                   "25-30 km", "30-35 km", "35-40 km", "40-45 km", 
                                   "45-50 km")) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# save
ggsave(plot = lasso_plot, filename = here("figA1.pdf"), height = 6, width = 12)

#### ADAPTIVE LASSO ####
# Calculate weights
# From MESS documentation:
# The weights returned are 1/abs(beta_hat)^nu where the beta-parameters are
# estimated from the corresponding linear model 
# (either multivariate or univariate)
lwt = adaptive.weights(y = g_y, x = g_x, weight.method = "multivariate")

# Fit
set.seed(20210301)
adapt_cv = cv.glmnet(y = g_y, x = g_x, alpha = 1, nfolds = 10, penalty.factor = lwt$weights)

adapt_fit = glmnet(y = g_y, x = g_x, alpha = 1, penalty.factor = lwt$weights)

# check coefs
coef(adapt_cv, s = "lambda.min") # 0-5, 10-15, 20-25, 30-35, 40-45.

coef(adapt_cv, s = "lambda.1se") 

# Visualize
# Create fit tidy output
adapt_df = tidy(adapt_fit) %>% 
  filter(term != "(Intercept)")

# Plot
adapt_plot = ggplot(adapt_df) +
  aes(x = lambda, y = estimate, color = term) +
  labs(x = expression(lambda), y = "Coefficient", color = "Predictor") +
  geom_hline(yintercept = 0) +
  geom_line(size = 1) +
  annotate("rect", xmin = adapt_cv$lambda.min, xmax = adapt_cv$lambda.1se,
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#0072B2") +
  geom_vline(xintercept = c(adapt_cv$lambda.min, adapt_cv$lambda.1se), size = 1,
             colour = "#0072B2", linetype = "dashed") +
  scale_colour_viridis_d(labels = c("0-5 km", "10-15 km","20-25 km",
                                    "30-35 km", "35-40 km", "40-45 km", 
                                    "45-50 km")) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# save
ggsave(plot = adapt_plot, filename = here("figA2.pdf"), height = 6, width = 12)

#### LASSO + OLS ####
# Based on Belloni and Chernozhukov (2013). Point is to use a goodness of fit threshold to reduce bias from irrelevant effects. 
# Ratkovic and Tingley (2017) call this Lasso + OLS. 
# Specifically, we are interested in what the original calls OLS post-fit Lasso. The procedure is the following:

   
# 1. Fit and tune lasso to select an initial set of effects (coefficients)
# 2. Compute residual variance of selected lasso fit
# 3. Fit OLS in all subsets of the initial fit (in this case, in all relevant subsets, as some combinations may not make sense in this application)
# 4. Select OLS fit with residual variance close to initial fit residual variance

# lasso optimal lambda fit
lasso_lmin = glmnet(y = g_y, x = g_x, alpha = 1, lambda = lasso_cv$lambda.min)

g_yhat = predict(lasso_lmin, newx = g_x)

lmin_res = g_y - g_yhat

qb = var(lmin_res)

# OLS implied by selected lasso
coef(lasso_cv, s = "lambda.min")

post_ols = lm(percchangeregELA0804 ~ exp_00_05 + exp_10_15 + exp_20_25 + exp_30_35 + exp_40_45,
              data = ghana_exp)

qb0 = var(post_ols$residuals)

# Compute gamma threshold
gamma = (qb0 - qb)/2

# make a grid of switches for each variable in model
reg_mat = expand.grid(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE), 
                      c(TRUE, FALSE))

# omit last row since it implies empty model
reg_mat = reg_mat[-32, ]

names(reg_mat) = c("exp_00_05", "exp_10_15", "exp_20_25", "exp_30_35", "exp_40_45")

regressors = c("exp_00_05", "exp_10_15", "exp_20_25", "exp_30_35", "exp_40_45")

model_list = apply(reg_mat, 1, function(x) as.formula(paste(c("percchangeregELA0804 ~ 1",
                                                              regressors[x]),
                                                            collapse = " + ")))

# Get results
post_ols_sigma = sapply(model_list, function(x) var(residuals(lm(x, data = ghana_exp))))

t = post_ols_sigma - qb

t < as.numeric(gamma)

reg_mat[which(t < as.numeric(gamma)), ]

# And it turns out that it selects the same model as lasso and adaptive lasso.

#### LASSOPLUS ####
# Based on Ratkovic and Tingley (2017). 
# Advantage is that simultaneously performs model selection, estimation, and produces a measure of uncertainty.
# Note the use of a different data set
covs = as.matrix(ghana_ctl[, c("block", "Tcon")])

blocks = ghana_ctl$block

set.seed(20210302)
lplus = sparsereg(y = ghana_ctl$percchangeregELA0804, X = g_x, treat = covs, id = blocks)

# These are the variables as they appear in order
colnames(lplus$X)

# This is the density distribution of the estimates
plot(density(lplus$beta.ci[, 1]))

# useful output that I will use to visualize results later
lplus_out = summary(lplus$beta.ci)

#### VISUALIZE RELATIONSHIP BETWEEN VARIABLES ####
cormat = round(cor(ghana_exp), 2)

colnames(cormat) = c("Outcome", "0-5\nkm", "5-10\nkm", "10-15\nkm", "15-20\nkm", "20-25\nkm", 
                     "25-30\nkm", "30-35\nkm", "35-40\nkm", "40-45\nkm", "45-50\nkm")

rownames(cormat) = c("Outcome", "0-5 km", "5-10 km", "10-15 km", "15-20 km", "20-25 km", 
                     "25-30 km", "30-35 km", "35-40 km", "40-45 km", "45-50 km")

cormat[upper.tri(cormat)] = NA  

melted_cormat = melt(cormat, na.rm = TRUE)

corplot = ggplot(melted_cormat) +
  aes(x = Var2, y = Var1, fill = value) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient2(mid = "white", low = "gray20", high = "gray20") +
  scale_y_discrete(limits = rev(levels(melted_cormat$Var1))) +
  theme(legend.position = "none") +
  labs(x = "", y = "")

# save
ggsave(plot = corplot, filename = here("fig1.pdf"), height = 6, width = 12)

#### VISUALIZE TUNING + ESTIMATES ###
# A reminder that I already did the original model in Ichino and Schundeln (2012)
summary(lm_0)

# do the lasso-related ones
# For this we need to adjust the confidence interval for false coverage-statement rate (FCR)
# See Benjamini and Yekutieli (2005): https://doi.org/10.1198/016214504000001907

# default significance level 
alpha = 0.05

# Number of selected parameter
R = 5

# number of candidate parameters
m = 10

# FCR alpha
fcr_alpha = alpha*R/m

fcr_alpha

1-fcr_alpha # 97.5% CI

ghana_lm = cbind(ghana_ctl, exposure)

lm_1 = lm_robust(percchangeregELA0804 ~ exp_00_05 + exp_10_15 + exp_20_25 + exp_30_35 +
                   exp_40_45 + Tcon, data = ghana_lm, 
                 clusters = block, 
                 fixed_effects = block,
                 se_type = "stata")

fcr_adjust_ci = as.numeric(confint(lm_1, level = 1-fcr_alpha))

# For LASSOplus we only need to compute relevant values from posterior distribution
lplus_out$quantiles[1,]

tidy(lm_1)

lm_2 = data.frame(term = "exp_00_05", 
                  estimate = lplus_out$statistics[1, 1], 
                  std.error = lplus_out$statistics[1, 3],
                  statistic = NA,
                  p.value = NA,
                  conf.low = lplus_out$quantiles[1, 1],
                  conf.high = lplus_out$quantiles[1, 5],
                  df = NA,
                  outcome = "percchangeregELA0804",
                  model = "LASSOplus")

# Merge and replace FCR adjusted CIs where appropriate
estimates = rbind(
  tidy(lm_0) %>% mutate(model = "Original"),
  tidy(lm_1) %>% mutate(model = "Lasso", 
                        conf.low = fcr_adjust_ci[1:6], conf.high = fcr_adjust_ci[7:12]),
  tidy(lm_1) %>% mutate(model = "Adaptive lasso", 
                        conf.low = fcr_adjust_ci[1:6], conf.high = fcr_adjust_ci[7:12]),
  tidy(lm_1) %>% mutate(model = "Lasso + OLS", 
                        conf.low = fcr_adjust_ci[1:6], conf.high = fcr_adjust_ci[7:12]),
  lm_2
)

estimates = estimates %>% 
  filter(term != "Tcon") %>% 
  mutate(term = ifelse(term == "assignedTin5C", "exp_00_05",
                       ifelse(term == "assignedTin0510C", "exp_05_10", term)))

estimates$model = fct_relevel(estimates$model, "Original", "Lasso", 
                              "Adaptive lasso", "Lasso + OLS", "LASSOplus")

#### VISUALIZE ####
est_plot = ggplot(estimates) +
  aes(x = reorder(term, desc(term)), y = estimate) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_linerange(aes(x = term, ymin = conf.low, ymax = conf.high),
                 size = 1) +
  facet_grid(model ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = c("exp_00_05" = "0-5 km",
                              "exp_05_10" = "5-10 km",
                              "exp_10_15" = "10-15 km", 
                              "exp_20_25" = "20-25 km",
                              "exp_30_35" = "30-35 km",
                              "exp_40_45" = "40-45 km")) +
  labs(x = "Treated ELAs in range", y = "Estimate") +
  coord_flip() +
  theme(strip.text.y = element_text(angle = 0))

# save
ggsave(plot = est_plot, filename = here("fig2.pdf"), height = 12, width = 12)
