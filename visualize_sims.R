#### OVERVIEW ####
# This document visualizes the results of the simulations created with `make_sims.R` through a compute cluster.
# The output is stored on `data/sim_df.rda`.

# Following a road connectivity network of electoral areas in Ghana from:
  
#  Bowers, Jake, Bruce A. Desmarais, Mark Frederickson, Nahomi Ichino, Hsuan-Wei Lee, and Simi Wang. 2018. "Models, methods, and network topology: Experimental design for the study of interference." *Social Networks* 54: 196-208

#### SETUP ####
library(tidyverse) # utils
library(GGally) # plot network
library(sna) # social network analysis
library(interference) # causal inference interference
library(here) # relative path management

## ggplot2 global options
theme_set(theme_bw(base_size = 20))

#### DATA ####
load(here("sim_df.rda"))

# I forgot to capitalize model names!
sim_df = sim_df %>% mutate(model = ifelse(model == "algorithm", "Protocol", "Oracle"))

## Load road adjacency network of electoral areas in Ghana
## From Bowers et al (2018 Social Networks)
net = as.matrix(read.csv(here("GHAadjmatrix.csv"),row.names=1))

# Remove row and col names so matrix is recognized as symmetric
rownames(net) = NULL
colnames(net) = NULL

# Set diagonal to zero so observations are not considered neighbors to themselves
diag(net) = 0

# sample size
nrow(net)

# node degree
deg = degree(net)

summary(deg)

sd(deg)



#### SAMPLE NETWORK EXPERIMENT ####
# Vector of treatment assignment
set.seed(20210330)
z = rbinom(n = nrow(net), size = 1, prob = 0.1)

# Visualize
net = network(net)

net %v% "condition" = ifelse(z == 1, "Treatment", "Control")

sample_net = ggnet2(net, node.size = 4, color = "condition", 
                    palette = c("Treatment" = "black", "Control" = "gray"),
                    color.legend = "Condition", legend.position = "bottom", 
                    legend.size = 20)

# save
ggsave(plot = sample_net, filename = here("fig3.pdf"), height = 6, width = 12)

#### COMPUTE PERFORMANCE METRICS ####
# Filter funky results (see explanation in paper)
sim_df_short = sim_df %>%
  filter(estimate < 1 & estimate > -1)

# Performance metrics
performance = sim_df_short %>% 
  group_by(model, kappa = k, alpha, gamma, tau) %>% 
  summarize(bias = mean(estimate - estimand, na.rm = TRUE),
            mse = mean((estimate - estimand)^2, na.rm = TRUE),
            mad = mean(abs(estimate - estimand), na.rm = TRUE),
            power = mean(p.value < 0.05, na.rm = TRUE),
            precision = 1/var(estimate, na.rm = TRUE),
            k_hat_mean = mean(k_hat, na.rm = TRUE),
            k_hat_lo = as.numeric(quantile(k_hat, 0.025)),
            k_hat_hi = as.numeric(quantile(k_hat, 0.975)),
            deg_mean = mean(avg_deg),
            deg_lo = as.numeric(quantile(avg_deg, 0.025)),
            deg_hi = as.numeric(quantile(avg_deg, 0.975)))

# create labels
performance = performance %>% 
  mutate(alpha = recode_factor(as.factor(alpha),
                               "0.1" = "alpha == 0.1",
                               "0.2" = "alpha == 0.2",
                               "0.3" = "alpha == 0.3",
                               "0.4" = "alpha == 0.4")) %>%
  mutate(gamma = recode_factor(as.factor(gamma),
                               "0.5" = "gamma == 0.5",
                               "0.6" = "gamma == 0.6",
                               "0.7" = "gamma == 0.7",
                               "0.8" = "gamma == 0.8",
                               "0.9" = "gamma == 0.9",
                               "1" = "gamma == 1"))

performance$model = fct_relevel(performance$model, "Protocol", "Oracle")

#### VISUALIZE MAIN RESULTS ####
bias_plot = ggplot(performance %>% 
                     filter(tau == 0.26) %>% 
                     filter(alpha == "alpha == 0.1" | alpha == "alpha == 0.4") %>% 
                     filter(gamma == "gamma == 0.5" | gamma == "gamma == 1")) +
  aes(x = as.factor(kappa), y = bias, shape = model, size = model) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Bias", shape = "Model", size = "Model")

kbar_plot = ggplot(performance %>% filter(tau == 0.26) %>% 
                     filter(alpha == "alpha == 0.1" | alpha == "alpha == 0.4") %>% 
                     filter(gamma == "gamma == 0.5" | gamma == "gamma == 1")) +
  aes(x = as.factor(kappa), y = k_hat_mean, shape = model, size = model) +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = expression(paste('Mean selected upper bound (',bar(k),')')),
       shape = "Model", size = "Model") +
  ylim(0,4)

# save
ggsave(plot = bias_plot, filename = here("fig4.pdf"), height = 6, width = 12)
ggsave(plot = kbar_plot, filename = here("fig5.pdf"), height = 6, width = 12)


#### EXTRA VIZ FOR APPENDIX ####

# Bias
bias_026 = ggplot(performance %>% filter(tau == 0.26)) +
  aes(x = as.factor(kappa), y = bias, shape = model, size = model) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Bias", shape = "Model", size = "Model")

bias_063 = ggplot(performance %>% filter(tau == 0.63)) +
  aes(x = as.factor(kappa), y = bias, shape = model, size = model) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Bias", shape = "Model", size = "Model")

# Consistency
mad_026 = ggplot(performance %>% filter(tau == 0.26)) +
  aes(x = as.factor(kappa), y = mad, shape = model, size = model) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Mean absolute deviation", shape = "Model", size = "Model")

mad_063 = ggplot(performance %>% filter(tau == 0.63)) +
  aes(x = as.factor(kappa), y = mad, shape = model, size = model) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Mean absolute deviation", shape = "Model", size = "Model")

# Power
power_026 = ggplot(performance %>% filter(tau == 0.26)) +
  aes(x = as.factor(kappa), y = power, shape = model, size = model) +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Power", shape = "Model", size = "Model")

power_063 = ggplot(performance %>% filter(tau == 0.63)) +
  aes(x = as.factor(kappa), y = power, shape = model, size = model) +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = "Power", shape = "Model", size = "Model")

# Upper bound
kbar_026 = ggplot(performance %>% filter(tau == 0.26)) +
  aes(x = as.factor(kappa), y = k_hat_mean, shape = model, size = model) +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = expression(paste('Mean selected upper bound (',bar(k),')')),
       shape = "Model", size = "Model") +
  ylim(0,4)

kbar_063 = ggplot(performance %>% filter(tau == 0.63)) +
  aes(x = as.factor(kappa), y = k_hat_mean, shape = model, size = model) +
  geom_point() +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(3, 5)) +
  facet_grid(gamma ~ alpha, label = "label_parsed") +
  theme(legend.position = "bottom") +
  labs(x = expression(paste('True upper bound (',kappa,')')),
       y = expression(paste('Mean selected upper bound (',bar(k),')')),
       shape = "Model", size = "Model") +
  ylim(0,4)

# save
ggsave(plot = bias_026, filename = here("figB1.pdf"), height = 12, width = 12)
ggsave(plot = mad_026, filename = here("figB2.pdf"), height = 12, width = 12)
ggsave(plot = power_026, filename = here("figB3.pdf"), height = 12, width = 12)
ggsave(plot = kbar_026, filename = here("figB4.pdf"), height = 12, width = 12)
ggsave(plot = bias_063, filename = here("figB5.pdf"), height = 12, width = 12)
ggsave(plot = mad_063, filename = here("figB6.pdf"), height = 12, width = 12)
ggsave(plot = power_063, filename = here("figB7.pdf"), height = 12, width = 12)
ggsave(plot = kbar_063, filename = here("figB8.pdf"), height = 12, width = 12)





