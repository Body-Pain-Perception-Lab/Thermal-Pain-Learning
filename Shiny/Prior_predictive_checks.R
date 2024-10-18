# These scripts are simulations of nsim trajectories for sampled parameter values of each learning model. 
# The script therefore produces prior predictive checks i.e. nsim of them with the parameter values settled on for the paper:

if (!require("pacman")) install.packages("pacman")
pacman::p_load("renv", "here", "tidyverse","shiny","brms")

# Load your utility functions
source(here::here("Shiny","utility_shiny.R"))

# get contngency space from one participant (two options of contingency1.csv or contingency2.csv)
u = read.csv(here::here("Shiny","Data","contingency1.csv"))
u = u$x


# Number of simulations:
nsim = 50

## parameters
# for all models:
e0 =  brms::inv_logit_scaled(rnorm(nsim,0,1))
zeta = exp(rnorm(nsim,1,1))

# HGF
pi2_0 = exp(rnorm(nsim,1,1))

omega = exp(rnorm(nsim,-4,3))

# RW
alpha = brms::inv_logit_scaled(rnorm(nsim,-1,2))

# Pearce hall:
S = brms::inv_logit_scaled(rnorm(nsim,0,2))

a_0 = brms::inv_logit_scaled(rnorm(nsim,0,1))

# Sutton k1
mu = exp((rnorm(nsim,log(3),3)))

Rhat = exp(rnorm(nsim,10,1))

h_0 = exp(rnorm(nsim,-1,1))

# plots
hgf = data.frame(e0 = e0,zeta = zeta, omega = omega,pi2_0=pi2_0) %>% mutate(draw = 1:n()) %>% rowwise() %>% 
  mutate(hgf = list(get_hgf(u,e0,zeta,omega,pi2_0))) %>% select(draw,hgf) %>% unnest() %>% 
  ggplot(aes(x = trial, y = expect, group = draw))+geom_line(alpha = 0.25)+
  theme_minimal()+ggtitle("Hierachical Gaussian Filter")



rw = data.frame(e0 = e0,zeta = zeta, alpha = alpha) %>% mutate(draw = 1:n()) %>% rowwise() %>% 
  mutate(rw = list(get_rw(u,e0,zeta,alpha))) %>% select(draw,rw) %>% unnest() %>% 
  ggplot(aes(x = trial, y = expect, group = draw))+geom_line(alpha = 0.25)+
  theme_minimal()+ggtitle("Rescorla wagner")



ph = data.frame(e0 = e0,zeta = zeta, S = S, a_0 = a_0) %>% mutate(draw = 1:n()) %>% rowwise() %>% 
  mutate(ph = list(get_ph(u,e0,zeta,S = S, a_0 = a_0))) %>% select(draw,ph) %>% unnest() %>% 
  ggplot(aes(x = trial, y = expect, group = draw))+geom_line(alpha = 0.25)+
  theme_minimal()+ggtitle("Pearce hall")



su1 = data.frame(e0 = e0,zeta = zeta,mu = mu,Rhat = Rhat, h_0 = h_0) %>% mutate(draw = 1:n()) %>% rowwise() %>% 
  mutate(su1 = list(get_su1(u,e0,zeta,mu = mu,Rhat = Rhat, h_0 = h_0))) %>% select(draw,su1) %>% unnest() %>% 
  ggplot(aes(x = trial, y = expect, group = draw))+geom_line(alpha = 0.25)+
  theme_minimal()+ggtitle("Sutton k1")


library(patchwork)
hgf+rw+ph+su1
