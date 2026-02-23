###### Fit foraging distance model to simulated data!
### With floral x landscape effects and ~multiple timepoints~
### February 10, 2026
### J. Melanson


##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(sf)
library(stringr)
library(boot)
library(filelock)

# Prep workspace
# local
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC 2/bombus_project"
#bombus_path = "/Users/jenna1/Documents/bombus_project"
#bombus_path = "/Users/jenna1/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/UBC/bombus_project"

# remote
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "~/projects/def-ckremen/melanson/"
source("src/ff_helper.R")

# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
veg2022 = read.csv(paste0(bombus_path, "/raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "/raw_data/2023vegetationdata.csv"))

### Make some plots of how we expect rho to change with rho1, rho2, rho3, rhomax
# rho0 = 0
# rho1 = 0.5
# rho2 = -0.5
# rho3 = -0.1
# predictors = expand.grid(fq = seq(-3,3, by = 0.01),
#                         lq = seq(-3,3, by = 0.01))
# 
# df_long = crossing(predictors) %>%
#   mutate(
#     rho0 = rho0,
#     rho1 = rho1,
#     rho2 = rho2,
#     rho3 = rho3,
#     rho = inv.logit(rho0+rho1*fq+rho2*lq+rho3*fq*lq))
# 
# ggplot(df_long, aes(x = fq, y = rho, color = lq)) +
#   geom_point(size = 0.01) +
#   theme_bw()

### Make params grid
rho1 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
rho2 = c(-0.5, -0.25, 0, 0.25, 0.5)
rho3 = c(-0.1, -0.05, 0, 0.05, 0.1)
params = expand.grid(rho1 = rho1, rho2 = rho2, rho3 = rho3)


### Simulate data
sim = draw_floral_true_sites(ncolonies = 6000, # positive integer, per year
                             colony_sizes= rep(20,6000), #vector of length ncolonies
                             rho0 = 0,
                             rho1 = params$rho1[task_id], # floral attractiveness
                             rho2 = params$rho2[task_id], #landscape effect
                             rho3 = params$rho3[task_id], # interactions of local & landscape scale
                             rhomax = 1,
                             effort1 = effort2022,
                             effort2 = effort2023,
                             spec1 = mixtus_sibs2022,
                             spec2 = mixtus_sibs2023
)
sim = sim %>%
  group_by(colonyid) %>%
  filter(sum(counts) > 0)
write.csv(sim, paste0("analysis/floral_foraging/floraltest/simulation_rho00_rho1", 
                      params$rho1[task_id], "_rho2", 
                      params$rho2[task_id], "_rho3", 
                      params$rho3[task_id], "_rhomax1.csv"))

### Get stan data for mixtus sim
sim = read.csv(paste0("analysis/floral_foraging/floraltest/simulation_rho00_rho1", params$rho1[task_id],
                      "_rho2", params$rho2[task_id], 
                      "_rho3", params$rho3[task_id], 
                      "_rhomax1.csv"))
mixtus_sim_data = prep_stan_floralforaging_simdata(sim, rhomax = 1)
mixtus_sim_data$deltaprior = 1
mixtus_sim_data$colonybounds = NULL
mixtus_sim_data$rhomax = NULL

### Fit floral foraging model
stanfile = "models/floral_foraging.stan"

fit = stan(file = stanfile,
          data = mixtus_sim_data, seed = 5838299,
          chains = 4, cores = 4,
          iter = 2000,
          refresh = 10,
          verbose = TRUE)
saveRDS(fit, paste0("analysis/floral_foraging/floraltest/simulation_rho00_rho1", params$rho1[task_id],
                    "_rho2", params$rho2[task_id], 
                    "_rho3", params$rho3[task_id], 
                    "_rhomax1.rds"))

# parallel way
myrow = params[task_id,]
summar = summary(fit)$summary
loglik_rows = grep("^loglik\\[", rownames(summar))
myrow$loglik2.5= mean(summar[loglik_rows, "2.5%"])
myrow$loglik50= mean(summar[loglik_rows, "50%"])
myrow$loglik97.5= mean(summar[loglik_rows, "97.5%"])

myrow$rho0_2.5 = summar["rho0", "2.5%"]
myrow$rho0_50 = summar["rho0", "50%"]
myrow$rho0_97.5 = summar["rho0", "97.5%"]

myrow$rho1_2.5 = summar["rho1", "2.5%"]
myrow$rho1_50 = summar["rho1", "50%"]
myrow$rho1_97.5 = summar["rho1", "97.5%"]

myrow$rho2_2.5 = summar["rho2", "2.5%"]
myrow$rho2_50 = summar["rho2", "50%"]
myrow$rho2_97.5 = summar["rho2", "97.5%"]

myrow$rho3_2.5 = summar["rho3", "2.5%"]
myrow$rho3_50 = summar["rho3", "50%"]
myrow$rho3_97.5 = summar["rho3", "97.5%"]

myrow$rhomax_2.5 = summar["rhomax", "2.5%"]
myrow$rhomax_50 = summar["rhomax", "50%"]
myrow$rhomax_97.5 = summar["rhomax", "97.5%"]

if(file.exists("analysis/floral_foraging/saved.csv")){
  
  lock_file = "analysis/floral_foraging/saved.lock"
  csv_file = "analysis/floral_foraging/saved.csv"
  lock = lock(lock_file, timeout = 60000)  
  
  if (!is.null(lock)) {
    data = read.csv(csv_file)
    
    data = rbind(data, myrow)
    
    write.csv(data, csv_file, row.names = FALSE)
    
    # release lock
    unlock(lock)
  }
} else{
  write.csv(myrow, "analysis/floral_foraging/saved.csv", row.names = FALSE)
}


# #plot comparisons
# result = read.csv(paste0(basepath, "modeltest/params_new.csv"))
# result$Rmax[is.na(result$Rmax)] = "none"
# result$Rmax = as.factor(result$Rmax)
# result$deltaprior[is.na(result$deltaprior)] = "none"
# result$deltaprior = as.factor(result$deltaprior)
# 
# 
# pd = position_dodge(width = 0.03)
# 
# ggplot(result[result$key != "multinomial_interceptonly",], aes(x = rho, y = rho50, color = Rmax, group = sim)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
#   geom_point(position = pd, size = 1) +
#   geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
#                 position = pd,
#                 width = 0) +
#   facet_grid(key~inclusion, scales = "free_y") +
#   scale_y_continuous(
#     limits = c(0.2, NA),
#     expand = expansion(mult = c(0, 0.05))) +
#   scale_x_continuous(
#     limits = c(0.2, NA),
#     expand = expansion(mult = c(0, 0.05))) +
#   theme_bw()
# 
# 
# pd = position_dodge(width = 0.3)
# ggplot(result, aes(x = key, y = loglik50, group = interaction(sim, Rmax), color = Rmax)) +
#   geom_point(position = pd, size = 1) +
#   geom_errorbar(aes(ymin = loglik2.5, ymax = loglik97.5),
#                 position = pd,
#                 width = 0) +
#   facet_grid(rho~inclusion) +
#   theme_bw()
