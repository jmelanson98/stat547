###### Test simple model performance
# With ~biologically UNrealistic~ simulated datasets
### February 18, 2026
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
library(terra)
library(stringr)
library(purrr)


##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC 2/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
source("src/analysis_functions.R")


################################################################################
###############               SIMULATE DATA                 ###################
################################################################################

# # Load in data
sibs1 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
sibs2 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
effort1 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))

# Get colony counts (real data)
data = prep_stan_simpleforaging_bothyears(sibs1,
                                          sibs2,
                                          effort1,
                                          effort2)

# Get trap data
traps_m = data[[3]]

##### Simulate datasets #####

for (i in 1:3){
  for (rho in c(0.25,0.375, 0.5, 0.625, 0.75, 0.875,1)){
    sim = draw_simple_true_sites(sample_size = 20000,
                                 number_colonies = 1002,
                                 colony_sizes = rep(1000,1002),
                                 trap_data = traps_m,
                                 rho = rho,
                                 theta = 0.5,
                                 distance_decay = "exponentiated_quadratic")
    sim = sim %>% group_by(colonyid) %>% filter(sum(counts)>0)
    write.csv(sim, paste0("analysis/kernel_locations/modeltest/perfectsimulation", i, "_rho", rho, ".csv"))


  }
}

################################################################################
###############                  FIT MODELS                  ###################
################################################################################
# make params table
rho = c(0.25,0.375, 0.5, 0.625, 0.75, 0.875, 1)
sim = c(1,2,3)
Rmax = c(1.3, 2.7, 4.15)
stanfile = c("models/simple_multinomial_uniformdisprior.stan",
             "models/simple_multinomial_uniformdisprior_floral.stan",
             "models/simple_multinomial_nointercepts.stan",
             "models/simple_multinomial_interceptonly.stan",
             "models/simple_multinomial.stan")
params = expand.grid(rho = rho, sim = sim, Rmax = Rmax, stanfile = stanfile)

stanfilekey = data.frame(stanfile = c("models/simple_multinomial_uniformdisprior.stan",
                                      "models/simple_multinomial_uniformdisprior_floral.stan",
                                      "models/simple_multinomial_nointercepts.stan",
                                      "models/simple_multinomial_interceptonly.stan",
                                      "models/simple_multinomial.stan"),
                         key = c("multinomial_notheta", 
                                 "multinomial_theta", 
                                 "multinomial_nointercept", 
                                 "multinomial_interceptonly", 
                                 "multinomial_noRmax"))
params = left_join(params, stanfilekey)
params$Rmax[params$stanfile %in% c("models/simple_multinomial_nointercepts.stan",
                                   "models/simple_multinomial_interceptonly.stan",
                                   "models/simple_multinomial.stan"
                                   )] = NA
params$inclusion = "all"
params = unique(params)


# load in dataset
dir = paste0("analysis/kernel_locations/modeltest/perfectsimulation",params$sim[task_id], "_rho", params$rho[task_id], ".csv")
dataset = read.csv(dir)

# make modifications to stan data for each model type
stan_data = prep_stan_simpleforaging_neutraldata(dataset)

if (params$key[task_id] %in% c("multinomial_theta", "multinomial_notheta")){
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$steepness = 50
}

if (params$key[task_id] %in% c("multinomial_theta")){
  stan_data$fq = dataset$fq
}

#fit and save model
fit = stan(file = as.character(params$stanfile[task_id]),
           data = stan_data, seed = 5838299,
           chains = 4, cores = 4,
           iter = 3000, warmup = 1000,
           control = list(max_treedepth = 15),
           verbose = TRUE)

stan_data$rho = 1
stanfile = "models/fixedrho.stan"
sm = stan_model(file = stanfile)
fit_opt = rstan::optimizing(sm,
                     data = stan_data,
                      hessian = FALSE)
saveRDS(fit, 
        paste0("analysis/kernel_locations/modeltest/perfectsimulation_",
               params$key[task_id], 
               "_rho", params$rho[task_id], 
               "_Rmax", params$Rmax[task_id], 
               "_sim", params$sim[task_id],
               "_", params$inclusion[task_id], ".rds"))


# get model summary statistics
params$rho2.5 = NA
params$rho50 = NA
params$rho97.5 = NA
params$theta2.5 = NA
params$theta50 = NA
params$theta97.5 = NA
params$loglik2.5 = NA
params$loglik50 = NA
params$loglik97.5 = NA

for (i in 1:nrow(params)){
  fit = readRDS(paste0("analysis/kernel_locations/modeltest/perfectsimulation_", 
                       params$key[i], 
                       "_rho", params$rho[i], 
                       "_Rmax", params$Rmax[i], 
                       "_sim", params$sim[i],
                       "_", params$inclusion[i], ".rds"))
  summar = summary(fit)$summary
  loglik_rows = grep("^loglik\\[", rownames(summar))
  params$loglik2.5[i]= mean(summar[loglik_rows, "2.5%"])
  params$loglik50[i]= mean(summar[loglik_rows, "50%"])
  params$loglik97.5[i]= mean(summar[loglik_rows, "97.5%"])
  
  if (params$key[i] %in% c("multinomial_theta", "multinomial_notheta", "multinomial_nointercept", "multinomial_noRmax")){
    params$rho2.5[i] = summar["rho", "2.5%"]
    params$rho50[i] = summar["rho", "50%"]
    params$rho97.5[i] = summar["rho", "97.5%"]
    
    if(params$key[i] == "multinomial_theta"){
      params$theta2.5[i] = summar["theta", "2.5%"]
      params$theta50[i] = summar["theta", "50%"]
      params$theta97.5[i] = summar["theta", "97.5%"]
    }
  }
  print(paste0("Done with ", i))
}
write.csv(params, "analysis/kernel_locations/modeltest/perfectsimulationparams.csv")


# plot comparisons
result = read.csv("analysis/kernel_locations/modeltest/perfectsimulationparams.csv")
result$Rmax[is.na(result$Rmax)] = "none"
result$Rmax = as.factor(result$Rmax)


pd = position_dodge(width = 0.03)

ggplot(result[result$key != "multinomial_interceptonly",], aes(x = rho, y = rho50, color = Rmax, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  facet_grid(key~inclusion, scales = "free_y") +
  scale_y_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw()


pd = position_dodge(width = 0.3)
ggplot(result, aes(x = key, y = loglik50, group = interaction(sim, Rmax), color = Rmax)) +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = loglik2.5, ymax = loglik97.5),
                position = pd,
                width = 0) +
  facet_grid(rho~inclusion) +
  theme_bw()


