###### Test simple model performance
# With ~biologically realistic~ simulated datasets
### January 24, 2026
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
library(filelock)


##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
basepath = "analysis/kernel_locations/modeltesting/"
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
source("src/analysis_functions.R")


################################################################################
###############               SIMULATE DATA                 ###################
################################################################################

# # # Load in data
# sibs1 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
# sibs2 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
# effort1 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
# effort2 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
# samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)
# 
# # Get colony counts (real data)
# data = prep_stan_simpleforaging_bothyears(sibs1,
#                                           sibs2,
#                                           effort1,
#                                           effort2,
#                                           samplepoints)
# 
# # Get trap data
# traps_m = data[[3]]
# 
# ##### Simulate datasets #####
# 
# for (i in 1:3){
#   for (rho in c(0.875,1)){
#     sim = draw_simple_true_sites(sample_size = 2000,
#                                  number_colonies = 11700,
#                                  colony_sizes = rep(20,11700),
#                                  trap_data = traps_m,
#                                  rho = rho,
#                                  theta = 0.5,
#                                  distance_decay = "exponentiated_quadratic")
#     sim = sim %>% group_by(colonyid) %>% filter(sum(counts)>0)
#     write.csv(sim, paste0(basepath, "modeltest/simulation", i, "_rho", rho, ".csv"))
# 
# 
#   }
# }

################################################################################
###############                  FIT MODELS                  ###################
################################################################################
# make params table
rho = c(0.25,0.375, 0.5, 0.625, 0.75, 0.875, 1)
sim = c(1,2,3)
Rmax = c(1.3, 2.7, 4.15)
deltaprior = c(0.5, 1, 2)
stanfile = c("models/simple_multinomial_uniformdisprior.stan",
          "models/simple_multinomial_normaldisprior.stan",
          "models/deltaprior.stan",
          "models/simple_multinomial.stan")
params = expand.grid(rho = rho, sim = sim, Rmax = Rmax, stanfile = stanfile, deltaprior = deltaprior)

stanfilekey = data.frame(stanfile = c("models/simple_multinomial_uniformdisprior.stan",
                                      "models/simple_multinomial_normaldisprior.stan",
                                      "models/deltaprior.stan",
                                      "models/simple_multinomial.stan"),
                         key = c("uniformdisprior", "normaldisprior", "deltaprior", "none"))
params = left_join(params, stanfilekey)
params$Rmax[params$stanfile %in% c("models/simple_multinomial.stan",
                              "models/deltaprior.stan")] = NA
params$deltaprior[params$stanfile %in% c("models/simple_multinomial.stan",
                                         "models/simple_multinomial_uniformdisprior.stan",
                                         "models/simple_multinomial_normaldisprior.stan")] = NA
params = unique(params)


# load in dataset
dir = paste0(basepath, "modeltest/simulation",params$sim[task_id], "_rho", params$rho[task_id], ".csv")
dataset = read.csv(dir)


# make modifications to stan data for each model type
stan_data = prep_stan_simpleforaging_neutraldata(dataset)

if (params$key[task_id] == "uniformdisprior"){
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$steepness = 50
  stan_data$colonycenters = NULL
} else if (params$key[task_id] == "normaldisprior") {
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$colonycenters = NULL
} else if (params$key[task_id] == "deltaprior"){
    stan_data$colonybounds = NULL
    stan_data$deltaprior = params$deltaprior[task_id]
} else if (params$key[task_id] == "none") {
    stan_data$colonycenters = NULL
  }



#fit and save model
# fit = stan(file = as.character(params$stanfile[task_id]),
#                   data = stan_data, seed = 5838299,
#                   chains = 4, cores = 4,
#                   iter = 3000, warmup = 1000,
#                   control = list(max_treedepth = 15),
#                   verbose = TRUE)
# saveRDS(fit,
#         paste0(basepath, "modeltest/",
#                params$key[task_id],
#                "_rho", params$rho[task_id],
#                "_Rmax", params$Rmax[task_id],
#                "_deltaprior", params$deltaprior[task_id],
#                "_sim", params$sim[task_id],
#                "_", params$inclusion[task_id], ".rds"))


# parallel way
myrow = params[task_id,]
fit = readRDS(paste0("~/projects/def-ckremen/melanson/fv_landscapeforaging/analysis/kernel_locations/modeltesting/modeltest/",
                        params$key[task_id],
                        "_rho", params$rho[task_id],
                        "_Rmax", params$Rmax[task_id],
                        "_deltaprior", params$deltaprior[task_id],
                       "_sim", params$sim[task_id],
                        "_.rds"))
summar = summary(fit)$summary
loglik_rows = grep("^loglik\\[", rownames(summar))
myrow$loglik2.5= mean(summar[loglik_rows, "2.5%"])
myrow$loglik50= mean(summar[loglik_rows, "50%"])
myrow$loglik97.5= mean(summar[loglik_rows, "97.5%"])

myrow$rho2.5 = summar["rho", "2.5%"]
myrow$rho50 = summar["rho", "50%"]
myrow$rho97.5 = summar["rho", "97.5%"]
print(myrow)

if(file.exists(paste0(basepath, "saved.csv"))){

  lock_file = paste0(basepath,"saved.lock")
  csv_file = paste0(basepath,"saved.csv")
  lock = lock(lock_file, timeout = 60000)  
  
  if (!is.null(lock)) {
    data = read.csv(csv_file)
    
    data = rbind(data, myrow)
    
    write.csv(data, csv_file, row.names = FALSE)
    
    # release lock
    unlock(lock)
}
print('hello!')
} else {
print('hello, creating file!')
  write.csv(myrow, paste0(basepath, "saved.csv"), row.names = FALSE)
}



#plot comparisons
result = read.csv(paste0(basepath, "saved.csv"))
result$Rmax[is.na(result$Rmax)] = "none"
result$Rmax = as.factor(result$Rmax)
result$deltaprior[is.na(result$deltaprior)] = "none"
result$deltaprior = as.factor(result$deltaprior)


pd = position_dodge(width = 0.03)
# Set color scheme for plots
faded_pale = "#D2E4D4"
faded_light = "#B6D3B8"
faded_medium = "#609F65"
faded_strong = "#4E8353"
faded_green = "#355938"
faded_dark = "#1B2D1C"


light_gold = "#F7EAC0"
lm_gold = "#F2DC97"
medium_gold = "#ECCF6F"
gold = "#E2B41D"
dark_gold = "#B99318"
darker_gold = "#907313"

uniform = ggplot(result[result$key == "uniformdisprior",], aes(x = rho, y = rho50, color = Rmax, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  xlab(expression(rho[G])) +
  ylab("") +
  scale_color_manual(
    name = expression(R[max]),
    values = c(faded_green, faded_medium, faded_light)) +
  scale_y_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw()

normal = ggplot(result[result$key == "normaldisprior",], aes(x = rho, y = rho50, color = Rmax, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  xlab(expression(rho[G])) +
  ylab(expression(rho[I])) +
  scale_color_manual(
    name = expression(R[max]),
    values = c(faded_green, faded_medium, faded_light)) +
  scale_y_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(legend.position = "none")

deltaprior = ggplot(result[(result$deltaprior == 1 | result$deltaprior == 2),], aes(x = rho, y = rho50, color = deltaprior, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  scale_color_manual(
    name = expression(sigma[d]),
    values = c(dark_gold, medium_gold)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw()

noprior = ggplot(result[result$key == "none",], aes(x = rho, y = rho50, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  xlab("") +
  ylab(expression(rho[I])) +
  scale_color_discrete(name = expression(sigma[d])) +
  scale_y_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0, 0.05))) +
  theme_bw()

grid.arrange(noprior, deltaprior,
             normal, uniform,
             ncol =2, widths = c(6.5,8),
             heights = c(8,9))



library(cowplot)

aligned <- align_plots(
  noprior, deltaprior,
  normal, uniform,
  align = "hv",      # horizontal + vertical
  axis = "tblr"      # align all axes
)
out = plot_grid(
  aligned[[1]], aligned[[2]],
  aligned[[3]], aligned[[4]],
  ncol = 2,
  labels = c("(A)", "(B)", "(C)", "(D)"),
  label_size = 14
)

ggsave("analysis/kernel_locations/modeltesting/report_figures/simplesim.png",
       out, height = 2000, width = 2500, units = "px")


