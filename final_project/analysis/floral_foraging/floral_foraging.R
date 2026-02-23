###### Fit foraging distance model real data!
### With floral effects and ~multiple timepoints~
### December 5, 2025
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

# Prep workspace
# local
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC 2/bombus_project"

# remote
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "~/projects/def-ckremen/melanson/"
source("src/ff_helper.R")

# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
veg2022 = read.csv(paste0(bombus_path, "/raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "/raw_data/2023vegetationdata.csv"))

### Get stan data for mixtus and impatiens
mixtus_data = prep_stan_floralforaging(mixtus_sibs2022,
                                       mixtus_sibs2023,
                                       effort2022,
                                       effort2023,
                                       veg2022,
                                       veg2023)

impatiens_data = prep_stan_floralforaging(impatiens_sibs2022,
                                       impatiens_sibs2023,
                                       effort2022,
                                       effort2023,
                                       veg2022,
                                       veg2023)


# Get task id
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

stanfile = "models/floral_foraging_zinb.stan"

if (task_id == 1){
  stanFit = stan(file = stanfile,
                 data = mixtus_data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 2000,
                 refresh = 10,
                 verbose = TRUE)
  saveRDS(stanFit, "analysis/foraging_modelfits/mixtusfloralfit.rds")
} else{
  stanFit = stan(file = stanfile,
                 data = impatiens_data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 2000,
                 refresh = 10,
                 verbose = TRUE)
  saveRDS(stanFit, "analysis/foraging_modelfits/impatiensfloralfit.rds")
}



# #Plot the posteriors of some colonies
# 
# plot_list = list()
# numplots = 9
# legends = list()
# 
# for (i in 1:numplots){
#   delta_draws = cbind(x = rstan::extract(stanFit, pars = "delta_x")$delta[, i],
#                       y = rstan::extract(stanFit, pars = "delta_y")$delta[, i])
# 
#   p = ggplot(delta_draws, aes(x = x, y = y)) +
#     geom_density_2d_filled(alpha = 0.8) +
# 
#     #plot trap locations / sizes / quality
#     geom_point(data = CKT[CKT$sibshipID ==i,], aes(x = trap_x, y = trap_y, size = counts, colour = "red")) +
#     scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
# 
#     #miscellaneous
#     labs(title = paste("Colony", i),
#          size = "Number of Captures",
#          level = "Colony Posterior") +
#     guides(colour = "none") +
#     #xlim(c(1220, 1225)) +
#     #ylim(c(455,460)) +
#     coord_equal() +
#     theme_bw()
# 
#   # save legend
#   g <- ggplotGrob(p)
#   legend_index <- which(g$layout$name == "guide-box-right")
#   legend <- g$grobs[[legend_index]]
# 
#   # remove legend from plot
#   p <- p + theme(legend.position = "none")
# 
#   #save plot
#   plot_list[[i]] = p
#   legends[[1]] = legend
# }
# 
# fig = grid.arrange(grobs = plot_list, ncol = 3)
# fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
# ggsave(paste0("analysis/colony_posteriors/foragingmodel_", task_id, ".jpg"), fig, height = 3000, width = 4000, units = "px")