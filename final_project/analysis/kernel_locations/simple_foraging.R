###### Fit foraging distance model to real data!
### Get estimates of cumulative colony posteriors
### Get individual colony foraging kernels
### November 14, 2025
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
library(grid)
library(gridExtra)
library(tibble)
library(sf)
library(terra)
library(stringr)

##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
source("src/analysis_functions.R")


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


# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))

# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Create params table
Rmax = c(1.3, 2.7, 4.15)
deltaprior = NA
species = c("mixtus", "impatiens")
model = c("normal", "uniform")
params1 = expand.grid(Rmax = Rmax, deltaprior = deltaprior, species = species, model = model)

Rmax = NA
deltaprior = c(0.5,1,2)
species = c("mixtus", "impatiens")
model = "deltaprior"
params2 = expand.grid(Rmax = Rmax, deltaprior = deltaprior, species = species, model = model)

params = rbind(params1, params2)
params$task_id = 1:nrow(params)

#################################################################
# Prepare data for Stan
#################################################################

if (params$species[task_id] == "mixtus"){
  data = prep_stan_simpleforaging_bothyears(mixtus_sibs2022,
                                            mixtus_sibs2023,
                                            effort2022,
                                            effort2023,
                                            samplepoints)
} else if (params$species[task_id] == "impatiens"){
  data = prep_stan_simpleforaging_bothyears(impatiens_sibs2022,
                                            impatiens_sibs2023,
                                            effort2022,
                                            effort2023,
                                            samplepoints)
}

stan_data = data[[1]]
CKT = data[[2]]


#################################################################
# Think about different likelihood distributions for distance
#################################################################
# probability density for log1p_exp
# x = seq(from = 0, to = 3, by = 0.01)
# Rmax = 2
# steepness = 50
# flog1p = function(x) {
#   return(exp(-log(1+exp((x-Rmax)*steepness))))
# }
# likelihood = flog1p(x)
# defintegral = integrate(flog1p, lower = 0, upper = Inf)
# prob_density = likelihood/Rmax
# plot(x, likelihood)
# plot(x,prob_density)
# 
# # likelihood for exponential
# loglik2 = -exp(steepness*(x-Rmax))
# lik2 = exp(loglik2)
# plot(x, lik2)
# 
# # likelihood for half normal
# library(fdrtool)
# prob_density = dhalfnorm(x, theta=sqrt(pi/2)/(Rmax/3), log = FALSE)
# plot(x, prob_density)
# hn = function(x) {
#   return(dhalfnorm(x, theta=sqrt(pi/2)/(Rmax/3), log = FALSE))
# }
# defintegral = integrate(hn, lower = 0, upper = Inf)

#################################################################
# Fit models
#################################################################

if (params$model[task_id] == "normal"){
  stanfile = "models/simple_multinomial_normaldisprior.stan"
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$colonycenters = NULL
  fit = stan(file = stanfile,
             data = stan_data, seed = 5838299,
             chains = 4, cores = 4,
             iter = 4000, warmup = 1000,
             verbose = TRUE)
  modelname = paste0(params$species[task_id], "_normal_Rmax", params$Rmax[task_id])
  saveRDS(fit, paste0("analysis/kernel_locations/foraging_modelfits/", modelname, ".rds"))
} else if (params$model[task_id] == "uniform") {
  stanfile = "models/simple_multinomial_uniformdisprior.stan"
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$steepness = 50
  stan_data$colonycenters = NULL
  fit = stan(file = stanfile,
             data = stan_data, seed = 5838299,
             chains = 4, cores = 4,
             iter = 4000, warmup = 1000,
             verbose = TRUE)
  modelname = paste0(params$species[task_id], "_uniform_Rmax", params$Rmax[task_id])
  saveRDS(fit, paste0("analysis/kernel_locations/foraging_modelfits/", modelname,".rds"))
} else if (params$model[task_id] == "deltaprior"){
  stanfile = "models/deltaprior.stan"
  stan_data$colonybounds = NULL
  stan_data$deltaprior = params$deltaprior[task_id]
  fit = stan(file = stanfile,
             data = stan_data, seed = 5838299,
             chains = 4, cores = 4,
             iter = 4000, warmup = 1000,
             verbose = TRUE)
  modelname = paste0(params$species[task_id], "_deltaprior", params$deltaprior[task_id])
  saveRDS(fit, paste0("analysis/kernel_locations/foraging_modelfits/", modelname,".rds"))
}


#############################################
# Extract values and makes some plots
#############################################
# Load in uniform models
mix.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_uniform_Rmax1.3.rds")
imp.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_uniform_Rmax1.3.rds")
mix.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_uniform_Rmax2.7.rds")
imp.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_uniform_Rmax2.7.rds")
mix.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_uniform_Rmax4.15.rds")
imp.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_uniform_Rmax4.15.rds")

# Load in normal models
normmix.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax1.3.rds")
normimp.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_normal_Rmax1.3.rds")
normmix.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax2.7.rds")
normimp.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_normal_Rmax2.7.rds")
normmix.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax4.15.rds")
normimp.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_normal_Rmax4.15.rds")

# Load in deltaprior models
deltamix.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_deltaprior0.5.rds")
deltaimp.fit.min = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_deltaprior0.5.rds")
deltamix.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_deltaprior1.rds")
deltaimp.fit.med = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_deltaprior1.rds")
deltamix.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/mixtus_deltaprior2.rds")
deltaimp.fit.max = readRDS("analysis/kernel_locations/foraging_modelfits/impatiens_deltaprior2.rds")


postm_min = as.data.frame(mix.fit.min)
postm_med = as.data.frame(mix.fit.med)
postm_max = as.data.frame(mix.fit.max)
posti_min = as.data.frame(imp.fit.min)
posti_med = as.data.frame(imp.fit.med)
posti_max = as.data.frame(imp.fit.max)
normpostm_min = as.data.frame(normmix.fit.min)
normpostm_med = as.data.frame(normmix.fit.med)
normpostm_max = as.data.frame(normmix.fit.max)
normposti_min = as.data.frame(normimp.fit.min)
normposti_med = as.data.frame(normimp.fit.med)
normposti_max = as.data.frame(normimp.fit.max)
deltapostm_min = as.data.frame(deltamix.fit.min)
deltapostm_med = as.data.frame(deltamix.fit.med)
deltapostm_max = as.data.frame(deltamix.fit.max)
deltaposti_min = as.data.frame(deltaimp.fit.min)
deltaposti_med = as.data.frame(deltaimp.fit.med)
deltaposti_max = as.data.frame(deltaimp.fit.max)

combinedpost =  bind_rows(
  data.frame(model = "B. mixtus", rho = postm_min$rho, Rmax = 1.3, distribution = "uniform"),
  data.frame(model = "B. mixtus", rho = postm_med$rho, Rmax = 2.7, distribution = "uniform"),
  data.frame(model = "B. mixtus", rho = postm_max$rho, Rmax = 4.15, distribution = "uniform"),
  data.frame(model = "B. impatiens", rho = posti_min$rho, Rmax = 1.3, distribution = "uniform"),
  data.frame(model = "B. impatiens", rho = posti_med$rho, Rmax = 2.7, distribution = "uniform"),
  data.frame(model = "B. impatiens", rho = posti_max$rho, Rmax = 4.15, distribution = "uniform"),
  data.frame(model = "B. mixtus", rho = normpostm_min$rho, Rmax = 1.3, distribution = "normal"),
  data.frame(model = "B. mixtus", rho = normpostm_med$rho, Rmax = 2.7, distribution = "normal"),
  data.frame(model = "B. mixtus", rho = normpostm_max$rho, Rmax = 4.15, distribution = "normal"),
  data.frame(model = "B. impatiens", rho = normposti_min$rho, Rmax = 1.3, distribution = "normal"),
  data.frame(model = "B. impatiens", rho = normposti_med$rho, Rmax = 2.7, distribution = "normal"),
  data.frame(model = "B. impatiens", rho = normposti_max$rho, Rmax = 4.15, distribution = "normal"),
  data.frame(model = "B. mixtus", rho = deltapostm_min$rho, Rmax = 1.3, distribution = "deltaprior"),
  data.frame(model = "B. mixtus", rho = deltapostm_med$rho, Rmax = 2.7, distribution = "deltaprior"),
  data.frame(model = "B. mixtus", rho = deltapostm_max$rho, Rmax = 4.15, distribution = "deltaprior"),
  data.frame(model = "B. impatiens", rho = deltaposti_min$rho, Rmax = 1.3, distribution = "deltaprior"),
  data.frame(model = "B. impatiens", rho = deltaposti_med$rho, Rmax = 2.7, distribution = "deltaprior"),
  data.frame(model = "B. impatiens", rho = deltaposti_max$rho, Rmax = 4.15, distribution = "deltaprior")
)

write.csv(combinedpost, "analysis/kernel_locations/real_rhoposteriors.csv")
mix = read.csv("analysis/kernel_locations/mixtusrhoposteriors.csv")
imp1 = read.csv("analysis/kernel_locations/impatiensrhoposteriors1.csv")
imp2 = read.csv("analysis/kernel_locations/impatiensrhoposteriors2wq.csv")
combinedpost = rbind(mix, imp1, imp2)
combinedpost$Rmax[combinedpost$Rmax == 1.3] = 0.5
combinedpost$Rmax[combinedpost$Rmax == 2.7] = 1
combinedpost$Rmax[combinedpost$Rmax == 4.15] = 2
combinedwide <- combinedpost %>%
  group_by(Rmax, model) %>%
  mutate(draw_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(Rmax, draw_id),
    names_from = model,
    values_from = rho
  )
combinedwide$diff = combinedwide$`B. impatiens` - combinedwide$`B. mixtus`

posteriors = ggplot(combinedpost, aes(x = rho, fill = model, color = model)) +
  scale_fill_manual(values = c(lm_gold, faded_strong)) +
  scale_color_manual(values = c(dark_gold, faded_dark)) +
  geom_histogram(aes(y = after_stat(density))) +
  facet_grid(~Rmax) +
  ylab("posterior density") +
  xlab(expression(rho)) +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

posteriordiff = ggplot(combinedwide, aes(x = diff)) +
  geom_histogram(aes(y = after_stat(density)), fill = "lavender", color = "darkorchid4") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(~Rmax) +
  ylab("posterior density") +
  xlab(expression(rho[imp] - rho[mix])) +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

grid.arrange(posteriors,posteriordiff, ncol = 1)
ggsave("analysis/kernel_locations/modeltesting/report_figures/impvsmix.png", posteriordiff, height = 1000, width = 2000, units = "px")
ggsave("docs/appendix_figures/foraging_sensitivity.jpg", posteriors, height = 1500, width = 1600, units = "px")

# Now another for the main body -- just 2.7 km
unif27 = ggplot(combinedpost[combinedpost$Rmax == 2.7 & combinedpost$distribution == "uniform",], aes(x = rho, fill = model, color = model)) +
  scale_fill_manual(values = c(lm_gold, faded_strong)) +
  scale_color_manual(values = c(dark_gold, faded_dark)) +
  geom_histogram(aes(y = after_stat(density))) +
  ylab("Posterior density") +
  xlab(expression(rho)) +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
unif27

# Make plots of foraging decay with distance
# Incorporate posterior uncertainty, for each Rmax!
x = seq(0, 2, by = 0.0001)
mix_df = cbind(data.frame(x = x),
               as.data.frame(matrix(NA, nrow = length(x), ncol = 1000)))
imp_df = cbind(data.frame(x = x),
               as.data.frame(matrix(NA, nrow = length(x), ncol = 1000)))
for (draw in 1:1000){
  mixrho = sample(combinedpost$rho[combinedpost$model == "B. mixtus" & combinedpost$Rmax == 2.7 & combinedpost$distribution == "uniform"], 1)
  mix_df[,draw+1] = exp(-0.5*(x/mixrho)^2)
  imprho = sample(combinedpost$rho[combinedpost$model == "B. impatiens" & combinedpost$Rmax == 2.7 & combinedpost$distribution == "uniform"], 1)
  imp_df[,draw+1] = exp(-0.5*(x/imprho)^2)
}

# get summary stats
y_mix = as.matrix(mix_df[ , -1])
mix_df$mean  <- rowMeans(y_mix, na.rm = TRUE)
mix_df$lower <- apply(y_mix, 1, quantile, probs = 0.025, na.rm = TRUE)
mix_df$upper <- apply(y_mix, 1, quantile, probs = 0.975, na.rm = TRUE)

y_imp = as.matrix(imp_df[ , -1])
imp_df$mean  <- rowMeans(y_imp, na.rm = TRUE)
imp_df$lower <- apply(y_imp, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_df$upper <- apply(y_imp, 1, quantile, probs = 0.975, na.rm = TRUE)

foraging_decay = ggplot() +
  geom_line(data = mix_df, aes(x = x, y = mean), color = faded_dark, size = 1) +
  geom_ribbon(data = mix_df, aes(x = x, ymin = lower, ymax = upper), fill = faded_medium, alpha = 0.2) +
  geom_vline(xintercept = 1.96*mean(
    combinedpost$rho[combinedpost$model == "B. mixtus" & combinedpost$Rmax == 2.7 & combinedpost$distribution == "uniform"]),
    color = faded_dark, linetype = "dashed") +
  geom_line(data = imp_df, aes(x = x, y = mean), color = darker_gold, size = 1) +
  geom_ribbon(data = imp_df, aes(x = x, ymin = lower, ymax = upper), fill = medium_gold, alpha = 0.2) +
  geom_vline(xintercept = 1.96*mean(
    combinedpost$rho[combinedpost$model == "B. impatiens" & combinedpost$Rmax == 2.7 & combinedpost$distribution == "uniform"]),
    color = darker_gold, linetype = "dashed") +
  ylab("Relative visitation") +
  xlab("Distance from nest (km)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(),
    panel.grid.minor = element_blank()
  )
foraging_decay

# make grid
grid = grid.arrange(unif27, foraging_decay, ncol = 2)
ggsave("figures/manuscript_figures/foragingmain.png", grid, width = 2000, height = 1000, units = "px")


# #Plot the posteriors of some colonies
# #colonies = c(5, 7, 11, 18)
# 
# delta_draws1 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1024],
#                       y = rstan::extract(fit, pars = "delta_y")$delta[, 1024]))
# delta_draws1$colony = "Colony 1024"
# delta_draws2 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1373],
#                      y = rstan::extract(fit, pars = "delta_y")$delta[, 1373]))
# delta_draws2$colony = "Colony 1373"
# delta_draws3 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1399],
#                      y = rstan::extract(fit, pars = "delta_y")$delta[, 1399]))
# delta_draws3$colony = "Colony 1399"
# 
# plot1 = ggplot(delta_draws1, aes(x = x, y = y)) +
#     geom_density_2d_filled(alpha = 0.8, bins = 9) +
#     scale_fill_brewer(palette = "Greens") +
# 
#     #plot trap locations / sizes / quality
#     geom_point(data = CKT[CKT$stansibkey ==1024,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
#     scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
#     xlim(c(min(CKT[CKT$stansibkey ==1024,]$trap_x)-1, max(CKT[CKT$stansibkey ==1024,]$trap_x) + 1)) +
#     ylim(c(min(CKT[CKT$stansibkey ==1024,]$trap_y)-1, max(CKT[CKT$stansibkey ==1024,]$trap_y) + 1)) +
# 
#     #miscellaneous
#     labs(x = "",
#          y = "Latitude",
#          title = "Colony 1024") +
#     guides(fill = "none", colour = "none", size = "none") +
#     theme(legend.position = "none") +
#     coord_equal() +
#     theme_bw()
# 
# plot2 = ggplot(delta_draws2, aes(x = x, y = y)) +
#   geom_density_2d_filled(alpha = 0.8, bins = 9) +
#   scale_fill_brewer(palette = "Greens") +
# 
#   #plot trap locations / sizes / quality
#   geom_point(data = CKT[CKT$stansibkey ==1373,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
#   scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
#   xlim(c(min(CKT[CKT$stansibkey ==1373,]$trap_x)-1, max(CKT[CKT$stansibkey ==1373,]$trap_x) + 1)) +
#   ylim(c(min(CKT[CKT$stansibkey ==1373,]$trap_y)-1, max(CKT[CKT$stansibkey ==1373,]$trap_y) + 1)) +
# 
#   #miscellaneous
#   labs(x = "Longitude",
#        y = "",
#        title = "Colony 1373") +
#   guides(fill = "none", colour = "none", size = "none") +
#   theme(legend.position = "none") +
#   coord_equal() +
#   theme_bw()
# 
# plot3 = ggplot(delta_draws3, aes(x = x, y = y)) +
#   geom_density_2d_filled(alpha = 0.8, bins = 9) +
#   scale_fill_brewer(palette = "Greens") +
# 
#   #plot trap locations / sizes / quality
#   geom_point(data = CKT[CKT$stansibkey ==1399,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
#   scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
#   xlim(c(min(CKT[CKT$stansibkey ==1399,]$trap_x)-1, max(CKT[CKT$stansibkey ==1399,]$trap_x) + 1)) +
#   ylim(c(min(CKT[CKT$stansibkey ==1399,]$trap_y)-1, max(CKT[CKT$stansibkey ==1399,]$trap_y) + 1)) +
# 
#   #miscellaneous
#   labs(x = "",
#        y = "",
#        title = "Colony 1399") +
#   guides(fill = "none", colour = "none", size = "none") +
#   theme(legend.position = "none") +
#   coord_equal() +
#   theme_bw()
# 
# #### Make a grid plot
# grid = grid.arrange(plot1, nullGrob(), plot2, nullGrob(), plot3,
#                     ncol = 5, widths = c(10, 1, 10,1,10))
# 
# ggsave(paste0("figures/colonyposteriors", modelname, ".jpg"), grid, height = 3000, width = 3200, units = "px")
# 


#################################################################
# Compute raster layer of cumulative colony posteriors
#################################################################

# # get posterior draws
# all_delta_mix = as.data.frame(cbind(lon = 1000*unlist(rstan::extract(mix.fit, pars = "delta_x")),
#                     lat = 1000*unlist(rstan::extract(mix.fit, pars = "delta_y"))))
# all_delta_imp = as.data.frame(cbind(lon = 1000*unlist(rstan::extract(imp.fit, pars = "delta_x")),
#                                     lat = 1000*unlist(rstan::extract(imp.fit, pars = "delta_y"))))
#
#
# # convert posterior draws to spatvector
# deltavector_mix = terra::vect(all_delta_mix[,], geom=c("lon", "lat"), crs=crs(landscape_raster), keepgeom=FALSE)
# deltavector_imp = terra::vect(all_delta_imp[,], geom=c("lon", "lat"), crs=crs(landscape_raster), keepgeom=FALSE)
#
#
# # create and fill raster
# rmix_empty = rast(
#   xmin = xmin(deltavector_mix),
#   xmax = xmax(deltavector_mix),
#   ymin = ymin(deltavector_mix),
#   ymax = ymax(deltavector_mix),
#   resolution = 30,
#   crs = crs(deltavector_mix))
# rimp_empty = rast(
#   xmin = xmin(deltavector_imp),
#   xmax = xmax(deltavector_imp),
#   ymin = ymin(deltavector_imp),
#   ymax = ymax(deltavector_imp),
#   resolution = 30,
#   crs = crs(deltavector_imp))
# values(rmix_empty) = 0
# value(rimp_empty) = 0
# deltavector_mix$count = 1
# deltavector_imp$count = 1
#
# rmix = rasterize(
#   deltavector_mix,
#   rmix_empty,
#   field = "count",
#   fun = "sum",
#   background = 0)
# plot(rmix)
#
#
# rimp = rasterize(
#   deltavector_imp,
#   rimp_empty,
#   field = "count",
#   fun = "sum",
#   background = 0)
# plot(rimp)


#################################################################
# Check out log-predictive density and perform LOOCV
#################################################################

# mixnormmin = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax1.23.rds")
# mixnormmed = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax2.68.rds")
# mixnormmax = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax4.14.rds")
#
# log_lik_draws1 = rstan::extract(mixnormmin, pars = "loglik")$loglik[,]
# log_lik_1 = data.frame(iter = 1:4000,
#                            log_lik = rowSums(log_lik_draws),
#                        model = "min")
# log_lik_draws2 = rstan::extract(mixnormmed, pars = "loglik")$loglik[,]
# log_lik_2 = data.frame(iter = 1:4000,
#                        log_lik = rowSums(log_lik_draws),
#                        model = "med")
# log_lik_draws3 = rstan::extract(mixnormmax, pars = "loglik")$loglik[,]
# log_lik_3 = data.frame(iter = 1:4000,
#                        log_lik = rowSums(log_lik_draws),
#                        model = "max")
#
# log_lik_df = rbind(log_lik_1, log_lik_2, log_lik_3)
#
#
# ggplot(log_lik_df, aes(x = iter, y = log_lik, color = model)) +
#   geom_line() +
#   theme_bw()
#
# min = loo(mixnormmin, pars = "loglik")
# med = loo(mixnormmed, pars = "loglik")
# max = loo(mixnormmax, pars = "loglik")
