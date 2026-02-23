###### Full posterior versus profile likelihood
### February 19, 2026
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
library(grid)
library(gtable)


setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
basepath = "analysis/kernel_locations/modeltesting/"
source("src/analysis_functions.R")

###############################################
# Load in a ! *totally perfect* ! dataset 
###############################################
#(meaning...waaaaay better than any data we could ever hope to achieve in the real world)

dataset1 = read.csv(paste0(basepath,"entropy/perfectdata.csv"))
stan_data = prep_stan_simpleforaging_neutraldata(dataset)
stan_data$colonycenters = NULL

#fit and save model (full posterior)
# stanfile = "models/simple_multinomial.stan"
# perfectfit = stan(file = stanfile,
#            data = stan_data, seed = 5838299,
#            chains = 4, cores = 4,
#            iter = 3000, warmup = 1000,
#            control = list(max_treedepth = 15),
#            verbose = TRUE)
# saveRDS(perfectfit, paste0(basepath, "entropy/perfectfit.rds"))
# draws = as.data.frame(rstan::extract(perfectfit, pars = "rho"))
# write.csv(draws, paste0(basepath, "entropy/perfectfit_draws.csv"))
draws = read.csv(paste0(basepath, "entropy/perfectfit_draws.csv"))

# perform profile likelihood
# rho = seq(0.15,2, by = 0.02)
# profile = data.frame(rho = rho,
#                      value = NA)
# stanfile = "models/fixedrho.stan"
# sm = stan_model(file = stanfile)
# 
# init_vals = NULL
# 
# for (i in 1:nrow(profile)) {
# 
#   stan_data$rho = profile$rho[i]
# 
#   init_fun = if (is.null(init_vals)) {"random"
#   } else {function() as.list(init_vals)}
# 
#   fit_opt = rstan::optimizing(
#     sm, data = stan_data,
#     init = init_fun,
#     algorithm = "LBFGS",
#     hessian = FALSE)
# 
#   profile$value[i] = fit_opt$value
#   init_vals <- fit_opt$par
# }
# profile = profile[4:nrow(profile),]
# write.csv(profile, paste0(basepath, "entropy/perfectprofile.csv"))
profile = read.csv(paste0(basepath, "entropy/perfectprofile.csv"))


# Now make some plots
xmax = max(draws$rho)
p1 = ggplot(draws, aes(x = rho)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, color = "darkorchid4", fill = "lavender") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab("") +
  ylab("posterior density") +
  xlim(c(0.25,xmax)) +
  theme_bw()

p2 = ggplot(profile, aes(x = rho, y = value)) +
  geom_line(color = "goldenrod", size = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab(expression(rho)) +
  xlim(c(0.25,xmax)) +
  ylab("profile log posterior") +
  theme_bw()

g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
g1$widths <- unit.pmax(g1$widths, g2$widths)
g2$widths <- unit.pmax(g1$widths, g2$widths)
col1 <- rbind(g1, g2, size = "first")



###############################################
# Load in a realistic simulated dataset 
###############################################

dataset2 = read.csv(paste0(basepath, "modeltest/simulation2_rho0.5.csv"))
stan_data = prep_stan_simpleforaging_neutraldata(dataset)
stan_data$colonycenters = NULL

#fit and save model (full posterior)
# stanfile = "models/simple_multinomial.stan"
# simfit = stan(file = stanfile,
#                   data = stan_data, seed = 5838299,
#                   chains = 4, cores = 4,
#                   iter = 3000, warmup = 1000,
#                   control = list(max_treedepth = 15),
#                   verbose = TRUE)
# saveRDS(simfit, paste0(basepath,"entropy/simfit.rds"))
# draws = as.data.frame(rstan::extract(simfit, pars = "rho"))
# write.csv(draws, paste0(basepath,"entropy/simfit_draws.csv"))
draws = read.csv(paste0(basepath,"entropy/simfit_draws.csv"))


# perform profile likelihood
rho = seq(0.15,2, by = 0.05)
profile = data.frame(rho = rho,
                     value = NA)
stanfile = "models/fixedrho.stan"
sm = stan_model(file = stanfile)

init_vals = NULL

for (i in 1:nrow(profile)) {

  stan_data$rho = profile$rho[i]

  init_fun = if (is.null(init_vals)) {"random"
  } else {function() as.list(init_vals)}

  fit_opt = rstan::optimizing(
    sm, data = stan_data,
    init = init_fun,
    algorithm = "LBFGS",
    hessian = FALSE)

  profile$value[i] = fit_opt$value
  init_vals <- fit_opt$par
}
profile = profile[4:nrow(profile),]
write.csv(profile, paste0(basepath, "entropy/simprofile.csv"))
profile = read.csv(paste0(basepath, "entropy/simprofile.csv"))

# Now make some plots
xmax = max(draws$rho)
p3 = ggplot(draws, aes(x = rho)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, color = "darkorchid4", fill = "lavender") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab("") +
  ylab("") +
  xlim(c(0.25,xmax)) +
  theme_bw()

p4 = ggplot(profile, aes(x = rho, y = value)) +
  geom_line(color = "goldenrod", size = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab(expression(rho)) +
  xlim(c(0.25,xmax)) +
  ylab("") +
  theme_bw()

g3 = ggplotGrob(p3)
g4 = ggplotGrob(p4)
g3$widths <- unit.pmax(g3$widths, g4$widths)
g4$widths <- unit.pmax(g3$widths, g4$widths)
col2 <- rbind(g3, g4, size = "first")

###############################################
# Load in a real dataset 
###############################################

# Load in data
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))

stan_data = prep_stan_simpleforaging_bothyears(impatiens_sibs2022,
                                          impatiens_sibs2023,
                                          effort2022,
                                          effort2023)
stan_data = stan_data[[1]]

#fit and save model (full posterior)
# stanfile = "models/simple_multinomial.stan"
# realfit = stan(file = stanfile,
#                   data = stan_data[[1]], seed = 5838299,
#                   chains = 4, cores = 4,
#                   iter = 3000, warmup = 1000,
#                   control = list(max_treedepth = 15),
#                   verbose = TRUE)
# saveRDS(realfit, "paste0(basepath, "entropy/realfit.rds"))
# draws = as.data.frame(rstan::extract(realfit, pars = "rho"))
# write.csv(draws, paste0(basepath, "entropy/realfit_draws.csv"))
draws = read.csv(paste0(basepath, "entropy/realfit_draws.csv"))

# perform profile likelihood
rho = seq(0.15,2, by = 0.05)
profile = data.frame(rho = rho,
                     value = NA,
                     hessian = NA
                     )
stanfile = "models/fixedrho.stan"
sm = stan_model(file = stanfile)

init_vals = NULL

for (i in 1:nrow(profile)) {
  
  stan_data$rho = profile$rho[i]
  
  init_fun = if (is.null(init_vals)) {"random"
  } else {function() as.list(init_vals)}
  
  fit_opt = rstan::optimizing(
    sm, data = stan_data,
    init = init_fun,
    algorithm = "LBFGS",
    hessian = FALSE)
  
  profile$value[i] = fit_opt$value
  init_vals <- fit_opt$par
  print(paste0("Done with ", profile$rho[i]))
}
profile = profile[4:nrow(profile),]
write.csv(profile, paste0(basepath, "entropy/realprofile.csv"))
profile = read.csv(paste0(basepath, "entropy/realprofile.csv"))

# Now make some plots
xmax = max(draws$rho)
p5 = ggplot(draws, aes(x = rho)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, color = "darkorchid4", fill = "lavender") +
  xlab("") +
  ylab("") +
  xlim(c(0.25,xmax)) +
  theme_bw()

p6 = ggplot(profile, aes(x = rho, y = value)) +
  geom_line(color = "goldenrod", size = 1.2) +
  xlab(expression(rho)) +
  xlim(c(0.25,xmax)) +
  ylab("") +
  theme_bw()

g5 = ggplotGrob(p5)
g6 = ggplotGrob(p6)
g5$widths <- unit.pmax(g5$widths, g6$widths)
g6$widths <- unit.pmax(g5$widths, g6$widths)
col3 <- rbind(g5, g6, size = "first")


# Make final plots
a = dataset1 %>% group_by(colonyid) %>% summarize(n=sum(counts))
b = dataset2 %>% group_by(colonyid) %>% summarize(n=sum(counts))
c = dataset3 %>% group_by(stansibkey) %>% summarize(n=sum(count))

plota = ggplot(a, aes(x = n)) +
  geom_histogram(color = "black", fill = "beige") +
  xlab("siblingship size") +
  ylab("frequency") +
  theme_bw()
plotb = ggplot(b, aes(x = n)) +
  geom_histogram(color = "black", fill = "beige") +
  xlab("siblingship size") +
  ylab("") +
  theme_bw()
plotc = ggplot(c, aes(x = n)) +
  geom_histogram(color = "black", fill = "beige") +
  xlab("siblingship size") +
  ylab("") +
  theme_bw()

final = cbind(col1, col2, col3, size = "first")
grid = grid.arrange(plota,plotb, plotc, ncol =3)
finalgrid = grid.arrange(grid, final, ncol = 1, heights = c(1,2))

finalgrid <- ggdraw() +
  draw_plot(finalgrid, 0, 0, 1, 0.97) +
  draw_plot_label(
    c("(A)", "(B)", "(C)",
      "(D)", "(E)", "(F)",
      "(G)", "(H)", "(I)"),
    
    x = c(0.02, 0.37, 0.70,
          0.02, 0.37, 0.70,
          0.02, 0.37, 0.70),
    
    y = c(0.98, 0.98, 0.98,
          0.65, 0.65, 0.65,
          0.33, 0.33, 0.33),
    
    hjust = 0,
    vjust = 1,
    size = 14,
    fontface = "bold"
  ) +
  theme(plot.background = element_rect(fill = "white", color = NA))


ggsave(paste0(basepath, "report_figures/fig1.png"), finalgrid, height = 2000, width = 2000, units = "px")
