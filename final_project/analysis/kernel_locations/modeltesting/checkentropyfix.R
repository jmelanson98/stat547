###### Fixing entropy issues
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
source("src/analysis_functions.R")
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

###############################################
# Possible "fixes" to test
###############################################

stanfile = c("models/simple_multinomial_uniformdisprior.stan",
             "models/simple_multinomial_normaldisprior.stan")
Rmax = c(1.3, 2.7, 4.15)
deltaprior = NA
params1 = expand.grid(stanfile = stanfile, Rmax = Rmax, deltaprior = deltaprior)

stanfile = c("models/simple_multinomial_nointercepts.stan",
             "models/simple_multinomial.stan")
Rmax= NA
deltaprior = NA
params2 = expand.grid(stanfile = stanfile, Rmax = Rmax, deltaprior = deltaprior)

stanfile = c("models/deltaprior.stan")
Rmax = NA
deltaprior = c(0.5,1,2)
params3 = expand.grid(stanfile = stanfile, Rmax = Rmax, deltaprior = deltaprior)

params = rbind(params1, params2, params3)

stanfilekey = data.frame(stanfile = c("models/simple_multinomial_uniformdisprior.stan",
                                      "models/simple_multinomial_normaldisprior.stan",
                                      "models/simple_multinomial_nointercepts.stan",
                                      "models/simple_multinomial.stan",
                                      "models/deltaprior.stan"),
                         key = c("simple_multinomial_uniformdisprior", 
                                 "simple_multinomial_normaldisprior", 
                                 "simple_multinomial_nointercepts", 
                                 "simple_multinomial", 
                                 "deltaprior"))
params = left_join(params, stanfilekey)


###############################################
# Load in and prepare data
###############################################

dataset = read.csv("analysis/kernel_locations/modeltesting/modeltest/simulation1_rho0.5.csv")
stan_data = prep_stan_simpleforaging_neutraldata(dataset)

if (params$stanfile[task_id] == "models/deltaprior.stan"){
  stan_data$colonybounds = NULL
  stan_data$deltaprior = params$deltaprior[task_id]
} else{
  stan_data$colonycenters = NULL
}

if (params$stanfile[task_id] %in% c("models/simple_multinomial_uniformdisprior.stan",
                               "models/simple_multinomial_normaldisprior.stan")){
  stan_data$Rmax = params$Rmax[task_id]
}

if (params$stanfile[task_id] == "models/simple_multinomial_uniformdisprior.stan"){
  stan_data$steepness = 50
}


###############################################
# Fit model
###############################################
stanfile = params$stanfile[task_id]
if (task_id %in% c(9,10,11)){
  fit = stan(file = stanfile,
                data = stan_data, seed = 5838299,
                chains = 4, cores = 4,
                iter = 3000, warmup = 1000,
                control = list(max_treedepth = 15),
                verbose = TRUE)
  draws = as.data.frame(rstan::extract(fit, pars = "rho"))
  write.csv(draws, paste0("analysis/kernel_locations/modeltesting/entropy/",
                          params$key[task_id],
                          "_Rmax", params$Rmax[task_id],
                          "_deltaprior", params$deltaprior[task_id],
                          "_draws.csv"))
}

draws = read.csv(paste0("analysis/kernel_locations/modeltesting/entropy/", 
                        params$key[task_id], 
                        "_Rmax", params$Rmax[task_id], 
                        "_deltaprior", params$deltaprior[task_id],
                        "_draws.csv"))

###############################################
# Perform profile likelihood
###############################################
rho = seq(0.15,5, by = 0.02)
profile = data.frame(rho = rho,
                     value = NA)
stanfile = paste0("models/fixedrho_", params$key[task_id], ".stan")
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
  print(paste0("Done with", profile$rho[i]))
}

profile = profile[5:nrow(profile),]
write.csv(profile, paste0("analysis/kernel_locations/modeltesting/entropy/", 
                                                   params$key[task_id], 
                                                   "_Rmax", params$Rmax[task_id], 
                                                   "_deltaprior", params$deltaprior[task_id],
                                                   "_profile.csv"))
profile = read.csv(paste0("analysis/kernel_locations/modeltesting/modeltesting/entropy/", 
                          params$key[task_id], 
                          "_Rmax", params$Rmax[task_id], 
                          "_deltaprior", params$deltaprior[task_id],
                          "_profile.csv"))

###############################################
# Now make some plots...
###############################################
xmax = max(draws$rho, 1)

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
col2 = rbind(g3, g4, size = "first")
task11=col2
saveRDS(col2, paste0("analysis/kernel_locations/modeltesting/modeltesting/entropy/", 
                                              params$key[task_id], 
                                              "_Rmax", params$Rmax[task_id], 
                                              "_deltaprior", params$deltaprior[task_id],
                                              "_grob.rds"))



toplist = list()
bottomlist = list()
for (i in 1:11){
  task_id = i
  profile = read.csv(paste0("analysis/kernel_locations/modeltesting/entropy/", 
                            params$key[task_id], 
                            "_Rmax", params$Rmax[task_id], 
                            "_deltaprior", params$deltaprior[task_id],
                            "_profile.csv"))
  draws = read.csv(paste0("analysis/kernel_locations/modeltesting/entropy/", 
                          params$key[task_id], 
                          "_Rmax", params$Rmax[task_id], 
                          "_deltaprior", params$deltaprior[task_id],
                          "_draws.csv"))
  xmax = max(draws$rho, na.rm = TRUE, 1)
  toplist[[i]] = ggplot(draws, aes(x = rho)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, color = "darkorchid4", fill = "lavender") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    xlab("") +
    ylab("") +
    xlim(c(0.25,xmax)) +
    theme_bw()
  bottomlist[[i]] = ggplot(profile, aes(x = rho, y = value)) +
    geom_line(color = "goldenrod", linewidth = 1.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    xlab(expression(rho)) +
    xlim(c(0.25,xmax)) +
    ylab("") +
    theme_bw()
}


toplist[[8]] = toplist[[8]] + 
  ggtitle("initial model") +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5)) +
ylab("posterior density")
bottomlist[[8]] = bottomlist[[8]] + 
  ylab("profile log posterior")
toplist[[7]] = toplist[[7]] + 
  ggtitle(expression("without "*epsilon[k])) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5))
toplist[[3]] = toplist[[3]] + 
  ggtitle(expression(R[max]*" = 2.7 (Uniform)")) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5))
toplist[[4]] = toplist[[4]] + 
  ggtitle(expression(R[max]*" = 2.7 (Normal)")) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5))
toplist[[10]] = toplist[[10]] + 
  ggtitle(expression(sigma[d]*" = 1")) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5))


# ---- Select plots ----
top <- toplist[c(8,7,3,4,10)]
bottom <- bottomlist[c(8,7,3,4,10)]

# ---- Tighten margins ----
top <- lapply(top, function(p) {
  p + theme(plot.margin = margin(2, 2, 2, 2))
})
bottom <- lapply(bottom, function(p) {
  p + theme(plot.margin = margin(2, 2, 2, 2))
})

# ---- Convert to grobs ----
top_grobs <- lapply(top, ggplotGrob)
bottom_grobs <- lapply(bottom, ggplotGrob)

# ---- Labels (A)–(J) ----
labels <- paste0("(", LETTERS[1:10], ")")

# ---- Function to insert label into PANEL ----
add_label <- function(g, label) {
  
  panel <- g$layout[g$layout$name == "panel", ]
  
  gtable_add_grob(
    g,
    grobs = textGrob(label,
                     x = unit(-0.05, "npc"),   # move left outside panel
                     y = unit(1.05, "npc"),    # move slightly above panel
                     just = c("right", "bottom"),
                     rot = 0,                 # diagonal
                     gp = gpar(fontface = "bold", fontsize = 12)),
    t = panel$t,
    l = panel$l,
    z = Inf,
    clip = "off"   # IMPORTANT: allow drawing outside panel
  )
}

# ---- Build aligned columns ----
combined <- Map(function(tg, bg, i) {
  
  # Equalize widths so x-axes align visually
  max_widths <- unit.pmax(tg$widths, bg$widths)
  tg$widths <- max_widths
  bg$widths <- max_widths
  
  # Add labels inside panels
  tg <- add_label(tg, labels[i])
  bg <- add_label(bg, labels[i + 5])
  
  # Stack tightly
  arrangeGrob(tg, bg,
              ncol = 1,
              heights = c(1, 1),
              padding = unit(0, "line"))
  
}, top_grobs, bottom_grobs, seq_along(top_grobs))

# ---- Final layout ----
final_plot <- arrangeGrob(grobs = combined, ncol = 5)

ggsave("analysis/kernel_locations/modeltesting/report_figures/realentropyfix.png",
       final_plot,
       width = 2800,
       height = 1200,,
       units = "px",
       dpi = 300)

