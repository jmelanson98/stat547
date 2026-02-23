###### Functions used in analysis
#### October 29, 2025
#### J. Melanson



###################################################
# Function for creating raster from shapefile
##################################################
create.raster <- function(sf_object,
                          lc_types,
                          resolution) {
  
  library(terra)
  library(sf)
  
  # remove empty geometries
  sf_object <- sf_object[!st_is_empty(sf_object), ]
  sf_object[[lc_types]] <- as.factor(sf_object[[lc_types]])
  
  # convert to SpatVector
  vect_object <- vect(sf_object)
  
  # create template raster
  template <- rast(
    ext(vect_object),
    resolution = resolution,
    crs = crs(vect_object)
  )
  
  # rasterize using factor field directly
  landscape_raster <- rasterize(
    vect_object,
    template,
    field = lc_types
  )
  
  # write raster
  writeRaster(
    landscape_raster,
    paste0("~/projects/def-ckremen/melanson/landscape/full_rasters/FValley_", resolution, "res.tif"),
    overwrite = TRUE
  )
  
  plot(landscape_raster)
  
  return(landscape_raster)
}




###################################################
# Functions for calculating landscape metrics
##################################################

calculateIJI <- function(landcover.raster, 
                         site.shapefile,
                         buffer.sizes){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  # load packages
  library(landscapemetrics)
  
  # initiate landscape.dat dataframe
  landscape.dat = as.data.frame(c(site.shapefile$site_id))
  colnames(landscape.dat) = c("sample_pt")
  
  
  for (buffer in buffer.sizes){
    res1_lsm_l_iji <- sample_lsm(landcover.raster,
                                  y = site.shapefile,
                                  plot_id = site.shapefile$site_id,
                                  shape = "circle",
                                  size = buffer,
                                  what = "lsm_l_iji")
    res1_iji = res1_lsm_l_iji[,c("plot_id", "value")]
    buffer_iji = paste("landscape_iji_", buffer, sep = "")
    colnames(res1_iji) = c("sample_pt", buffer_iji)

    
    #add prop_seminat and iji for given buffer size to landscape.dat
    landscape.dat = left_join(landscape.dat, res1_iji, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(landscape.dat)
}



calculateSN <- function(landcover.raster, 
                         site.shapefile, 
                         landcover.classification, 
                         buffer.sizes){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  # load packages
  library(landscapemetrics)
  
  # initiate landscape.dat dataframe
  landscape.dat = as.data.frame(c(site.shapefile$site_id))
  colnames(landscape.dat) = c("sample_pt")
  
  
  for (buffer in buffer.sizes){
    #area of each class in each buffer zone
    res1_lsm_c_ca <- sample_lsm(landcover.raster,
                                y = site.shapefile,
                                plot_id = site.shapefile$site_id,
                                shape = "circle",
                                size = buffer,
                                return_raster = TRUE,
                                what = "lsm_c_ca")

    #set class IDs
    res1_lsm_c_ca = merge(res1_lsm_c_ca, landcover.classification[,c("lc_type", "class")], by = "class")
    colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] = "sample_pt"

    # get buffer areas
    res1_lsm_c_ca %>%
      group_by(sample_pt) %>%
      summarise(buffer_area=sum(value)) -> buffer_areas
    res1_lsm_c_ca = merge(res1_lsm_c_ca, buffer_areas, by = "sample_pt")
    
    # #calculate landscape seminatural area (hedgerow + blackberry + low edge + forest + wetland)
    res1_seminat = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow" | lc_type == "forest" | lc_type == "wetland")
    res1_seminat$split_prop = res1_seminat$value/res1_seminat$buffer_area
    res1_seminat %>%
      group_by(sample_pt) %>%
      summarize(prop_seminat = sum(split_prop)) -> res1_seminat_summary
    buffer_edge = paste("prop_seminat_", buffer, sep = "")
    colnames(res1_seminat_summary) = c("sample_pt", buffer_edge)
    
    #add prop_seminat for given buffer size to landscape.dat
    landscape.dat = left_join(landscape.dat, res1_seminat_summary, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(landscape.dat)
}



# function to compute IDW natural area for a single point

compute_idw_area = function(points, raster, buffer, rho) {
  # create buffer around point
  buf = vect(st_buffer(points, dist = buffer))
  
  # extract values and coordinates for points inside buffer
  vals_coords = terra::extract(raster, buf, cells = TRUE, xy = TRUE)
  vals = vals_coords[, "lc_type"]
  coords = vals_coords[, c("x", "y")]
  
  # get focal point coords
  pt_coords = st_coordinates(points)
  
  # calculate vector of distances from focal point
  d = sqrt((coords[,1] - pt_coords[1])^2 + (coords[,2] - pt_coords[2])^2)
  
  # compute weights
  w = exp(-0.5 * (d/rho)^2)
  
  # weighted sum
  pixel_area = res(raster)[1] * res(raster)[2]
  dw_area = sum(vals * w, na.rm = TRUE) * pixel_area
  
  return(dw_area)
  
}





#############################################
# Functions for foraging distance analysis
#############################################

# function to prep data for stan --- both years together
prep_stan_simpleforaging_bothyears = function(sibships1,
                                              sibships2,
                                              effort1,
                                              effort2,
                                              samplepoints){
  # Remove queens and males from sibships
  sibships1 = filter(sibships1, !str_detect(notes, "queen") & 
                       !str_detect(notes, "male"))
  sibships2 = filter(sibships2, !str_detect(notes, "queen") & 
                       !str_detect(notes, "male"))

  
  # Adjust sibship IDs
  sibships1$sibshipID = as.numeric(sibships1$sibshipID)
  sibships2$sibshipID = as.numeric(sibships2$sibshipID)
  sibships2$sibshipID = sibships2$sibshipID + max(sibships1$sibshipID)
  stansibkeys = data.frame(sibshipID = unique(c(sibships1$sibshipID, sibships2$sibshipID)),
                           stansibkey = 1:length(unique(c(sibships1$sibshipID, sibships2$sibshipID))))
  sibships1 = left_join(sibships1, stansibkeys)
  sibships2 = left_join(sibships2, stansibkeys)
  
  # Make trap locations into correct format for calculating distance
  # convert to meter based crs
  fv_points = st_read(paste0(bombus_path, "/landscape/fvbombus/fvbombus_points.shp"))
  fv_points$site_name[fv_points$site_name == "pit_meadows"] = "pitt_meadows"
  traps_sf_m = st_transform(fv_points, 32610)
  traps_m = data.frame(sample_pt = traps_sf_m$site_id,
                       site_name = traps_sf_m$site_name,
                       trap_x = st_coordinates(traps_sf_m)[,1],
                        trap_y = st_coordinates(traps_sf_m)[,2])
  traps_m$site_name = as.factor(traps_m$site_name)
  
  # Create site ids
  sibships1$site = as.factor(sibships1$site)
  sibships2$site = as.factor(sibships2$site)
  site_keys = data.frame(site = levels(sibships1$site),
                         site_name = levels(traps_m$site_name),
                         site_id = 1:length(unique(sibships1$site)))
  
  # calculate sample effort per trap and year
  effort_sum1 = effort1 %>%
    filter(round %in% sibships1$round | round %in% sibships2$round) %>%
    group_by(sample_point, year) %>%
    summarize(total_effort = 5*n())
  effort_sum2 = effort2 %>%
    filter(round %in% sibships1$round | round %in% sibships2$round) %>%
    group_by(sample_point, year) %>%
    summarize(total_effort = 5*n())
  effort_sum = rbind(effort_sum1, effort_sum2)
  traps_m = traps_m %>%
    left_join(effort_sum, by = c("sample_pt" = "sample_point")) %>%
    left_join(site_keys, by = "site_name") %>%
    arrange(site_id)
  traps_m = traps_m[!is.na(traps_m$total_effort),]
  traps_m$trap_x = as.numeric(traps_m$trap_x)/1000
  traps_m$trap_y = as.numeric(traps_m$trap_y)/1000
  traps_m$trap_id = 1:nrow(traps_m)
  
  # Make sibships x sites df
  cxl1 = sibships1 %>%
    distinct(site, stansibkey, year)
  cxl2 = sibships2 %>%
    distinct(site, stansibkey, year)
  cxl = rbind(cxl1, cxl2)
  
  # Join sibships x sites with traps
  cxl = cxl %>% full_join(traps_m, relationship = "many-to-many")
  
  # Get counts of bees in traps
  counts1 = sibships1 %>%
    filter(notes != "male") %>%
    group_by(site, sample_pt, year, stansibkey) %>%
    summarize(count=n()) %>%
    ungroup()
  counts2 = sibships2 %>%
    filter(notes != "male") %>%
    group_by(site, sample_pt, year, stansibkey) %>%
    summarize(count=n()) %>%
    ungroup()
  counts = rbind(counts1, counts2)
  
  
  filled_counts = cxl %>%
    full_join(counts) %>%
    mutate(count = replace_na(count, 0)) %>%
    arrange(stansibkey)
  filled_counts$yn = ifelse(filled_counts$count > 0, 1, 0)
  filled_counts$pointyear = paste0(filled_counts$year, filled_counts$sample_pt)
  
  # Get the start indices for observations of each siblingship
  lengths = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(length = n()) %>%
    arrange(stansibkey)
  starts = cumsum(c(1, lengths$length))[1:length(lengths$length)]
  
  # Get the site id for each colony
  colony_sites = distinct(cxl[c("stansibkey", "site_id")]) %>%
    arrange(stansibkey)
  
  # Get spatial bounds for each site
  colonybounds = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(lower_x = min(trap_x) -5,
              upper_x = max(trap_x) +5,
              lower_y = min(trap_y) -5,
              upper_y = max(trap_y) +5) %>%
    arrange(stansibkey)
  
  # Get colony centers for each site
  colonycenters = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(center_x = mean(trap_x),
              center_y = mean(trap_y)) %>%
    arrange(stansibkey)
  

  
  # Make data list for stan
  # DISTANCES IN KM, NOT METERS
  stan_data <- list(
    C = length(unique(filled_counts$stansibkey)),
    K = length(unique(filled_counts$trap_id)),
    O = nrow(filled_counts),
    starts = starts,
    lengths = lengths$length,
    trap_pos = cbind(filled_counts$trap_x, filled_counts$trap_y),
    colony_id = filled_counts$stansibkey,
    trap_id = filled_counts$trap_id,
    sample_effort = filled_counts$total_effort,
    y_obs = filled_counts$count,
    yn = filled_counts$yn,
    colonybounds = colonybounds[,2:5],
    colonycenters = colonycenters[,2:3]
  )
  
  out = list(stan_data, filled_counts, traps_m)
  
  return(out)
}



# function to neutral simulated data for stan
prep_stan_simpleforaging_neutraldata = function(filled_counts){

  # Adjust sibship IDs
  stansibkeys = data.frame(colonyid = unique(filled_counts$colonyid),
                           stansibkey = 1:length(unique(filled_counts$colonyid)))
  filled_counts = left_join(filled_counts, stansibkeys)
  
  # Make indicator for whether counts_ik > 0
  filled_counts$yn = ifelse(filled_counts$counts > 0, 1, 0)
  filled_counts = filled_counts %>%
    arrange(stansibkey)
  
  # Get the start indices for observations of each siblingship
  lengths = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(length = n()) %>%
    arrange(stansibkey)
  starts = cumsum(c(1, lengths$length))[1:length(lengths$length)]
  
  # Get spatial bounds for each site
  colonybounds = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(lower_x = min(trap_x) -5,
              upper_x = max(trap_x) +5,
              lower_y = min(trap_y) -5,
              upper_y = max(trap_y) +5) %>%
    arrange(stansibkey)
  
  # Get colony centers for each site
  colonycenters = filled_counts %>%
    group_by(stansibkey) %>%
    summarize(center_x = mean(trap_x),
              center_y = mean(trap_y)) %>%
    arrange(stansibkey)

  # Make data list for stan
  # DISTANCES IN KM, NOT METERS
  stan_data <- list(
    C = length(unique(filled_counts$stansibkey)),
    K = length(unique(filled_counts$trap_id)),
    O = nrow(filled_counts),
    starts = starts,
    lengths = lengths$length,
    trap_pos = cbind(filled_counts$trap_x, filled_counts$trap_y),
    colony_id = filled_counts$stansibkey,
    trap_id = filled_counts$trap_id,
    sample_effort = filled_counts$total_effort,
    y_obs = filled_counts$counts,
    yn = filled_counts$yn,
    colonybounds = colonybounds[,2:5],
    colonycenters = colonycenters[,2:3]
  )
  
  return(stan_data)
}

#############################################
# Functions to simulate neutral datasets
#############################################

##### Simulate multiple landscapes, simple kernel #####
draw_simple_true_sites = function(sample_size, # number of bees to sample
                                  number_colonies, # positive integer, total
                                  colony_sizes, #vector of length number_colonies
                                  trap_data, # data frame of true traps
                                  rho, # set param value
                                  theta, # floral attractiveness
                                  distance_decay # current options: "exponentiated_quadratic" and "exponential"
){
  ##### Define colony characteristics #####
  colonyid = 1:number_colonies
  N = sum(colony_sizes)
  colony_data = as.data.frame(cbind(colonyid, colony_sizes))
  colony_data$colony_x = NA
  colony_data$colony_y = NA
  colony_data$site = NA
  
  
  ##### Simulate colony locations #####
  count = 0
  for(s in unique(trap_data$site)){
    start = (count*number_colonies/6)+1
    end = (count+1)*number_colonies/6
    colony_data$colony_x[start:end] = runif(number_colonies/6, min(trap_data$trap_x[trap_data$site == s]) -2, max(trap_data$trap_x[trap_data$site == s]) +2)
    colony_data$colony_y[start:end] = runif(number_colonies/6, min(trap_data$trap_y[trap_data$site == s]) -2, max(trap_data$trap_y[trap_data$site == s]) +2)
    colony_data$site[start:end] = s
    count = count + 1
  }
  print("Colony simulation complete.")
  
  
  # optional plotting step to visualize traps and colonies
  a=ggplot() +
    geom_point(data = colony_data, aes(x = colony_x, y = colony_y), colour = "lightblue", size = 0.2) +
    geom_point(data = trap_data, aes(x = trap_x, y = trap_y), colour = "black", size = 1) +
    ylab("Northing") +
    xlab("Easting") +
    theme_minimal()
  
  # Add floral data to traps
  trap_data$fq = rnorm(nrow(trap_data), 0, 1)
  
  ##### Compute distance-based visitation #####
  matrix2022 = full_join(colony_data[seq(nrow(colony_data)) %% 2 != 0,], trap_data[trap_data$year == 2022,], relationship = "many-to-many")
  matrix2023 = full_join(colony_data[seq(nrow(colony_data)) %% 2 == 0,], trap_data[trap_data$year == 2023,], relationship = "many-to-many")
  full_matrix = rbind(matrix2022, matrix2023)
  full_matrix$dist_ik = sqrt((full_matrix$trap_x-full_matrix$colony_x)^2 + (full_matrix$trap_y-full_matrix$colony_y)^2)
  
  # compute visitation rate of colony at each trap
  if (distance_decay == "exponentiated_quadratic"){
    full_matrix$lambda_ik <- exp(-0.5 * (full_matrix$dist_ik / rho)^2 + theta*full_matrix$fq + log(full_matrix$total_effort))
  } else if (distance_decay == "exponential"){
    full_matrix$lambda_ik <- exp((-full_matrix$dist_ik / rho) + theta*full_matrix$fq + log(full_matrix$total_effort))
  } else {
    print("Sorry, not a valid decay function.")
  }
  
  ##### Start sampling #####
  full_matrix$counts = 0
  
  while(sum(full_matrix$counts) < sample_size){
    
    # set weights based on colony sizes
    full_matrix$w_i = full_matrix$colony_sizes/N
    
    # calculate Pr(s = k | s in kappa)
    
    #  first compute Pr(s = k)
    # to do this, sum over C:
    #### lambda_ik * w_i / D_i # without underlying resource landscape, D_i is a constant and can be dropped
    full_matrix$lambda_ik_scaled <- full_matrix$lambda_ik*full_matrix$w_i
    trap_probs = full_matrix %>%
      group_by(trap_id) %>%
      summarize(intensity = sum(lambda_ik_scaled))
    
    # calculate intensity at all traps
    kappa = sum(trap_probs$intensity)
    
    # prob of sampling from a particular trap given k in kappa
    trap_probs$prob = trap_probs$intensity/kappa
    
    # sample trap
    trap = sample(trap_probs$trap_id, size = 1, replace = TRUE, prob = trap_probs$prob)
    
    # calculate Pr(c = i | s = k)
    partial = full_matrix[full_matrix$trap_id==trap,]
    partial$colony_prob = partial$lambda_ik_scaled/trap_probs$intensity[trap_probs$trap_id==trap]
    
    colony = sample(partial$colonyid, size = 1, prob = partial$colony_prob)
    
    # record visitation event to yik
    full_matrix$counts[(full_matrix$colonyid == colony & full_matrix$trap_id == trap)] = full_matrix$counts[(full_matrix$colonyid == colony & full_matrix$trap_id == trap)] + 1
    
    
    #update ni and N
    N = N-1
    full_matrix$colony_sizes[full_matrix$colonyid == colony] = full_matrix$colony_sizes[full_matrix$colonyid == colony] - 1
    print(paste0("Now sampled ", sum(full_matrix$counts), " of ", sample_size, " bees."))
  }
  
  return(full_matrix)
}
  






####################################
# Load packages
###################################
load.packages <- function(){
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("here",
                "readr",
                "readxl",
                "writexl", 
                "openxlsx",
                "janitor",
                "filesstrings",
                "textshaping",
                "extrafont",
                "qrcode", 
                "lubridate", 
                "naniar", 
                "Hmisc",
                "sf",
                "maps",
                "broom",
                #"rgdal",
                "raster",
                "dismo",
                "terra",
                "spData",
                "tmap",
                "leaflet",
                #"rgeos",
                #"maptools",
                "mapview",
                "stars",
                "landscapemetrics",
                "ggrepel",
                "vctrs",
                "measurements", 
                "arsenal", 
                "vegan",
                "taxize",
                "iNEXT", 
                "tidyverse",
                "stringr",
                "ggplot2", 
                "purrr",
                "magrittr",
                "dplyr",
                "tidylog",
                "data.table",
                "diffr",
                "RColorBrewer")
  
  ipak(packages)
}
