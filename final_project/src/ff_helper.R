###### Helper functions for floral foraging paper
#### February 9, 2026
#### J. Melanson


#############################################
##### Simulate data
#############################################

draw_floral_true_sites = function(ncolonies, # positive integer, per year
                                  colony_sizes, #vector of length ncolonies
                                  rho0, # baseline foraging
                                  rho1, # floral attractiveness
                                  rho2, #landscape effect
                                  rho3, # interactions of local & landscape scale
                                  rhomax,
                                  effort1,
                                  effort2,
                                  spec1,
                                  spec2
){
  df_list = list()
  for (year in 1:2){
  ##### Define colony characteristics #####
    if (year == 1){
      colonyid = 1:ncolonies
    } else if (year == 2){
      colonyid = (ncolonies+1):(2*ncolonies)
    }
  
  N = sum(colony_sizes)
  colony_data = as.data.frame(cbind(colonyid, colony_sizes))
  colony_data$colony_x = NA
  colony_data$colony_y = NA
  colony_data$site = NA
  colony_data$year = year
  
  ##### Get trap data #####
  fv_points = st_read(paste0(bombus_path, "/landscape/fvbombus/fvbombus_points.shp"))
  traps_sf_m = st_transform(fv_points, 32610)
  trap_data = data.frame(sample_point = traps_sf_m$site_id,
                       trap_x = st_coordinates(traps_sf_m)[,1]/1000,
                       trap_y = st_coordinates(traps_sf_m)[,2]/1000,
                       site = traps_sf_m$site_name)
  trap_data$site[trap_data$site == "pit_meadows"] = "pitt_meadows"
  
  
  ##### Simulate colony locations #####
  count = 0
  for(s in unique(trap_data$site)){
    start = (count*ncolonies/6)+1
    end = (count+1)*ncolonies/6
    colony_data$colony_x[start:end] = runif(ncolonies/6, min(trap_data$trap_x[trap_data$site == s]) -2, max(trap_data$trap_x[trap_data$site == s]) +2)
    colony_data$colony_y[start:end] = runif(ncolonies/6, min(trap_data$trap_y[trap_data$site == s]) -2, max(trap_data$trap_y[trap_data$site == s]) +2)
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
  
  
  ##### Add landscape data to traps #####
  landscape_df = data.frame(site = unique(trap_data$site),
                            lq = rnorm(length(unique(trap_data$site)), 0, 1))
  trap_data = left_join(trap_data, landscape_df)
  
  ##### Compute distance between colonies and traps #####
  full_matrix = full_join(colony_data, trap_data, relationship = "many-to-many")
  full_matrix$dist_ik = sqrt((full_matrix$trap_x-full_matrix$colony_x)^2 + (full_matrix$trap_y-full_matrix$colony_y)^2)
  
  ##### Make long format wrt timepoints #####
  if (year == 1){
    effort = effort1[,c("sample_point", "round")]
    nsamples  = spec1 %>% 
      filter(!str_detect(notes, "queen") & !str_detect(notes, "male")) %>% 
      group_by(round) %>%
      summarize(n = n())
  } else if (year == 2){
    effort = effort2[,c("sample_point", "round")]
    nsamples  = spec2 %>% 
      filter(!str_detect(notes, "queen") & !str_detect(notes, "male")) %>% 
      group_by(round) %>%
      summarize(n = n())
  }
  
  long_matrix = right_join(full_matrix, effort, relationship = "many-to-many")
  
  ##### Add floral quality values #####
  long_matrix$fq = rnorm(nrow(long_matrix), 0, 1)
  long_matrix$counts = 0
  
  
  
  ##### Start sampling! #####
  for (t in unique(nsamples$round[nsamples$n > 0])){
    print(paste0("Sampling round ", t))
    # take a submatrix for this timepoint
    sub_matrix = long_matrix[long_matrix$round == t,]
    
    # compute visitation rate of colony at each trap
    sub_matrix$lambda_ik = exp(-sub_matrix$dist_ik^2 / 
                                   (rhomax*inv.logit(rho0 + rho1*sub_matrix$fq + rho2*sub_matrix$lq + rho3*sub_matrix$fq*sub_matrix$lq)))
    
    for (bee in 1:nsamples$n[nsamples$round == t]){
      
      # set weights based on colony sizes
      sub_matrix$w_i = sub_matrix$colony_sizes/N
      
      # calculate Pr(s = k | s in kappa)
      
      #  first compute Pr(s = k)
      # to do this, sum over C:
      #### lambda_ik * w_i / D_i # without underlying resource landscape, D_i is a constant and can be dropped
      sub_matrix$lambda_ik_scaled = sub_matrix$lambda_ik*sub_matrix$w_i
      trap_probs = sub_matrix %>%
        group_by(sample_point) %>%
        summarize(intensity = sum(lambda_ik_scaled))
      # calculate intensity at all traps
      kappa = sum(trap_probs$intensity)
      if (kappa == 0 || is.na(kappa)) {
        stop("kappa is zero or NA")
      }
      
      
      # prob of sampling from a particular trap given k in kappa
      trap_probs$prob = trap_probs$intensity/kappa

      
      # sample trap
      trap = sample(trap_probs$sample_point, size = 1, replace = TRUE, prob = trap_probs$prob)
      
      # calculate Pr(c = i | s = k)
      partial = sub_matrix[sub_matrix$sample_point==trap,]
      partial$colony_prob = partial$lambda_ik_scaled/trap_probs$intensity[trap_probs$sample_point==trap]
      
      colony = sample(partial$colonyid, size = 1, prob = partial$colony_prob)
      
      # record visitation event to yik in full matrix
      long_matrix$counts[(long_matrix$colonyid == colony & 
                            long_matrix$sample_point == trap & 
                            long_matrix$round == t)] = 
        long_matrix$counts[(long_matrix$colonyid == colony & 
                              long_matrix$sample_point == trap & 
                              long_matrix$round == t)] + 1
      
      
      #update ni and N
      N = N-1
      sub_matrix$colony_sizes[sub_matrix$colonyid == colony] = sub_matrix$colony_sizes[sub_matrix$colonyid == colony] - 1
      long_matrix$colony_sizes[long_matrix$colonyid == colony] = long_matrix$colony_sizes[long_matrix$colonyid == colony] - 1
      print(paste0("Now sampled ", sum(long_matrix$counts), " of ", nrow(spec1) + nrow(spec2), " bees."))
    }
  }
  df_list[[year]] = long_matrix
  }
  
  bothyears = rbind(df_list[[1]], df_list[[2]])
  return(bothyears)
}




#############################################
##### Prep data for stan
#############################################

# Function to prep *simulated data* for Stan
prep_stan_floralforaging_simdata = function(simdataset, rhomax){

  # Adjust sibship IDs
  simdataset$colonyid = as.numeric(simdataset$colonyid)
  stansibkeys = data.frame(colonyid = unique(c(simdataset$colonyid)),
                           stansibkey = 1:length(unique(c(simdataset$colonyid))))
  simdataset = left_join(simdataset, stansibkeys)
  
  # For each colony, remove rounds where we do not observe it
  simdataset = simdataset %>%
    group_by(stansibkey, round) %>%
    filter(sum(counts) > 0) %>%
    arrange(round, stansibkey, sample_point)
  
  # Make trap x year column
  trapyear = data.frame(trapyear = unique(paste0(simdataset$year, simdataset$sample_point)),
                        trapid = 1:length(unique(paste0(simdataset$year, simdataset$sample_point))))
  simdataset$trapyear = paste0(simdataset$year, simdataset$sample_point)
  simdataset = left_join(simdataset, trapyear) %>%
    arrange(round, stansibkey, trapid)
  
  # Get start of observations and number of observations per colony x timepoint combo
  # Get the start indices for observations of each siblingship
  lengths = simdataset %>%
    group_by(round, stansibkey) %>%
    summarize(length = n()) %>%
    arrange(round, stansibkey)
  starts = cumsum(c(1, lengths$length))[1:length(lengths$length)]
  
  # Get upper and lower bounds per site
  sitekey = data.frame(site = unique(simdataset$site),
                       siteid = 1:length(unique(simdataset$site)))
  simdataset = left_join(simdataset, sitekey)
  
  colonybounds = simdataset %>%
    group_by(stansibkey) %>%
    summarize(lower_x = min(trap_x) - 5,
              upper_x = max(trap_x) + 5,
              lower_y = min(trap_y) - 5,
              upper_y = max(trap_y) + 5) %>%
    arrange(stansibkey)
  
  # Get colony centers for each site
  colonycenters = simdataset %>%
    group_by(stansibkey) %>%
    summarize(center_x = mean(trap_x),
              center_y = mean(trap_y)) %>%
    arrange(stansibkey)
  
  
  # Get stan data!
  data = list(
    C = length(unique(simdataset$stansibkey)),
    K = length(unique(simdataset$trapid)),
    O = nrow(simdataset),
    CT = nrow(lengths),
    CTstarts = starts,
    CTlengths = lengths$length,
    trap_pos = cbind(simdataset$trap_x, simdataset$trap_y),
    colony_id = simdataset$stansibkey,
    trap_id = simdataset$trapid,
    fq = simdataset$fq,
    lq = simdataset$lq,
    yobs = simdataset$counts,
    colonybounds = colonybounds[,2:5],
    colonycenters = colonycenters[,2:3],
    rhomax = rhomax
    )
  
  return(data)
}

# Function to prep *real data* for Stan
prep_stan_floralforaging = function(sibships1,
                                    sibships2,
                                    effort1,
                                    effort2,
                                    veg1,
                                    veg2){
  
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
  
  # Get colonies x trap x timepoint matrix
  cxl1 = sibships1 %>%
    distinct(site, stansibkey, year)
  cxl2 = sibships2 %>%
    distinct(site, stansibkey, year)
  cxl = rbind(cxl1, cxl2)
  
  sereduced1 = effort1[(effort1$round %in% sibships1$round | effort1$round %in% sibships2$round),
                       colnames(effort1) %in% c("site", "sample_point", "round", "sample_id", "year")]
  sereduced2 = effort2[(effort2$round %in% sibships1$round | effort2$round %in% sibships2$round), 
                       colnames(effort2) %in% c("site", "sample_point", "round", "sample_id", "year")]
  sereduced = rbind(sereduced1, sereduced2)
  
  CKT = left_join(cxl, sereduced, by = c("site", "year"), relationship = "many-to-many")
  CKT$round = as.numeric(CKT$round)
  
  # Get site keys
  sitekey = data.frame(site = c("ED", "HR", "NR", "PM", "SD", "W"),
                       siteid = 1:6)
  
  # Get trap keys
  # here we want a different random intercept for each trap x year combo
  # trap "quality" may differ between year so we treat them as unique traps
  CKT$trapyear = paste0(CKT$sample_point, CKT$year)
  trapkey = data.frame(trapyear = sort(unique(CKT$trapyear)),
                       trap_id = 1:length(unique(CKT$trapyear)))
  
  
  # Get floral abundance estimates
  specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
  specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
  beeflowers = unique(c(specs2022$active_flower, specs2023$active_flower))
  beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]
  
  columnlist = c("sample_id", beeflowers)
  floral_df_wide = veg1[,colnames(veg1) %in% columnlist]
  floral_df_long1 = floral_df_wide %>%
    pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
    mutate(across(-c(sample_id, flower),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(-c(sample_id, flower),
                  ~ 10^(.-1))) %>%
    group_by(sample_id) %>%
    summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
    mutate(across(-c(sample_id),
                  ~ log10(.)))
  
  floral_df_wide = veg2[,colnames(veg2) %in% columnlist]
  floral_df_long2 = floral_df_wide %>%
    pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
    mutate(across(-c(sample_id, flower),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(-c(sample_id, flower),
                  ~ 10^(.-1))) %>%
    group_by(sample_id) %>%
    summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
    mutate(across(-c(sample_id),
                  ~ log10(.)))
  
  floral_df_long = rbind(floral_df_long1, floral_df_long2)
  
  
  # Get trap coordinates
  fv_points = st_read(paste0(bombus_path, "/landscape/fvbombus/fvbombus_points.shp"))
  traps_sf_m = st_transform(fv_points, 32610)
  traps_m = data.frame(sample_point = traps_sf_m$site_id,
                       trap_x = st_coordinates(traps_sf_m)[,1],
                       trap_y = st_coordinates(traps_sf_m)[,2],
                       site = traps_sf_m$site_name)
  traps_m$site[traps_m$site == "pit_meadows"] = "pitt_meadows"
  
  # Get Julian dates
  #add julian date to sample effort data frame
  effort1$date = paste(effort1$day, effort1$month, effort1$year)
  effort1$date = gsub(" ", "", effort1$date, fixed = TRUE)
  effort1$date <- as.POSIXlt(effort1$date, format = "%d%b%y")
  effort1$julian_date = effort1$date$yday
  
  effort2$date = paste(effort2$day, effort2$month, effort2$year)
  effort2$date = gsub(" ", "", effort2$date, fixed = TRUE)
  effort2$date <- as.POSIXlt(effort2$date, format = "%d%b%y")
  effort2$julian_date = effort2$date$yday
  
  effort = rbind(effort1[,colnames(effort1) %in% c("sample_id", "julian_date")], effort2[,colnames(effort2) %in% c("sample_id", "julian_date")])
  
  # Get nonzero counts
  counts1 = sibships1 %>%
    group_by(round, sample_pt, stansibkey) %>%
    summarize(count=n()) %>%
    ungroup()
  colnames(counts1) = c("round", "sample_point", "stansibkey", "counts")
  
  counts2 = sibships2 %>%
    group_by(round, sample_pt, stansibkey) %>%
    summarize(count=n()) %>%
    ungroup()
  colnames(counts2) = c("round", "sample_point", "stansibkey", "counts")
  
  counts = rbind(counts1, counts2)
  
  
  # Join!
  CKT = CKT %>%
    left_join(counts) %>%
    left_join(sitekey) %>%
    left_join(trapkey) %>%
    left_join(floral_df_long) %>%
    left_join(traps_m) %>%
    left_join(effort)
  
  # For each colony, remove rounds where we do not observe it
  CKT = CKT %>%
    group_by(stansibkey, round) %>%
    filter(sum(counts) > 0) %>%
    arrange(round, stansibkey, sample_point)
  
  # Clean up!
  CKT$counts[is.na(CKT$counts)] = 0 # these are true zeroes
  CKT$floral_abundance[CKT$floral_abundance == -Inf] = 0 #questionable -- not true zeroes
  CKT$floral_abundance[is.na(CKT$floral_abundance)] = 0 #questionable -- not true zeroes
  CKT$trap_x = as.numeric(CKT$trap_x)/1000
  CKT$trap_y = as.numeric(CKT$trap_y)/1000
  
  # Get start of observations and number of observations per colony x timepoint combo
  # Get the start indices for observations of each siblingship
  lengths = CKT %>%
    group_by(round, stansibkey) %>%
    summarize(length = n()) %>%
    arrange(round, stansibkey)
  starts = cumsum(c(1, lengths$length))[1:length(lengths$length)]
  
  
  # Get stan data!
  data = list(
    C = length(unique(CKT$stansibkey)),
    K = length(unique(CKT$trapyear)),
    O = nrow(CKT),
    CT = nrow(lengths),
    CTstarts = starts,
    CTlengths = lengths$length,
    trap_pos = cbind(CKT$trap_x, CKT$trap_y),
    colony_id = CKT$stansibkey,
    trap_id = CKT$trap_id,
    fq = scale(CKT$floral_abundance),
    yobs = CKT$counts,
    lower_x = (min(CKT$trap_x) - 5),
    upper_x = (max(CKT$trap_x) + 5),
    lower_y = (min(CKT$trap_y) - 5),
    upper_y = (max(CKT$trap_y) + 5)
  )
  
  return(list(data, CKT, traps_m))
}