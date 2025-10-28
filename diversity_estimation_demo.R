# Author: Murray SA Thompson
# Contact: murray.thompson@cefas.gov.uk
# Version: 1
# November 2024

# initial r scripts developed for: https://doi.org/10.1111/gcb.15443 and subsequently refined for:
# “Mapping benthic biodiversity to facilitate future sustainable development”. 
# Cooper, K.M., Thompson, M.S.A., Bolam, S.G., Peach, C.M., Webb, T.J., Downie, A-L. 2025. Ecosphere (in press)
# and 
# “From plankton to fish: 21st-century redistribution of marine biodiversity and the changing role of rare species” 
# Couce, E., Greig, L., Engelhard, G., Pinnegar, J., Cooper K.M., Hélaouët, P., Pecuchet, L., Peck, M.A., Lindegren, M., Thompson, M.S.A. (in prep). 

# example data originally derived from: https://doi.org/https://doi.org/10.14466/CefasDataHub.126 

# Clear global environment (if needed)
#rm(list=ls()); gc()

# set your working directory
setwd('your file path')

# Loading required packages
pkgs = c("geosphere", "lubridate", "sf", "tidyverse", "mapplots", "purrr", "iNEXT", "ggpubr", "parallel")
invisible(lapply(pkgs, library, character.only = TRUE))


#######################################################
# setting up the functions:

## data tidying function for fish survey data
tidy_fish_dat <- function(spp, firstyear, lastyear){
  
  ## Tidying the data (renaming variables, grouping, filtering, adding dates)
  dat <- spp %>%
    rename(month=Month, day=Day) %>%
    mutate(day=as.numeric(day),
           sample = paste(HaulID, round(latitude,2), round(longitude,2),
                          day, month, year, sep='_'),
           date = as.Date(paste(day, month, year, sep='/'), "%d/%m/%Y"),
           # raw count for these fish data generated as follows:   
           count = round(DensAbund_N_Sqkm/(1/SweptArea_KM2))) %>%
    group_by(sample, latitude, longitude, year, month, day, date, species) %>%
    summarise(count = sum(count)) %>%
    filter(count >= 1,
           !is.na(species), 
           year %in% (firstyear-1):(lastyear+1)) %>%
    ungroup()
}

get_samples <- function(tidy_sppdf, yr){
  
  # select samples within timeframe
  sub_yr = tidy_sppdf %>%
    filter(year == yr) 
  
  min_date = as.Date(min(sub_yr$date))-n_days
  max_date = as.Date(max(sub_yr$date))+n_days
  
  samples_yr = tidy_sppdf %>%
    filter(date >= min_date & date <= max_date)
  
  # If data can be used
  if(nrow(samples_yr) > 0) {
    
    # Pulling samples 
    samp_samps <- samples_yr %>% 
      distinct(sample) %>%
      pull()
    
    # Defining focal locations for year of interest
    samp_locs <- tidy_sppdf %>% 
      filter(year == yr) %>%
      distinct(sample) %>%
      pull()
    
    # Define all points to generate info from
    comp_1_dist <- tidy_sppdf %>% 
      filter(sample %in% unique(c(samp_locs, samp_samps))) %>%
      dplyr::select(longitude, latitude, sample, date) %>% 
      distinct() 
    
    # Making a list 
    comp_1_lon_lat <- comp_1_dist %>%
      rownames_to_column() %>%  
      dplyr::select(-c(date)) %>%
      pivot_longer(c(-rowname, -sample)) %>% 
      dplyr::select(-rowname) %>% 
      pivot_wider(names_from = sample, values_from = value) %>%
      dplyr::select(-name) %>% 
      as.list()
    
    # Creating sample dates dataframe
    sample_dates <- comp_1_dist %>%
      distinct(sample, date) 
    
    # List of dataframes to return from the function
    ret <- list(
      focal_year = yr,
      sample_dates = sample_dates,
      samp_locs = samp_locs,
      comp_1_dist = comp_1_dist,
      comp_1_lon_lat = comp_1_lon_lat,
      samp_samps = samp_samps
    )
    
  } else {
    # Creating a NULL object to return if samples weren't found
    ret <- NULL
    print(paste(yr, 'missing temporal data', sep=' '))
  }
  
  ret
}

# Pairwise differences with sub-sampling option----
get_pairwise_diffs <- function(tidy_sppdf, yr, n_days, radius, n_samp, sub_sample){
  
  
  # Run the samples_yr function to get the necessary set-up for the focal year
  samples_yr <- get_samples(tidy_sppdf, yr = yr)
  
  # Create a distance matrix
  dist_days <- difftime(samples_yr$sample_dates$date,
                        as.Date("2000/01/01"), units = 'days') %>% 
    dist(diag = TRUE, upper = TRUE) %>% 
    as.matrix() %>% 
    as_tibble() %>% 
    setNames(samples_yr$sample_dates$sample) %>% 
    mutate(sample = samples_yr$sample_dates$sample) %>% 
    dplyr::select(sample, everything())
  
  # Remove duplicates 
  dist_days[upper.tri(dist_days)]  <- NA
  
  # Pivot data
  samps_days  <- dist_days %>%
    pivot_longer(cols = -c(sample),
                 names_to = "sample2",
                 values_to = "time_days") %>%
    filter(!is.na(time_days)) %>%
    mutate(time_days = abs(time_days))
  
  # Check locs are in sample and samps in sample2
  samps_in_date1 <- samps_days %>%
    filter(sample %in% samples_yr$samp_locs &
             sample2 %in% samples_yr$samp_samps) %>%
    mutate(sample_locs = sample,
           sample_samps = sample2)
  samps_in_date2 <- samps_days %>%
    filter(sample %in% samples_yr$samp_samps &
             sample2 %in% samples_yr$samp_locs) %>%
    mutate(sample_locs = sample2,
           sample_samps = sample)
  
  # Identify samps in required timeframe
  samps_in_date <-  samps_in_date1 %>%
    bind_rows(samps_in_date2) %>%
    filter(time_days <= n_days) %>%
    dplyr::select(-c(sample, sample2))
  
  # List of all samps within timeframe to generate spatial info
  samps <- unique(c(samps_in_date$sample_locs, samps_in_date$sample_samps))
  
  sub_comp_1_lon_lat  <- samples_yr$comp_1_lon_lat[samps]
  
  # Pairwise distances of all samps within timeframe 
  pairwise_distances <- map_df(sub_comp_1_lon_lat, function(x) {
    purrr::map(sub_comp_1_lon_lat, distCosine, x)}) %>%
    mutate(sample = names(sub_comp_1_lon_lat)) %>%
    dplyr::select(sample, everything())
  
  # Removing duplicates 
  pairwise_distances[upper.tri(pairwise_distances)] <- NA
  
  # Pivoting data
  samps_dists <-  try(pairwise_distances %>% 
                        pivot_longer(cols= -c(sample),
                                     names_to = "sample2",
                                     values_to = "distance_m") %>%
                        filter(!is.na(distance_m)), silent=T)
  
  if(is(samps_dists, "try-error") == FALSE) {
    
    # Checking locs are in sample and samps in sample2
    samps_in_region1 <- samps_dists %>%
      filter(sample %in% samples_yr$samp_locs &
               sample2 %in% samples_yr$samp_samps) %>%
      mutate(sample_locs = sample,
             sample_samps = sample2)
    samps_in_region2 <- samps_dists %>%
      filter(sample %in% samples_yr$samp_samps &
               sample2 %in% samples_yr$samp_locs) %>%
      mutate(sample_locs = sample2,
             sample_samps = sample)
    
    
    # Subset samples within radius
    samps_in_dist <- samps_in_region1 %>%
      bind_rows(samps_in_region2) %>%
      filter(distance_m < radius,
             distance_m > 0) %>%
      dplyr::select(-c(sample, sample2)) 
    
    # Identifying all samps in time-space window
    samps_in_dist_days <- samps_in_dist %>%
      left_join(samps_in_date, c('sample_locs', 'sample_samps')) %>%
      arrange(sample_locs) %>%
      filter(!is.na(time_days))  %>%
      distinct()  
    
    ret <- samps_in_dist_days
    
    # To randomly sub-sample based on n_samp
    if(sub_sample == 'Y') {
      
      # Create data.frame for sub-sampling using n_samp
      sub_samps_in_dist_days <- data.frame(sample_locs=NA,
                                           sample_samps=NA,
                                           distance_m=NA,
                                           time_days=NA)
      
      for(loc in unique(samps_in_dist_days$sample_locs)){ 
        
        df <-  subset(samps_in_dist_days, sample_locs == loc)
        
        samps <-  df$sample_samps
        
        set.seed(21)
        
        samps <- samps %>% sample(min(length(unique(df$sample_samps)), n_samp)) 
        
        sub_samps_in_dist_days <- sub_samps_in_dist_days %>%
          bind_rows(subset(df, sample_locs == loc & sample_samps %in% samps)) %>%
          filter(!is.na(sample_locs))
      }
      
      samps_in_dist_days <- sub_samps_in_dist_days
      
    } else{ 
      
      samps_in_dist_days <- samps_in_dist_days
      
    }
    
    # You may want to subset samples further based on n_samp in area using sub_samps_in_dist_days
    samps_rad_info <- samps_in_dist_days %>%
      group_by(sample_locs) %>%
      summarise(n_samps = length(unique(sample_samps))) %>%
      as.data.frame()
    
    # All samples to estimate from inc. years either side
    samples_yr <- tidy_sppdf %>%
      subset(year==yr | 
               (year==yr+1) |
               (year==yr-1)) %>%
      mutate(month = month(date)) %>%
      ungroup() 
    
    # Format data, and include loc to estimate gamma diversity
    divdf <- samples_yr %>%
      filter(sample %in% unique(c(samps_in_dist_days$sample_locs, samps_in_dist_days$sample_samps))) %>%
      group_by(sample, species) %>%
      summarise(count = round(sum(count))) %>%
      filter(count >0) %>%
      pivot_wider(names_from = sample,
                  values_from = count,
                  values_fn = sum,
                  values_fill = 0) %>%
      as.data.frame()
    rownames(divdf) = divdf[,1]
    divdf = divdf[,-1]
    
    # List of dataframes to return from function
    returnlist <- list("divdf" = divdf, 
                       "samples_yr" = samples_yr, 
                       "samps_rad_info" = samps_rad_info,
                       "samps_in_dist_days" = samps_in_dist_days) 
    
    return(returnlist)
    
    
  } else {
    ret <- NULL
  }
  
  ret
}

# diversity estimation function
div_est <-function(tidy_sppdf, bt, gamma_n_samp, yr, n_days, radius, n_samp, sub_sample){
  
  # select samples within timeframe
  sub_yr = tidy_sppdf %>%
    filter(year == yr) 
  
  min_date = as.Date(min(sub_yr$date))-n_days
  max_date = as.Date(max(sub_yr$date))+n_days
  
  samples_yr = tidy_sppdf %>%
    filter(date >= min_date & date <= max_date)
  
  ### Alpha diversity estimation for all samples 
  
  # Format data, and include loc to estimate gamma diversity
  a_divdf <- samples_yr %>%
    group_by(sample, species) %>%
    summarise(count = round(sum(count))) %>%
    filter(count >0) %>%
    pivot_wider(names_from = sample,
                values_from = count,
                values_fn = sum,
                values_fill = 0) %>%
    as.data.frame()
  rownames(a_divdf) = a_divdf[,1]
  a_divdf = a_divdf[,-1]
  
  # alpha diversity estimation using individual-based rarefaction and extrapolation
  myest_yr <- iNEXT(a_divdf, q=c(0,1,2), size=c(bt$b,bt$t), 
                    datatype="abundance", nboot = 50) 
  
  # sample-level alpha estimates
  a_samp_dat <- myest_yr$iNextEst$size_based %>%
    filter(m == bt$b,
           qD < bt$b) %>% ## new (replaces sample_div_yr)
    rename(sample=Assemblage) %>%
    mutate(hill_n = case_when(Order.q==0 ~ 'sample_a_q0',
                              Order.q==1 ~ 'sample_a_q1',
                              Order.q==2 ~ 'sample_a_q2')) %>%
    dplyr::select(sample, qD, hill_n) %>%
    pivot_wider(names_from = hill_n,
                values_from = qD)
  
  samples <- get_pairwise_diffs(tidy_sppdf, 
                                n_days=n_days, 
                                radius=radius, 
                                n_samp=n_samp, 
                                yr=yr, 
                                sub_sample=sub_sample)
  
  
  # Creating a dataframe to store the values
  samps_rad_info_yr <- tibble(year = yr)
  
  # Gather samples within radius, calculate sum of distances and days
  for(i in 1:nrow(samples$samps_rad_info)){#i=1
    loc <- unique(samples$samps_rad_info$sample_locs)[i]
    samps_in_region <- unique(subset(samples$samps_in_dist_days, sample_locs == loc)$sample_samps)
    # Sample distance and time clustering
    dddat <- samples$samps_in_dist_days %>%
      filter(sample_locs %in% loc,
             sample_samps %in% samps_in_region)
    
    # Sample-level abundance
    dat <- samples_yr %>%
      filter(sample %in% c(loc, samps_in_region)) %>%
      group_by(sample) %>%
      summarise(count = sum(count)) %>%
      ungroup() %>%
      summarise(av_count = mean(count),
                sd_count = sd(count),
                # Coefficient of variation (a type of beta-diversity for abundance)
                cv_count = sd(count)/mean(count),
                tot_count = sum(count)) %>%
      mutate(sample = loc)
    
    # Alpha diversity estimation
    a_dat <- myest_yr$iNextEst$size_based %>%
      filter(m == bt$b,
             qD < bt$b) %>% ## new (replaces sample_div_yr)
      filter(Assemblage %in% c(loc, samps_in_region)) %>%
      group_by(Order.q) %>%
      summarise(av_alpha = mean(qD)) %>%
      mutate(hill_n = case_when(Order.q==0 ~ 'q0',
                                Order.q==1 ~ 'q1',
                                Order.q==2 ~ 'q2')) %>%
      dplyr::select(-Order.q) %>%
      pivot_wider(names_from = hill_n,
                  values_from = av_alpha) 
    
    # Add to created dataframe
    samples$samps_rad_info[i, 'sum_km_between_samps'] = round(sum(dddat$distance_m)/1000, 1)
    samples$samps_rad_info[i, 'sum_time_between_samps'] = sum(dddat$time_days)
    samples$samps_rad_info[i, 'av_count'] = dat$av_count
    samples$samps_rad_info[i, 'sd_count'] = dat$sd_count
    samples$samps_rad_info[i, 'cv_count'] = dat$cv_count 
    samples$samps_rad_info[i, 'tot_count'] = dat$tot_count
    samples$samps_rad_info[i, 'a_q0'] = ifelse(nrow(a_dat)<1|is.null(a_dat$q0), NA, a_dat$q0)
    samples$samps_rad_info[i, 'a_q1'] = ifelse(nrow(a_dat)<1|is.null(a_dat$q1), NA, a_dat$q1)
    samples$samps_rad_info[i, 'a_q2'] = ifelse(nrow(a_dat)<1|is.null(a_dat$q2), NA, a_dat$q2)
    samples$samps_rad_info[i, 'alpha_n'] = bt$b ## new
    samples$samps_rad_info[i, 'gamma_n_samp'] = gamma_n_samp ## new
    
    if(length(samps_in_region)>1) {
      
      # Incidence matrix for sample-based rarefaction
      g_divdf <- samples$divdf[,c(loc, samps_in_region)]
      g_divdf <- g_divdf[rowSums(g_divdf) > 0, ]
      g_divdf[g_divdf>0] = 1
      
      samples$samps_rad_info[i, 'g_taxa_count'] = nrow(g_divdf) ##
      samples$samps_rad_info[i, 'year'] = yr ##
      
      # Gamma-diversity estimation using sample-based rarefaction and extrapolation
      gam_samp_yr <- try(iNEXT(list(g_divdf), q=c(0,1,2), datatype="incidence_raw", 
                               size=c(5, n_samp+1, gamma_n_samp), nboot = 50), silent = T)
      
      
      if(is(gam_samp_yr,"try-error")==F) {
        
        g_dat <- gam_samp_yr$iNextEst$size_based %>%
          filter(t == gamma_n_samp) %>%
          mutate(hill_n = case_when(Order.q == 0 ~ 'q0',
                                    Order.q == 1 ~ 'q1',
                                    Order.q == 2 ~ 'q2')) %>%
          dplyr::select(hill_n, qD) %>%
          pivot_wider(names_from = hill_n,
                      values_from = qD)
        
        if(nrow(g_dat)>0) {
          
          samples$samps_rad_info[i, 'g_q0'] = g_dat$q0
          samples$samps_rad_info[i, 'g_q1'] = g_dat$q1
          samples$samps_rad_info[i, 'g_q2'] = g_dat$q2 
          
        } else {
          
          samples$samps_rad_info[i, 'g_q0'] = NA
          samples$samps_rad_info[i, 'g_q1'] = NA
          samples$samps_rad_info[i, 'g_q2'] = NA 
          
          print('missing spatial samples for gamma estimates')
          
        }
        
      } else {
        
        samples$samps_rad_info[i, 'g_q0'] = NA
        samples$samps_rad_info[i, 'g_q1'] = NA
        samples$samps_rad_info[i, 'g_q2'] = NA 
        
        print('insufficient species occurrences for gamma estimate')
        
      }
      
      # Beta diversity estimation
      samples$samps_rad_info <- samples$samps_rad_info %>%
        mutate(b_q0 = if_else(is.na(g_q0)| is.na(a_q0), NA, g_q0/a_q0),
               b_q1 = if_else(is.na(g_q1)| is.na(a_q1), NA, g_q1/a_q1),
               b_q2 = if_else(is.na(g_q2)| is.na(a_q2), NA, g_q2/a_q2)) 
      
      samps_rad_info_yr <- samples$samps_rad_info
      
    }
    
  }
  
  if (exists("gam_samp_yr") && is.data.frame(gam_samp_yr)) {
    
    save(gam_samp_yr, myest_yr, file=paste0(outpath, 'biodiversity_estimates_additional_info_',yr,'_',file_name,'.Rdata'))
    
  } else {
    
    save(myest_yr, file=paste0(outpath, 'biodiversity_estimates_additional_info_',yr,'_',file_name,'.Rdata'))
    
  }
  
  # Store the data we want to use
  save(samps_rad_info_yr, a_samp_dat, file=paste0(outpath, 'biodiversity_estimates_',yr,'_',file_name,'.Rdata')) 
  
  return(glimpse(samps_rad_info_yr))
  
}


all_func <- function(yr){
  
  # Biodiversity estimations function
  div_est(tidy_sppdf, bt, gamma_n_samp, yr=yr, n_samp=n_samp, n_days=n_days, radius=radius, sub_sample=sub_sample)
  
}


# Parallel processing function
parallel_process <- function(){
  
  # Calculate the number of cores
  no_cores <- no_cores
  
  # Initiate cluster:
  cl <- makeCluster(no_cores)
  # Export packages:
  clusterEvalQ(cl, {   library(iNEXT); library(tidyverse); library(purrr); library(ggpubr);
    library(dplyr); library(lubridate); library(geosphere); library(janitor)})  
  # Export variables:
  clusterExport(cl, varlist=c('tidy_sppdf', 'years', 'all_func', 'get_samples',  
                              'get_pairwise_diffs', 'div_est', 'samp_info', 
                              'n_samp', 'bt', 'gamma_n_samp', 'sub_sample', 
                              'radius', 'n_days', 
                              'no_cores', 'outpath', 'file_name', 'years'))  
  ll <- parLapply(cl, years, function(yr) all_func(yr)) #tryCatch(all_func(yr), error = function(e) e) 
  stopCluster(cl)
  
}

combine_all_years <- function(outpath, years, file_name) {
  
  #### Now need to read all the annual files and rbind them together:
  setwd(outpath)
  
  # List all .RData files in the folder
  rdata_files <- list.files(pattern = "\\.Rdata$")
  
  # Check folder to see which years are available to load and exclude additional info files to load
  exclude_patterns <- c("additional_info", "combined", "final")
  rdata_files_to_load <- rdata_files[!grepl(paste(exclude_patterns, collapse = "|"), basename(rdata_files))]
  rdata_files_to_load <- rdata_files_to_load[grepl(paste(file_name, collapse = "|"), basename(rdata_files_to_load))]
  
  #### Now need to read all the annual files and rbind them together:
  load(rdata_files_to_load[1])
  
  all_df <- a_samp_dat %>%
    left_join(samps_rad_info_yr, c('sample'='sample_locs'))
  
  for(df in rdata_files_to_load[-1]) {
    
    load(df)
    
    yr_df = try(a_samp_dat %>%
                  left_join(samps_rad_info_yr, c('sample'='sample_locs')), silent = T)
    
    if(is(yr_df,"try-error")==F) {
      all_df <- all_df %>%
        bind_rows(yr_df)
    }
    
  }
  
  
  all_df <- all_df %>%
    dplyr::select(-year) %>%
    left_join(samp_info, c('sample')) 
  
  save(all_df, file=paste0(outpath, 'combined_biodiversity_estimates_', 
                           min(years), '_', max(years), '_', file_name, '.Rdata'))
  load(paste0(outpath, 'combined_biodiversity_estimates_', 
              min(years), '_', max(years), '_', file_name, '.Rdata'))
  
  return(all_df)
  
}


#######################################
# applying those diversity estimation functions to survey data

# load example fish survey data from otter trawls in the North Sea between 2010-2015 called 'spp' 
load('example_fish_survey_data.Rdata')

# check % identified to species
spp %>%
  mutate(taxa_level = case_when(!is.na(family) & is.na(genus) ~ 'family',
                                !is.na(genus) & is.na(species) ~ 'genus',
                                !is.na(species) ~ 'species',
                                TRUE ~ 'other'))  %>%
  count(taxa_level) %>% 
  mutate(prop = n / sum(n))

# name to be used when saving results 
file_name='GOV_fish'

# structure data correctly
tidy_sppdf <- tidy_fish_dat(spp, firstyear = min(spp$year), lastyear = max(spp$year))

##### set parameters for diversity assessment
# number of samples to use for gamma estimation
n_samp = 9
# whether to restrict to n_samp for gamma and beta 
#('Y' recommended because increased n samples tends to add novel beta-diversity, 
# systematically increasing the asymptote, see: https://doi.org/10.1111/gcb.15443)
sub_sample = 'Y'
# the number of days between samples
n_days = 182 # i.e. 6 months
# spatial radius in meters used to estimate gamma 
radius = 75000
# years to be assessed
years=sort(unique(tidy_sppdf$year)) 

# path to store results (ensure you have a folder to store the results in your working directory)
outpath = paste(getwd(), 'diversity_estimates', sep = '/')

# Calculate the number of cores for parallel processing
no_cores <- detectCores() - 2

# sample information and counts
samp_info = tidy_sppdf %>%
  group_by(sample, year, month, date, latitude, longitude) %>%
  summarise(sample_count = sum(count),
            log10_sample_count = log10(sum(count))) 

# Sampling depth for rarefaction = 2*mean sum of sample count across all data
bt = samp_info %>%
  ungroup() %>%
  summarise(b = round(mean(sample_count)*2),
            t = round(max(sample_count)))

# number of samples to extrapolate to for gamma-diversity estimates
gamma_n_samp <- (n_samp+1)*2

################################################################################

# Running the parallel_process function using the tidy data (can take a few minutes) 
parallel_process()

# load and combine all years together
all_df <- combine_all_years(outpath, years, file_name)
summary(all_df)

# justification for subsetting:
# alpha cutoff close to alpha_n (i.e. where all individuals in a sample are singletons), yielding unreliable estimates
# estimate beta- and gamma-diversity only used where n_samp == 9
# beta diversity >20 would represent >100% turnover, hence removed
# hill number diversity values should be 0>=1>=2 

final_df = all_df %>%
  mutate(sample_a_q0 = case_when(sample_a_q0 < alpha_n-1 ~ sample_a_q0, TRUE ~ NA),
         sample_a_q1 = case_when(sample_a_q1 < alpha_n-1 & sample_a_q1 <= sample_a_q0 ~ sample_a_q1,
                                 TRUE ~ NA),
         sample_a_q2 = case_when(sample_a_q2 < alpha_n-1 &
                                   sample_a_q2 <= sample_a_q1 ~ sample_a_q2,
                                 TRUE ~ NA),
         g_q0 = case_when(!is.na(sample_a_q0) & n_samp == 9 & b_q0 <= 20 ~ g_q0, TRUE ~ NA),
         g_q1 = case_when(!is.na(sample_a_q1) & n_samp == 9 & b_q1 <= 20 & g_q0 >= g_q1 ~ g_q1, TRUE ~ NA),
         g_q2 = case_when(!is.na(sample_a_q2) & n_samp == 9 & b_q2 <= 20 & g_q1 >= g_q2 ~ g_q2, TRUE ~ NA),
         # update beta-diversity estimates following NAs in alpha and beta above
         b_q0 = case_when(b_q0 <= 20 & n_samp == 9 & !is.na(g_q0) ~  b_q0, TRUE ~ NA),
         b_q1 = case_when(b_q1 <= 20 & n_samp == 9 & !is.na(g_q1) ~  b_q1, TRUE ~ NA),
         b_q2 = case_when(b_q2 <= 20 & n_samp == 9 & !is.na(g_q2) ~  b_q2, TRUE ~ NA),
         assemblage = file_name)


##### Producing basic maps using the outputs
# World map
world_shp = sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# Average gamma at one degree lat-lon level just for look see 
av_spatial_div = final_df %>% 
  mutate(log10_av_count = log10(av_count),
         log10_tot_count = log10(tot_count)) %>%
  dplyr::select(sample, latitude, longitude, year, month, date,
                log10_sample_count, cv_count, log10_tot_count, 
                sample_a_q0, sample_a_q1, sample_a_q2,
                b_q0, b_q1, b_q2,
                g_q0, g_q1, g_q2) %>%
  pivot_longer(cols= -c(sample:date),
               names_to = "metric",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(lat=round(latitude),
         lon=round(longitude),
         metric = factor(metric, levels = c('sample_a_q0','b_q0', 'g_q0', 
                                            'sample_a_q1', 'b_q1','g_q1',  
                                            'sample_a_q2','b_q2','g_q2',
                                            'log10_sample_count','cv_count','log10_tot_count'))) %>%
  group_by(lat, lon, metric) %>%
  summarise(av_metric = mean(value, na.rm = T))

labels_list<-list("sample_a_q0"= expression(alpha*", q"[0]),
                  "sample_a_q1"= expression(alpha*", q"[1]),
                  "sample_a_q2"= expression(alpha*", q"[2]),
                  "b_q0"= expression(beta*", q"[0]),
                  "b_q1"= expression(beta*", q"[1]),
                  "b_q2"= expression(beta*", q"[2]),
                  "g_q0"= expression(gamma*", q"[0]),
                  "g_q1"= expression(gamma*", q"[1]),
                  "g_q2"= expression(gamma*", q"[2]),
                  "log10_sample_count" = "Log10(sample count)",
                  "cv_count" = "Cv count",
                  "log10_tot_count" = "Log10(total count)")

plots_sp_div = list()
for(met in levels(av_spatial_div$metric)) { # met=levels(av_spatial_div$metric)[10] 
  
  p = av_spatial_div %>%
    filter(metric == met,
           !is.na(av_metric),
           !is.na(lat)) %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = av_metric)) +
    scale_fill_viridis_c(name = '',
                         breaks = scales::breaks_pretty(n = 3)) +
    labs(x='Longitude', y='Latitude', title=labels_list[[met]]) +
    guides(fill = guide_colourbar(barwidth = 1)) +
    scale_x_continuous(
      name = "Longitude",
      breaks = seq(-5, 10, by = 5)   # fewer ticks along x-axis
    ) +
    scale_y_continuous(
      name = "Latitude",
      breaks = seq(50, 62, by = 5)
    ) +
    geom_sf(data = world_shp, 
            fill = 'black', 
            color = 'black',
            size = 0.1) +
    coord_sf(xlim = c(-5, 10), ylim = c(49, 62)) +
    theme(axis.text.x = element_text(size = 8),  # smaller labels
          axis.text.y = element_text(size = 8),
          plot.margin = margin(2, 2, 2, 2),
          panel.background = element_rect(fill = 'grey80'),
          panel.border = element_rect(colour='black', fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plots_sp_div[[met]] = p
}

sp_div_plts = do.call(ggarrange, c(plots_sp_div, align='hv', ncol=3, nrow=4))

# Final plot for all metrics
#pdf(paste0(outpath,"/", "spatial_diverity_", file_name, ".pdf"), width = 8, height = 10, pointsize = 12)
sp_div_plts
#dev.off()

