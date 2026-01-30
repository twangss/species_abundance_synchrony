library(sf)
library(gdistance)
library(tidyverse)
library(stars)
library(terra)
library(raster)

#==== Get abundance coordinates and taxa data from model selection file
pinn_df=read.csv("data/Pinniped database.csv")
turt_df_long <- pinn_df %>% arrange(ID) %>%
  filter(Use == "yes") %>%
  mutate(CommonName = Common.Name, Site = Location.of.population, RMU = Subspecies) %>%
  mutate(ts_id = paste(CommonName,"#",Site)) %>% 
  pivot_longer(
    cols = starts_with("y"),
    names_to = "Year",
    names_prefix = "y",
    values_to = "Count",
    values_drop_na = TRUE) %>%
  mutate(Count=as.numeric(Count),Year=as.numeric(Year)) %>%
  group_by(ID) %>%
  drop_na(Count) %>%
  mutate(n=n()) %>% filter(n>3) # MUST BE THREE 

turt_IDs = turt_df_long %>% ungroup() %>% 
  dplyr::select(ID, Latitude,Longitude,CommonName,RMU) %>% distinct() 

turt_IDs = turt_df_long %>% 
  filter(Indiv.sum==1) %>% # only look at counted populations (states), avoid double counting multiple time series for a single location
  ungroup() %>% # turt_df_long
  # changed these selects for turtles 2025?
  dplyr::select(ID, Latitude,Longitude,CommonName,RMU) %>% distinct() #%>%

#==== Set up the distance shapefiles & rasters
sf_use_s2(T)
ocean = st_read("data/Shapefile/ne_10m_ocean/ne_10m_ocean.shp") %>%
  st_set_crs(4326) %>%
  mutate(value = 1)
# Rasterize the ocean coverage for gdistance package
r <- rast(ext(ocean), nrow = 1000, ncol = 1000)
ocean_conductive_rast <- raster::rasterize(vect(ocean), r, field=1,background = 0) 
ocean_rastlayer = raster(ocean_conductive_rast)
plot(ocean_rastlayer)
ocean_transition <- gdistance::transition(ocean_rastlayer, transitionFunction = mean, directions = 8,
                                          symm= T)
ocean_transition <- geoCorrection(ocean_transition, type="c")

#==== Calculate shortest ocean path length
distance_df = NULL
length_shapefiles = NULL
Name_RMU = turt_IDs %>% dplyr::select(CommonName,RMU) %>%
  distinct()
i=20
for (i in 1:nrow(Name_RMU)){ 
  CommonName_tmp = Name_RMU$CommonName[i]
  RMU_tmp = Name_RMU$RMU[i]
  cat(CommonName_tmp,"\n",RMU_tmp,"\n")
  species_i_coords = turt_IDs %>% 
    filter(CommonName == CommonName_tmp,
           RMU == RMU_tmp)
  
  nsites = nrow(species_i_coords)
  
  turt_df_species_IDs_tmp = species_i_coords %>% pull(ID) %>% 
    # sort for easier matching in for loop later
    sort() %>% unique()
  nsites = length(turt_df_species_IDs_tmp)
  
  # distance_matrix = matrix(data = NA, nrow = nsites, ncol = nsites)
  
  for(site_i in 1:nsites){
    # distance_matrix[site_i,site_i] = 0
    if(site_i == nsites){next}
    for (site_j in (site_i + 1):nsites) {
    
      site_i_ID = turt_df_species_IDs_tmp[site_i]
      site_j_ID = turt_df_species_IDs_tmp[site_j]
      
      turt_df_long_i = species_i_coords %>% filter(ID == site_i_ID) 
      turt_df_long_j = species_i_coords %>% filter(ID == site_j_ID) 
    
      sf_use_s2(T)
      site_i_coords = turt_df_long_i %>% dplyr::select(Longitude,Latitude) %>% as.matrix()
      site_j_coords = turt_df_long_j %>% dplyr::select(Longitude,Latitude) %>% as.matrix()
      
      distance <- shortestPath(ocean_transition, 
                               site_i_coords,site_j_coords, 
                               output = "SpatialLines")
      
      distance_sf = st_as_sf(distance) %>% st_set_crs(4326) %>%
        mutate(row_id = site_i_ID,col=site_j_ID)
      distance_sf = st_wrap_dateline(distance_sf) 
      
      distance_tmp = st_length(distance_sf)
      
      distance_df_tmp = data.frame(row_id = site_i_ID, col_id=site_j_ID, 
                              CommonName =  CommonName_tmp, RMU = RMU_tmp, distance_m = distance_tmp)
      
      length_shapefiles = rbind(length_shapefiles,distance_sf)
      
      distance_df=rbind(distance_df,distance_df_tmp)
      
    }
  }
}

#==== Save out intermediate step for distance data & distance shapefiles for plotting
# saveRDS(distance_df, "data/site_distance/seal_distance_long_subsp_2025_1025.RDS")
# st_write(length_shapefiles, "data/site_distance/seal_length_shapefiles_subsp_2025_1025.shp", append=F)

#==== Set up Distance Matrix format for MARSS AND set up simulation covariance values
# distance_df = readRDS("data/site_distance/seal_distance_long_subsp_2025_1025.RDS")

# Simulation spatial synchrony true parameters
alpha_hat = 0.9
rho_hat = 200000
sd_hat = 0.11

CommonName_vec = sort(unique(distance_df$CommonName))
distance_matrix_list = NULL
covariance_matrix_list = NULL
CommonName_i=4
for (CommonName_i in 1:length(CommonName_vec)){
  CommonName_tmp = CommonName_vec[CommonName_i]
  cat(CommonName_tmp,"\n")
  
  species_i_coords = turt_IDs %>% 
    filter(CommonName == CommonName_tmp)
  nsites = nrow(species_i_coords)
  
  species_i_IDs = species_i_coords %>% pull(ID) %>% sort()
  covariance_matrix = matrix(data = 0, nrow = nsites, ncol = nsites)
  rownames(covariance_matrix) = species_i_IDs
  colnames(covariance_matrix) = species_i_IDs
  
  distance_matrix = matrix(data = -99, nrow = nsites, ncol = nsites)
  rownames(distance_matrix) = species_i_IDs
  colnames(distance_matrix) = species_i_IDs
  
  for(site_i in 1:nsites){
    covariance_matrix[site_i,site_i] = 1
    distance_matrix[site_i,site_i] = 0
    if(site_i == nsites){next}
    for (site_j in (site_i + 1):nsites) {
  
      site_i_ID = species_i_IDs[site_i]
      site_j_ID = species_i_IDs[site_j]
      
      distance_m_tmp = distance_df %>% filter(row_id == min(site_i_ID,site_j_ID),
                             col_id == max(site_i_ID,site_j_ID)) %>%
        pull(distance_m)
      
      distance_km_tmp = as.numeric(distance_m_tmp/1000)
      if(length(distance_m_tmp) == 0){next} #if theres no distance data skip
      covariance_matrix[site_i,site_j] = alpha_hat*exp((-1/2)*(1/rho_hat)*distance_km_tmp^2)
      covariance_matrix[site_j,site_i] = covariance_matrix[site_i,site_j]
      
      distance_matrix[site_i,site_j] = distance_km_tmp
      distance_matrix[site_j,site_i] = distance_matrix[site_i,site_j]
  
    }
  }
  distance_matrix[site_i,site_i] = 0
  # turn the correlation matrix into a covariance matrix
  covariance_matrix_list[[CommonName_i]] = covariance_matrix*sd_hat*sd_hat
  names(covariance_matrix_list)[CommonName_i] = CommonName_tmp
  
  distance_matrix_list[[CommonName_i]] = distance_matrix
  names(distance_matrix_list)[CommonName_i] = CommonName_tmp
}

# save out files: covariance matrix for simulation experiment, distance matrix for realworld model fitting 

# saveRDS(distance_matrix_list, "data/site_distance/seal_distance_matrix_subsp_2025_1025.RDS")
# distance_matrix_list = readRDS("data/site_distance/seal_distance_matrix_subsp_2025_1025.RDS") # TURTLES

# saveRDS(covariance_matrix_list, "data/site_distance/seal_covariance_matrix_subsp_2025_1025.RDS")
# covariance_matrix_list = readRDS("data/site_distance/seal_covariance_matrix_subsp_2025_1025.RDS")


#=== Simulate abundance trends for simulation experiment
start_pops_year1 = turt_df_long %>% 
  group_by(ID) %>%
  summarise(log_Countmean_y1 = log(mean(Count)))

N_sims = 10
sim_i = 1
GP= F
true_state_name = "Ueq" # dependent on GP = F or T 
for (sim_i in 1:N_sims){
  set.seed(123+sim_i)
Year_seq = min(turt_df_long$Year):max(turt_df_long$Year)
sim_pop_matrix = NULL
for ( i in 1:length(covariance_matrix_list)){
  print(i)
  covariance_matrix = covariance_matrix_list[[i]]
  covariance_matrix_r = covariance_matrix
  
  # TOGGLE ON FOLLOWING TO TURN OFF COVARIANCE
  if (GP==F){
  covariance_matrix_r = diag(diag(covariance_matrix))
  }
  
  # year 1
  row_IDs = as.numeric(rownames(covariance_matrix)) %>% sort()
  sim_pops_year_tmp = start_pops_year1 %>% filter(ID %in% row_IDs) %>% pull(log_Countmean_y1) %>%
    data.frame()
  
  growth_rates_tmp_tminus1 = as.vector(mvtnorm::rmvnorm(n=1,sigma = covariance_matrix_r))  
  # year 2 and on
  for (y_tmp in 1:(length(Year_seq)-1)){
    growth_rates_tmp = as.vector(mvtnorm::rmvnorm(n=1,
                                                  mean = rep(0,nrow(covariance_matrix_r)),
                                                  sigma = covariance_matrix_r))  
  
  sim_pops_year_tmp[,y_tmp+1] = sim_pops_year_tmp[,y_tmp]+growth_rates_tmp
  growth_rates_tmp_tminus1 = growth_rates_tmp 
  }
  
  sim_pop_matrix_tmp = as.matrix(sim_pops_year_tmp)
  rownames(sim_pop_matrix_tmp) = row_IDs
  colnames(sim_pop_matrix_tmp) = Year_seq
  
  sim_pop_matrix = rbind(sim_pop_matrix, sim_pop_matrix_tmp)
}

sim_pop_matrix_long = sim_pop_matrix %>% as.data.frame() %>%
  rownames_to_column("ID") %>%
  pivot_longer(!ID, names_to = "Year", values_to = "Count_log") %>%
  mutate(ID = as.numeric(ID),
         Year = as.numeric(Year),
         Count_log = as.numeric(Count_log))

turt_df_long_attribute = turt_df_long %>% 
  dplyr::select(ID, Common.Name, Subspecies,ts_id,Site,Indiv.pop.map,Indiv.sum) %>% distinct()
sim_pop_matrix_long_true = sim_pop_matrix_long %>% left_join(turt_df_long_attribute)

# plot it to see all the RW MA1
sim_pop_matrix_long_true %>% 
  ggplot(aes(x=Year, y=Count_log,color=ID,group=ID)) + 
  geom_line(alpha=0.8) +
  facet_wrap(~Common.Name,scales="free")

# join real observation years to simulation data
Obs_years_df = turt_df_long %>% dplyr::select(ID, Year) %>%
  mutate(Obs_year = T) 

sim_pop_matrix_long_true_obs = sim_pop_matrix_long_true %>% 
  left_join(Obs_years_df, by = join_by(ID, Year)) %>%
  ungroup()%>%
  rowwise()%>%
  mutate(Count_log_obs = rnorm(n=1,mean=Count_log,sd=0.3))

# plot it to see all the observation error and RW MA1
p = sim_pop_matrix_long_true_obs %>% 
  ggplot(aes(x=Year, y=Count_log_obs,color=ID,group=ID)) + 
  geom_line(alpha=0.8) +
  facet_wrap(~Common.Name,scales="free")

print(p)

save_file_name = paste("data/simulation/seal_sim_true_obsV6_",true_state_name,"_sim",sim_i,".RDS",sep="")
saveRDS(sim_pop_matrix_long_true_obs,save_file_name)
}
