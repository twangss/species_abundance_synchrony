library(tidyverse)
library(ggplot2)

### TRUE STATE SIMULATION PROCESSING
# read in the true simulated trajectories 
N_sims = 1
model_vec = c("Ueq", "GP")
all_tru_sim_df = NULL
for (sim_i in 1:10){
  for(model_i in model_vec){
    sim_file_name = paste("data/simulation/seal_sim_true_obsV6_",model_i,"_sim",sim_i,".RDS",sep="")
    sim_model_i = readRDS(sim_file_name) %>% mutate(sim = sim_i, true_state = model_i) 
    all_tru_sim_df=rbind(all_tru_sim_df,sim_model_i)
  }}

# just look at one sim if you want
model_i = "GP"
sim_i=1
sim_file_name = paste("data/simulation/seal_sim_true_obsV6_",model_i,"_sim",sim_i,".RDS",sep="")
sim_model_i = readRDS(sim_file_name) 
sim_model_i$Count_log_obs = sim_model_i$Count_log_obs[sim_model_i$Obs_year]

species_year_ranges = sim_model_i %>% group_by(Common.Name) %>%
  filter(Obs_year == T) %>% 
  mutate(min_year = min(Year),max_year = max(Year)) %>%
  dplyr::select(Common.Name, min_year,max_year) %>% distinct()

sim_model_i_trim = sim_model_i %>% left_join(species_year_ranges) %>%
  filter(Year >= min_year & Year <= max_year)
  
p = sim_model_i_trim %>% ggplot() + 
  geom_line(aes(y=Count_log,x=Year,group=ID),alpha=0.5) + 
  geom_point(aes(y=Count_log_obs,x=Year,group=ID),alpha=0.4,color="red") + 
  facet_wrap(~Common.Name, scales = "free",ncol=4) + theme_bw() + 
  labs(y= "Log Abundnace")
p
# pdf("model/plots/sim_trends_obs_example.pdf",width = 8, height = 10)
# print(p)
# dev.off()

  
#=== get true index 

# this has to be the real dataset
# Indiv.sum_table = turt_df_long %>% dplyr::select(ID,Indiv.sum) %>% distinct()

true_index = all_tru_sim_df %>% 
  # left_join(Indiv.sum_table) %>%
  # filter(Indiv.sum>0)%>%
  group_by(Common.Name, Year, sim, true_state) %>%
  summarise(total_n = sum(exp(Count_log))) %>%
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                                true_state == "GP" ~ "Synchrony Scenario"))
  

true_index %>% filter(Common.Name == "California sea lion") %>% ggplot() + 
  geom_line(aes(x = Year, y = total_n, group = as.factor(sim), color = true_state)) + 
  facet_wrap(~Common.Name+true_state)

# get coverage
# get the distances from nearest year of observation year
all_tru_sim_df_dist = all_tru_sim_df %>% mutate(distance = NA)  %>% filter(sim == 1 & true_state == "Ueq")
all_tru_sim_df_dist$Count_log_obs = all_tru_sim_df_dist$Count_log_obs[all_tru_sim_df_dist$Obs_year]
all_tru_sim_df_distance = NULL
ts_vec_loop = unique(all_tru_sim_df_dist$ID)
for (sim_i in 1){
  for (i in 1:length(ts_vec_loop)){
    ts_tmp = ts_vec_loop[i]
    ts_fits_tmp = all_tru_sim_df_dist %>% 
      filter(ID == ts_tmp, sim  == sim_i)
    a <- which(!is.na(ts_fits_tmp$Count_log_obs)) 
    b <- which(is.na(ts_fits_tmp$Count_log_obs))
    for (i in 1:length(b)){
      dist_tmp = min(abs(a - b[i]),na.rm=T)
      ts_fits_tmp$distance[b[i]] = dist_tmp}
    all_tru_sim_df_distance = rbind(all_tru_sim_df_distance,ts_fits_tmp)
  }}

all_tru_sim_df_distance$distance[is.na(all_tru_sim_df_distance$distance)] = 0

distance_tbl_to_join = all_tru_sim_df_distance %>% dplyr::select(ID,Year,distance)

# turt_sim_true_obs_coverage = all_tru_sim_df %>% head(2000000) %>%
#   mutate(Obs = case_when(Obs_year == T ~ exp(Count_log),
#                          .default = 0)) %>%
#   group_by(Common.Name, Year, sim, true_state) %>%
#   summarise(coverage = sum((Obs))/sum(exp(Count_log))) %>%
#   mutate(species = Common.Name, year = Year)

turt_sim_true_obs_coverage = all_tru_sim_df %>% 
  left_join(distance_tbl_to_join,by=c("Year","ID")) %>%
  dplyr::select(ID,Year,Common.Name,sim,true_state,Count_log,Obs_year,distance) %>%
  mutate(Obs = case_when(Obs_year == T ~ exp(Count_log),
                         .default = NA)) %>%
  group_by(ID, sim) %>%
  fill(Obs, .direction="downup") %>%
  group_by(Common.Name, Year, sim, true_state) %>%
  summarise(coverage = sum((Obs/(distance+1)))/sum(exp(Count_log))) %>%
  mutate(species = Common.Name, year = Year) %>%
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                                true_state == "GP" ~ "Synchrony Scenario"))

turt_sim_true_obs_n_pts = all_tru_sim_df %>% 
  dplyr::select(-c("ts_id","Subspecies","Site","Indiv.pop.map")) %>%
  filter(Obs_year == T ) %>% 
  group_by(Common.Name, sim, true_state,ID) %>%
  mutate(consec = case_when(lead(Year)-Year == 1 ~ 1,
                            .default = 0)) %>%
  group_by(Common.Name, sim, true_state) %>%
  summarise(n_data_pts=n(),
            n_ts = n_distinct(ID),
            n_consec_y = sum(consec)) %>%
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                                true_state == "GP" ~ "Synchrony Scenario"))

turt_sim_true_obs_n_pairs = all_tru_sim_df %>% 
  dplyr::select(-c("ts_id","Subspecies","Site","Indiv.pop.map")) %>%
  filter(Obs_year == T ) %>%
  group_by(Common.Name, sim, true_state,ID) %>%
  mutate(consec = case_when(lead(Year)-Year == 1 ~ 1,
                            .default = 0)) %>%
  filter(consec>0) %>% ungroup()%>%
  group_by(Common.Name, sim, true_state,Year) %>%
  summarise(n_pairs_pre=n(),consec_sum = sum(consec)) %>% #consec_sum is just a manual check you got calc right
  filter(n_pairs_pre>1) %>%
  mutate(n_pairs = n_pairs_pre*(n_pairs_pre-1)/2) %>%
  group_by(Common.Name, sim, true_state) %>%
  summarise(n_n_pairs = sum(n_pairs)) %>%
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                                true_state == "GP" ~ "Synchrony Scenario"))

# get variance of distance matrix
distance_matrix_list = readRDS("data/site_distance/seal_distance_matrix_subsp_2025_1025.RDS")
species_names = names(distance_matrix_list)
distance_sd = NULL
for (i in 1:length(distance_matrix_list)){
  species_i = names(distance_matrix_list)[i]
  distance_matrix_i = distance_matrix_list[[i]]
  vector_distance_i = distance_matrix_i[distance_matrix_i>0]
  species_sd_i = sd(vector_distance_i) + 1
  if(is.na(species_sd_i)){species_sd_i=0}
  distance_sd_tmp = data.frame(Common.Name = species_i,distance_sd = species_sd_i)
  distance_sd = rbind(distance_sd,distance_sd_tmp)
}

species_ts_attributes = turt_sim_true_obs_n_pts %>% left_join(turt_sim_true_obs_n_pairs) %>% 
  left_join(distance_sd) %>% ungroup() %>%
  mutate(gp_metric = n_consec_y*distance_sd) %>%
  dplyr::select(-c("true_state","sim")) %>% distinct() # metrics are for species level not true_state or sim
species_ts_attributes[is.na(species_ts_attributes)] <- 0

# SIMULATION OUTPUTS
model_vec = c("Ueq","GP")
N_sims = 10
pred_full = NULL
error_full = NULL
grid_full = NULL
# tag = "2x2_2024-09-14"
tag = "GPcorr_prior_2025-10-25"
for (model_i in model_vec){
for (sim_i in 1:N_sims){
  list_files = list.files("data/simulation/sim_outputs")
  list_files = paste("data/simulation/sim_outputs/",list_files,sep="")
  pred_name_tmp = paste("data/simulation/sim_outputs/seal_SIM",model_i,sim_i,"_pred_data_",tag,".RDS",sep="")
  error_name_tmp = paste("data/simulation/sim_outputs/seal_SIM",model_i,sim_i,"_error_params_",tag,".RDS",sep="")
  grid_name_tmp = paste("data/simulation/sim_outputs/seal_SIM",model_i,sim_i,"_mod_sel_grid_",tag,".RDS",sep="")

  if (!(pred_name_tmp %in% list_files)) {next}
  
  pred_tmp = readRDS(pred_name_tmp) %>% mutate(sim = sim_i, true_state = model_i)
  error_tmp = readRDS(error_name_tmp) %>% mutate(sim = sim_i, true_state = model_i)  %>% 
    mutate(species = paste(species,"- redo")) %>% dplyr::select(-q50)
  grid_tmp = readRDS(grid_name_tmp) %>% mutate(sim = sim_i, true_state = model_i)
  
  pred_full =  rbind(pred_full,pred_tmp)
  error_full =  rbind(error_full,error_tmp)
  grid_full =  rbind(grid_full,grid_tmp)
}}


# GRID SELECTION TABLE
library(loo)
loo_tbl = NULL
for (model_i in model_vec){
  for (sim_i in 1:10){
    print(sim_i)
    # list_files = list.files("data/simulation_GP_MA1/sim_outputs")
  loo_tmp = paste("data/simulation/sim_outputs/seal_SIM",model_i,sim_i,"_loo_",tag,".RDS",sep="")
  # if (!(loo_tmp %in% list_files)) {next}

loo_all = readRDS(loo_tmp)
for (i in seq(1,length(loo_all),by=2)){
  speciescut  =  names(loo_all)[i]
  species_i = sub("\\-.*", "", speciescut)

  species_i = str_replace(speciescut, "-GP", "")
  species_i = str_replace(species_i, "-Ueq", "")
  
  loo_compare_tbl = loo_compare(loo_all[[i]],loo_all[[i+1]])
  loo_compare_tbl_sort = loo_compare_tbl[order((row.names(loo_compare_tbl))),]
  
  loo_upper = loo_compare_tbl[2,"elpd_diff"] + 2*loo_compare_tbl[2,"se_diff"] 
  if (loo_upper < -4) {
    if(row.names(loo_compare_tbl)[1] == "model2"){best_model_i =  "GP"}
    if(row.names(loo_compare_tbl)[1] == "model1"){best_model_i =  "Ueq"}
  } else {best_model_i = c("both GP and Ueq")}
  
  loo_tbl_i = data.frame(species = species_i,model_type = best_model_i,sim=sim_i,true_state = model_i)
  loo_tbl = rbind(loo_tbl, loo_tbl_i)
}
}}
loo_tbl = loo_tbl %>% mutate(keep = 1) 
grid_full_loo = grid_full %>% left_join(loo_tbl)


## MODEL SELECTION PROCESS
best_models = loo_tbl %>% arrange(true_state,sim,species) %>%
  # filter(model_type %in% c("Ueq","GP")) %>% 
  # filter(true_state %in% c("Ueq","GP")) %>% 
  group_by(true_state,sim,species) %>% 
  # filter(div_p < 0.01 & rhat1.1 < 0.05) %>%
  # filter(div_p < 0.001) %>%
  # mutate(delta_loo = loo - min(loo)) %>%
  # filter(delta_loo < 4) %>%
  filter(keep == 1) %>%
  left_join(species_ts_attributes,by = c("species"="Common.Name")) %>%
  mutate(label=paste(species,", ","*P*=",n_ts,sep=""))

cols <- c("Synchrony Model" = "#186689", "Both Independent & Synchrony Models" = "#b6aa0d", 
          "Independent Model" = "#e23b0d")

p = best_models %>% 
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                                true_state == "GP" ~ "Synchrony Scenario")) %>%
  mutate(true_state = factor(true_state, levels = c("Synchrony Scenario","Independent Scenario"))) %>%
  mutate(model_type = case_when(model_type == "both GP and Ueq" ~ "Both Independent & Synchrony Models",
                                model_type == "GP" ~ "Synchrony Model")) %>%
  ggplot(aes(x=true_state, fill = model_type)) +
  geom_bar(position = "fill") + coord_flip()+
  scale_fill_manual(name = "Selected Model(s)",values = cols)+
  labs(y = "Proportion",x="True Scenario") +
  facet_wrap(ncol=4,facets = ~reorder(label, n_ts)) +
  # theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = ggfacet::element_textbox_highlight(
          size = 8,
          fill = "white", box.color = "white", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)),
        panel.border = element_rect(colour = "gray80", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "gray10"),
        axis.text.y = element_text(color = "gray10"),
        axis.title.x = element_text(color = "gray10"),
        axis.title.y = element_text(color = "gray10"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "none",
        plot.title = element_text(hjust = .01, vjust=-1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = c(.8,.01), legend.direction = "vertical") 
  
pdf("model/plots/sim_model_selection.pdf",width = 8, height = 8)
print(p)
dev.off()


# ERROR 
error_full = error_full %>%
  mutate(true_state = factor(true_state, levels = c("Ueq","MA1","GP","GP_MA1"))) %>%
  mutate(model_type = factor(model_type, levels = c("Ueq","MA1","GP","GP_MA1"))) %>%
mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                              true_state == "GP" ~ "Synchrony Scenario")) %>%
  mutate(model_type = case_when(model_type == "Ueq" ~ "Independent Model",
                                model_type == "GP" ~ "Synchrony Model"))

#obs error
error_full %>%
  filter(units == "sigma_obs") %>% 
  mutate(mean = case_when(model_type %in% c("Ueq","MA1") ~ mean, 
                          .default = mean)) %>%
  # only include models that were filtered through divergence and looic 
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state), 
  #            by = c("species","model_type","sim","true_state")) %>%
  ggplot() + 
  geom_histogram(aes(x=mean,y=after_stat(ncount))) + geom_vline(xintercept=0.3)+
  facet_wrap(~true_state+model_type,ncol=2)

 p = error_full %>%
  filter(units %in% c("sigma_obs")) %>% 
  mutate(mean = case_when(model_type %in% c("Ueq","MA1") ~ mean, 
                          .default = mean)) %>%
  mutate(species = str_replace(species, " - redo", "")) %>%
  # only include models that were filtered through divergence and looic shit
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state),
  #            by = c("sim"="sim","species"="species","true_state"="true_state","model_type"="model_type")) %>%
  left_join(species_ts_attributes,
             by = c("species"="Common.Name")) %>%
  mutate(n_data_pts = case_when(n_data_pts > -1 ~ n_data_pts,
                                .default = 7800)) %>%
  mutate(n_consec_y = case_when(n_consec_y > -1 ~ n_consec_y,
                                .default = 7822)) %>%
  ggplot() + 
  geom_linerange(aes(x=(n_data_pts),ymin=q025,ymax=q975),alpha=0.5) + 
  geom_point(aes(x=(n_data_pts),y=mean),size=1.5,alpha=0.5) + 
  geom_hline(yintercept=0.30,color="red")+
  scale_x_log10()+ #coord_cartesian(y=c(0,1))+
  facet_wrap(~true_state+model_type,ncol=2) + 
  labs(y = expression("Observation error standard dev"~(sigma[k])),
       x="# of data points") +
  theme_minimal() +
   theme(strip.background = element_blank(),
         strip.text = ggfacet::element_textbox_highlight(
           size = 10,
           fill = "white", box.color = "white", color = "gray10",
           halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
           padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)),
         panel.border = element_rect(colour = "gray80", fill=NA, size=1))

pdf("model/plots/sim_error_sigma_obs.pdf",width = 6, height = 6)
print(p)
dev.off()

error_full %>%
  filter(units %in% c("sigma_process_real")) %>% 
  mutate(hline = case_when(model_type %in% c("GP","Ueq") ~ 0.11,
                          .default = 0.11)) %>%
  # only include models that were filtered through divergence and looic shit
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state), 
  #            by = c("species","model_type","sim","true_state")) %>%
  ggplot() + 
  geom_histogram(aes(x=mean,y=after_stat(ncount))) + geom_vline(xintercept=0.1)+
  facet_wrap(~true_state+model_type,ncol=4)

p=error_full %>%
  filter(units %in% c("sigma_process_real")) %>% 
  mutate(hline = case_when(model_type %in% c("GP","Ueq") ~ 0.11,
                           .default = 0.11)) %>%
  mutate(species = str_replace(species, " - redo", "")) %>%
  # only include models that were filtered through divergence and looic shit
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state), 
  #            by = c("sim"="sim","species"="species","true_state"="true_state","model_type"="model_type")) %>%
  left_join(species_ts_attributes,
            by = c("species"="Common.Name")) %>%
  mutate(n_data_pts = case_when(n_data_pts > -1 ~ n_data_pts,
                                .default = 7800)) %>%
  mutate(n_consec_y = case_when(n_consec_y > -1 ~ n_consec_y,
                                .default = 1669)) %>%
  ggplot() + 
  geom_linerange(aes(x=(n_consec_y),ymin=q025,ymax=q975),alpha=0.5) + 
  geom_point(aes(x=(n_consec_y),y=mean),size=1.5,alpha=0.5) + 
  geom_hline(yintercept=0.11,color="red")+
  scale_x_log10()+
  # coord_cartesian(y=c(0,0.25))+
  facet_wrap(~true_state+model_type,ncol=2)  + 
  labs(y = expression("Process error standard dev"~(omega[k])),
       x="# of consecutive years of data") + 
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = ggfacet::element_textbox_highlight(
          size = 10,
          fill = "white", box.color = "white", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)),
        panel.border = element_rect(colour = "gray80", fill=NA, size=1))

pdf("model/plots/sim_error_sigma_proc.pdf",width = 6, height = 6)
print(p)
dev.off()

error_full %>%
  filter(units %in% c("alpha_gp")) %>% 
  # only include models that were filtered through divergence and looic shit
  inner_join(best_models_tbl%>% select(species,model_type,sim,true_state), 
             by = c("species","model_type","sim","true_state")) %>%
  ggplot() + 
  geom_histogram(aes(x=mean,y=after_stat(ncount))) + geom_vline(xintercept=0.1)+
  facet_wrap(~true_state+model_type,ncol=4)

p = error_full %>%
  filter(units %in% c("alpha_real")) %>% 
  mutate(hline = case_when(true_state %in% c("Synchrony Scenario") ~ 0.9, 
                           .default = 0)) %>%
  mutate(species = str_replace(species, " - redo", "")) %>%
  # only include models that were filtered through divergence and looic shit
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state), 
  #            by = c("sim"="sim","species"="species","true_state"="true_state","model_type"="model_type")) %>%
  left_join(species_ts_attributes,
            by = c("species"="Common.Name")) %>%
  ggplot() + 
  geom_linerange(aes(x=(n_ts),ymin=q025,ymax=q975),alpha=0.5) + 
  geom_point(aes(x=n_ts,y=mean),size=1.5,alpha=0.5) + 
  geom_hline(aes(yintercept=hline),color="red")+
  # scale_x_log10()+
  facet_wrap(~true_state+model_type,ncol=2) + 
  labs(y = expression("Spatial synchrony intercept "~(alpha[k])),
       x="# of populations in species") + 
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = ggfacet::element_textbox_highlight(
          size = 10,
          fill = "white", box.color = "white", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)),
        panel.border = element_rect(colour = "gray80", fill=NA, size=1))


pdf("model/plots/sim_error_sigma_alpha.pdf",width = 6, height = 3)
print(p)
dev.off()

p =error_full %>%
  filter(units %in% c("rho_gp")) %>% 
  mutate(hline = case_when(true_state %in% c("Synchrony Scenario") ~ 50000, 
                           .default = 0)) %>%
  mutate(species = str_replace(species, " - redo", "")) %>%
  # only include models that were filtered through divergence and looic shit
  # inner_join(best_models%>% dplyr::select(species,model_type,sim,true_state), 
  #            by = c("sim"="sim","species"="species","true_state"="true_state","model_type"="model_type")) %>%
  left_join(species_ts_attributes,
            by = c("species"="Common.Name")) %>% 
  ggplot() + 
  geom_linerange(aes(x=(n_ts),ymin=q025,ymax=q975),alpha=0.2) + 
  geom_point(aes(x=n_ts,y=mean),size=1.5,alpha=0.2) + 
  geom_hline(aes(yintercept=hline),color="red")+
  scale_x_log10()+
  facet_wrap(~true_state+model_type,ncol=2) + labs(y = "Synchrony decay (rho)",x="# of populations in species")+ 
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = ggfacet::element_textbox_highlight(
          size = 10,
          fill = "white", box.color = "white", color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)),
        panel.border = element_rect(colour = "gray80", fill=NA, size=1))

pdf("model/plots/sim_error_synchrony_decayrate.pdf",width = 6, height = 3)
print(p)
dev.off()


# ESTIMATED TREND
pred_full = pred_full %>%
  mutate(true_state = case_when(true_state == "Ueq" ~ "Independent Scenario",
                              true_state == "GP" ~ "Synchrony Scenario")) %>%
  mutate(model_type = case_when(model_type == "Ueq" ~ "Independent Model",
                                model_type == "GP" ~ "Synchrony Model"))

est_index = pred_full %>% 
  mutate(ID = as.numeric(word(ts, 2, sep = "_"))) %>% 
  # left_join(Indiv.sum_table) %>%
  # filter(Indiv.sum>0)%>%
  group_by(species,sim,model_type,true_state,year) %>%
  summarise(est_total_n = mean(tot.q50),
            est_total_n_lo = mean(tot.q10),
            est_total_n_hi = mean(tot.q90),
            sum_data = sum(exp(data),na.rm=T)) 

p = pred_full %>%
  filter((species == "Northern elephant seal" & sim == 2 & true_state == "Synchrony Scenario")) %>%
  ggplot() + 
  geom_ribbon(aes(x=year,ymin = mean-1.96*sd, ymax = mean+1.96*sd,fill=model_type), alpha=0.3) +
  # geom_line(aes(x=year,y=state_fit)) +
  geom_line(aes(x=year,y=mean,color=model_type)) + 
  geom_point(aes(x=year,y=data),color="grey10",alpha=0.3)+
  facet_wrap(~ts,scales="free",ncol=4) + 
  labs(y="Log abundance", x= "Year")+
  scale_color_manual(name = "Selected Model(s)",values = cols)+
  scale_fill_manual(name = "Selected Model(s)",values = cols)+
  theme_classic() +
  theme(legend.position="bottom")
print(p)

# pdf("model/plots/sim_model_type_example_NES.pdf",width = 8, height = 8)
# print(p)
# dev.off()

# joing est trend to true tend AND MAPE and coverage

est_n_true_index = est_index %>% 
  left_join(true_index, by = c("sim"="sim","species"="Common.Name","year"="Year","true_state"="true_state")) %>%
  mutate(prop_obs = sum_data/total_n)

p = est_n_true_index %>% 
  filter(sim==10, true_state=="Synchrony Scenario") %>%
  filter(sim==10, true_state=="Synchrony Scenario") %>%
  ungroup() %>% rowwise() %>%
  mutate(prop_obs = min(prop_obs,1)) %>%
  ggplot() + 
  # geom_vline(aes(xintercept=year, alpha = prop_obs),color="forestgreen") +
  geom_line(aes(x=year,y=total_n,color=true_state),color="black",alpha=0.5) + 
  geom_point(aes(x=year,y=total_n, size = prop_obs),color="black",alpha=0.5) +

  geom_ribbon(aes(x=year,ymax=est_total_n_hi,ymin=est_total_n_lo,fill=model_type),alpha=0.3) + 
  geom_line(aes(x=year,y=est_total_n,color=model_type,linetype= model_type),lwd=1.2) + 
  scale_size(range=c(0,1.5),breaks=c(0,.25,.50,.75,1.0),guide="legend")+
  expand_limits(y = c(0)) + 
  facet_wrap(~species,scales = "free",ncol=4) + 
  labs(y="Total Species Abundance",color="Model Type",fill="Model Type",linetype = "Model Type", size = "% observed") + 
  scale_color_manual(values = cols)+scale_fill_manual(values = cols)+
  theme_classic() + 
  theme(legend.position = c(.8,.03), legend.box = "horizontal") +
  guides(color=guide_legend(ncol=1), size=guide_legend(ncol=2),legend.spacing.y = unit(0.2, "cm"))

p
pdf("model/plots/sim_model_trend_example.pdf",width = 8.5, height = 8)
print(p)
dev.off()

# Mean Absolute Percent Error from EStimated index to true index
sim_MAPE_coverage = est_n_true_index %>%
  filter(!(species == "New Zealand sea lion")) %>%
  group_by(species,sim,model_type,true_state,year) %>%
  left_join(turt_sim_true_obs_coverage) %>%
  mutate(APE = 100*abs(total_n-est_total_n)/total_n, coverage = coverage*100) %>% 
  filter(APE < 100) %>%
  filter(coverage < 100) %>%
  mutate(model_type=as.factor(model_type)) 

library(lme4)

sim_MAPE_coverage$coverage_scaled <- scale(sim_MAPE_coverage$coverage)

final_model <- lmer(
  data = sim_MAPE_coverage, 
  log(APE) ~ coverage_scaled * true_state * model_type + 
    (1 + coverage_scaled | sim:species),
  control = lmerControl(optimizer = "bobyqa"))

summary(final_model)

plot_dat <- expand.grid(
  coverage = seq(min(sim_MAPE_coverage$coverage), 
                 max(sim_MAPE_coverage$coverage), 
                 length.out = 100),
  true_state = unique(sim_MAPE_coverage$true_state),
  model_type = unique(sim_MAPE_coverage$model_type))

cov_mean <- mean(sim_MAPE_coverage$coverage)
cov_sd  <- sd(sim_MAPE_coverage$coverage)

plot_dat$coverage_scaled <- (plot_dat$coverage - cov_mean) / cov_sd

plot_dat$pred_log_APE <- predict(final_model, newdata = plot_dat, re.form = NA)

plot_dat$pred_APE <- exp(plot_dat$pred_log_APE)

p=ggplot() +
  geom_point(data = sim_MAPE_coverage, 
             aes(x = coverage, y = APE, color = model_type), 
             alpha = 0.15, size = 1) +
  geom_line(data = plot_dat, 
            aes(x = coverage, y = pred_APE, color = model_type), 
            linewidth = 1.2) +
  facet_wrap(~true_state) +
  coord_cartesian(ylim=c(0, 20)) +
  scale_color_manual(values = c("Independent Model" = "#e23b0d", "Synchrony Model" = "#186689"), 
                     name = "Model Type",
                     labels = c("Independent Model", "Synchrony Model")) +
  theme_bw() +
  labs(y = "Absolute % Error", 
       x = "% Coverage") +
  theme(legend.position = "bottom")
pdf("model/plots/sim_model_type_APE_regression.pdf",width = 8, height = 4)
print(p)
dev.off()

library(emmeans)
library(dplyr)

# overall mean APE
emmeans::emmeans(final_model, ~ 1, type = "response")

# overall mean APE at coverage = 100%
emmeans::emmeans(final_model, ~ 1, 
                          at = list(coverage_scaled = full_cov_scaled), 
                          type = "response")


# difference btween synchrony and independent model in syncrhony secnario
library(emmeans)
emmeans(final_model, ~ model_type | true_state, 
                           at = list(coverage_scaled = avg_cov_scaled),
                           type = "response")

### Evaluate how well 80% posterior intervals captures TRUE ABUNDANCE
coverage_analysis = est_n_true_index %>% 
  mutate(caught = total_n >= est_total_n_lo & total_n <= est_total_n_hi) 

coverage_summary = coverage_analysis %>%
  group_by(true_state, model_type) %>%
  summarise(
    prop_covered = mean(caught, na.rm = TRUE),
    n_years = n())


p = ggplot(coverage_summary, aes(x = true_state, y = prop_covered, fill = model_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "True Abundance in 80% Posterior Interval", x = "True Scenario",) +
  scale_fill_manual(values = c("Independent Model" = "#e23b0d", "Synchrony Model" = "#186689"), 
                     name = "Model Type",
                     labels = c("Independent Model", "Synchrony Model")) +
  theme_bw()

pdf("model/plots/sim_model_type_80predictive_interval.pdf",width = 8, height = 4)
print(p)
dev.off()
