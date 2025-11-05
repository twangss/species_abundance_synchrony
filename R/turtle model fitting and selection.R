library(tidyverse)
library(ggplot2)
library(ggforce)
library(atsar)
library(rstan)
library(loo)
library(magrittr)
library(atsar)

source("R/marss_stan_functions.R")
turt_df = read.csv("data/All Turtles Data_latlon.csv") 

turt_df_long = turt_df %>% dplyr::select(-Count) %>% 
  mutate(ts_id = paste(CommonName,"#",Site)) %>% 
  filter(Use == "yes") %>%
  group_by(ID) %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Year",
    names_prefix = "X",
    values_to = "Count",
    values_drop_na = TRUE) %>%
  mutate(Count=as.numeric(Count),Year=as.numeric(Year))  %>%
  group_by(ID) %>%arrange(ID) %>%
  drop_na(Count) %>%
  mutate(n=n()) %>% filter(n>2) #%>% mutate(Scaling = 1)

turt_df_long %>% ungroup() %>% select(CommonName,RMU) %>%
  distinct()

# Prep model fitting loop
species_names <- unique(turt_df_long$CommonName)
N_species <- length(species_names)

distance_matrix_list = readRDS("data/site_distance/turtle_distance_matrix_2025_0925.RDS")

mod_sel_grid=NULL
error_params=NULL
pred_data=NULL
loo_all = NULL
gaussian_rel = NULL

taxa = "turtle_real"
start_date = Sys.Date()
mod_sel_file_name = paste("data/outputs/",taxa,"_mod_sel_grid_",start_date,".RDS",sep="")
pred_data_file_name = paste("data/outputs/",taxa,"_pred_data_",start_date,".RDS",sep="")
error_params_file_name = paste("data/outputs/",taxa,"_error_params_",start_date,".RDS",sep="")
loo_file_name = paste("data/outputs/",taxa,"_loo_",start_date,".RDS",sep="")
gaussian_file_name = paste("data/outputs/",taxa,"_gp_",start_date,".RDS",sep="")

options(mc.cores = parallel::detectCores())

for (i in 1:N_species){
  species_i <- species_names[i]
  cat(species_i,"\n")
  species_i_pre_df <- turt_df_long %>% 
    filter(CommonName == species_i) %>% arrange(ID) %>%
    mutate(log.spawner = log(Count+1)) %>% 
    group_by(ID) %>% 
    mutate(n=sum(!is.na(log.spawner))) 
  
  species_i_df <- species_i_pre_df %>% ungroup() %>% 
    dplyr::select(ID,Year,Site, log.spawner) %>%  
    # drop_na(log.spawner) %>%
    group_by(ID,Site) %>% 
    complete(Year = min(Year):max(Year),
             fill = list(log.spawner = NA)) %>%
    ungroup() %>%
    pivot_wider(names_from = c("Site","ID"), values_from = "log.spawner") %>% 
    arrange(Year) %>% 
    column_to_rownames(var = "Year") %>% 
    as.matrix() %>% # turn into a matrix with year down the rows
    t() # make time across the columns
  
  # set up state hypotheses
  id_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>%
    distinct(ID,Indiv.pop.map) %>% pull(Indiv.pop.map)
  RMU_vec = species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>% filter(Indiv.sum==1) %>% distinct(ID,RMU) %>%
    pull(RMU)
  species_i_pre_df = species_i_pre_df %>% mutate(Sampling.method = Units_clean)
  add_up_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>%
    distinct(ID,Indiv.sum,Scaling) %>% mutate(add_up = Indiv.sum*Scaling) %>% pull(add_up)
  
  #== Set up for Stan
  id_states = as.numeric(factor(id_vec))
  stocks <- (as.numeric(factor(RMU_vec))) 
  n_states <- max(as.numeric(factor(id_vec)))
  r_unequal_units = (as.numeric(fct_inorder(units_vec)))
  unique_units_vec = unique(units_vec)
  
  # covariance matrices
  r_unequal <- seq(1, nrow(species_i_df))
  r_equal <- rep(1, nrow(species_i_df))
  uq_unequal <- seq(1, n_states)
  uq_equal <- rep(1, n_states)
  

  # set up mcmc options
    mcmc_list  = list(n_mcmc = 4000, n_burn = 1000, n_chain = 3, n_thin = 1,step_size=0.4,adapt_delta=0.9)
      
  # obs error specifications
    method_vec_in = r_equal
    method_vec = "sigma_obs"
      
  # MARSS specifications
    marss_list = list(states = id_states, 
                        obsVariances = method_vec_in, 
                        proVariances = uq_equal,
                        trends = uq_equal,
                        stocks = uq_unequal) 
      
  # set up data_list
    data_list_tmp <- setup_data(y = species_i_df,
                                  est_nu = FALSE,
                                  est_trend = FALSE,
                                  family = "gaussian", 
                                  marss = marss_list)
      
      
  # Spatial covariance matrix
    data_list_tmp$data$Dist = distance_matrix_list[[species_i]]
    data_list_tmp$data$add_up = add_up_vec
      
    if (nrow(species_i_df) == 1) {data_list_tmp$data$Dist = matrix(0,1,1) } # if only 1 population

    model_vec = c("Ueq","GP") # candidate models: Ueq - independent processes, GP - spatial synchrony process

    for (model_i in model_vec){
      cat("Model",model_i,"\n")
      model_type = model_i
      
      if (model_type == "Ueq") {cmd_file_name = "marss_cmd_prior.stan"}
      if (model_type == "GP") {cmd_file_name = "marss_gp_cmd_prior.stan"}
      
      # Priors for process and observation error params
      # Inverse gamma priors on sigmas: process error (alpha, beta), observation error (alpha, beta)
      # half normal prior on spatial synchrony decay
      # beta prior on spatial synchrony intercept
      data_list_tmp$data$sigma_obs_priors = c(0.05, 0.2)
      data_list_tmp$data$sigma_proc_priors = c(0.05, 0.05)
      data_list_tmp$data$rho_gp_priors = c(0, 80000)
      data_list_tmp$data$alpha_gp_priors = c(1, 2)
        
        
      
      library(cmdstanr)
      file <- file.path(cmdstan_path(), "pinnipeds", cmd_file_name) # user needs to specificy CmdStan folder
      mod <- cmdstan_model(file)
      
      fit2 <- mod$sample(
        data = data_list_tmp$data,
        seed = 124,
        chains = mcmc_list$n_chain,
        parallel_chains = mcmc_list$n_chain,
        iter_warmup = mcmc_list$n_burn,
        iter_sampling = mcmc_list$n_mcmc,
        adapt_delta=0.99,
        step_size=0.4,
        refresh = 500 # print update every 500 iters
      )
      cat("fitting is done \n")

      # save model diags
      fit2_sum = fit2$summary()

      diag_summ = fit2$diagnostic_summary()
      divergent_perc = sum(diag_summ$num_divergent)/(mcmc_list$n_chain*mcmc_list$n_mcmc)

      eff_sample_size = fit2_sum$ess_tail/(mcmc_list$n_chain*mcmc_list$n_mcmc)
      n_eff0.01_tmp = sum(eff_sample_size<0.01,na.rm=T)/length(eff_sample_size)
      
      rhat1.1_tmp = sum(fit2_sum$rhat>1.01,na.rm=T)/length(fit2_sum$rhat)
      
      model_diags = data.frame(species = species_i, #rmu = rmu_i,proc_state = Q, 
                               rhat1.1=rhat1.1_tmp,n_eff0.01=n_eff0.01_tmp,div_p=divergent_perc,
                               loo = loo_tmp$estimates["looic","Estimate"],
                               loo_se = loo_tmp$estimates["looic","SE"]) %>% 
        mutate(model_type = model_type,prior = prior_i)
      mod_sel_grid = rbind(mod_sel_grid,model_diags) 
      
      # save loos for model selection
      loo_tmp = loo(fit2$draws("log_lik"))
      
      loo_list_name = paste(species_i,model_type,sep="-")
      loo_all[[loo_list_name]] = loo_tmp
      
      saveRDS(loo_all,loo_file_name)

      # Save Process/Obs Sigmas
      if (model_type == "Ueq") {proc_variable_names = c("sigma_process_real")}
      if (model_type == "GP") {proc_variable_names = c("sigma_process_real","alpha_real","rho_gp")}
      
      sigma_proc_draws <- fit2$draws(proc_variable_names,format = "matrix")
      sigma_proc_tmp = sigma_proc_draws %>% data.frame() %>% pivot_longer(everything(), names_to = "units") %>%
        group_by(units) %>%
        summarise(mean = mean(value),
                  sd = sd(value),
                  q025 = quantile(value, probs = c(0.025)),
                  q50 = quantile(value, probs = c(0.5)),
                  q975 = quantile(value, probs = c(0.975)),
                  rhat = fit2$summary(proc_variable_names)$rhat) %>%
        mutate(var="sigma_proc")

      sigma_obs_draws <- fit2$draws("sigma_obs",format = "matrix")
      sigma_obs_tmp = data.frame(units = levels(fct_inorder(method_vec)), #"sigma_obs", #unique_units_vec, #levels(fct_inorder(method_vec))
                                 mean = apply(sigma_obs_draws, 2, mean),
                                 sd = apply(sigma_obs_draws, 2, sd),
                                 q025 = apply(sigma_obs_draws,2,quantile, probs = c(0.025)),
                                 q50 = quantile(sigma_obs_draws, probs = c(0.5)),
                                 q975 = apply(sigma_obs_draws,2,quantile, probs = c(0.975)),
                                 rhat = fit2$summary("sigma_obs")$rhat,
                                 var="sigma_obs")
      
      error_params_tmp = rbind(sigma_proc_tmp,sigma_obs_tmp) %>% mutate(model_type = model_type,prior = prior_i)
      error_params_tmp$species = species_i 
      error_params = rbind(error_params,error_params_tmp)
      
      # Spatial synchrony relationship
      if (model_type == "GP") {
        gaussian_rel_out = fit2$draws("gaussian_rel",format = "matrix")
        gaussian_rel_tmp = data.frame(tot.q10 = apply(gaussian_rel_out,2,quantile, probs = c(0.10)),
                                      tot.q50 = apply(gaussian_rel_out,2,quantile, probs = c(0.50)),
                                      tot.q90 = apply(gaussian_rel_out,2,quantile, probs = c(0.90)),
                                      var="gaussian_rel",species = species_i) %>%
          mutate(x = 1:10000)
        gaussian_rel = rbind(gaussian_rel,gaussian_rel_tmp) 
        saveRDS(gaussian_rel,gaussian_file_name)
      }
      
      # Save out raw data, population-specific predicted abundances, and total abundance estimates
      # Get predicted data and raw data
      raw_data = species_i_df %>% data.frame() %>% 
        set_colnames(colnames(species_i_df)) %>%
        rownames_to_column("ts") %>% 
        pivot_longer(!ts,
                     names_to="year", 
                     values_to = "data") %>% 
        mutate(year=as.numeric(year)) 
      
      # total species abundance 
      tot_sum_x = fit2$draws("tot_sum_x",format = "matrix")
      
      tot_sum = data.frame(tot.mean = apply(tot_sum_x, 2, mean),
                           tot.sd = apply(tot_sum_x, 2, sd),
                           tot.q05 = apply(tot_sum_x,2,quantile, probs = c(0.05)),
                           tot.q10 = apply(tot_sum_x,2,quantile, probs = c(0.10)),
                           tot.q25 = apply(tot_sum_x,2,quantile, probs = c(0.25)),
                           tot.q50 = apply(tot_sum_x,2,quantile, probs = c(0.50)),
                           tot.q75 = apply(tot_sum_x,2,quantile, probs = c(0.75)),
                           tot.q90 = apply(tot_sum_x,2,quantile, probs = c(0.9)),
                           tot.q95 = apply(tot_sum_x,2,quantile, probs = c(0.95)),
                           var="total_abundance_norm") %>%
        mutate(year = as.numeric(colnames(species_i_df)))
      
      # extract individual population abundances
      preds = fit2$summary(variables = "pred")
      preds_matrix = cmdstanr_preds_sd(preds,species_i_df)
      pred_data_tmp = preds_matrix %>% left_join(raw_data) %>% 
        mutate(species = species_i, model_type = model_type) %>%
        left_join(tot_sum)
      
      p = pred_data_tmp %>% ggplot() + 
        geom_ribbon(aes(x=year,ymin = mean-1.96*sd, ymax = mean+1.96*sd), alpha=0.2) +
        geom_line(aes(x=year,y=(mean))) +
        geom_point(aes(x=year,y=(data)),color="red")+
        facet_wrap(~ts,scales="free") + theme_bw()
      print(p)
      
      pred_data=rbind(pred_data,pred_data_tmp)
    
      # save the files
      saveRDS(mod_sel_grid,mod_sel_file_name)
      saveRDS(error_params,error_params_file_name)
      saveRDS(pred_data,pred_data_file_name)
      
      
      
  }
}


## Read in model outputs
mod_sel_file_name = "data/outputs/turtle_real_mod_sel_grid_2025-09-17.RDS"
mod_sel_grid = readRDS(mod_sel_file_name) 

error_params_file_name = "data/outputs/turtle_real_error_params_2025-09-17.RDS"
error_params = readRDS(error_params_file_name) 

gaussian_file_name = "data/outputs/turtle_real_gp_2025-09-17.RDS"
gaussian_rel = readRDS(gaussian_file_name) 

loo_file_name = "data/outputs/turtle_real_loo_2025-09-17.RDS"
loo_all = readRDS(loo_file_name) 

pred_data_file_name = "data/outputs/turtle_real_pred_data_2025-09-17.RDS"
pred_data = readRDS(pred_data_file_name) 

## Model selection using LOO comparisons
loo_tbl_pre = NULL
i=1
for (i in seq(1,length(loo_all),by=2)){ 
  speciescut  =  names(loo_all)[i]
  species_i = sub("\\-.*", "", speciescut)
  
  species_i = str_replace(speciescut, "-GP", "")
  species_i = str_replace(species_i, "-Ueq", "")
  
  loo_compare_tbl = loo_compare(loo_all[[i]],loo_all[[i+1]])
  loo_compare_tbl_sort = loo_compare_tbl[order((row.names(loo_compare_tbl))),]
  
  loo_upper = loo_compare_tbl[2,"elpd_diff"] + 1.282*loo_compare_tbl[2,"se_diff"] 
  loo_lower = loo_compare_tbl[2,"elpd_diff"] - 1.282*loo_compare_tbl[2,"se_diff"] 
  if (loo_upper < -2) {
    if(row.names(loo_compare_tbl)[1] == "model2"){best_model_i =  "GP"}
    if(row.names(loo_compare_tbl)[1] == "model1"){best_model_i =  "Ueq"}
  } else {best_model_i = c("GP","Ueq")}
  
  loo_tbl_i = data.frame(species = species_i,model_type = best_model_i)
  loo_tbl_pre = rbind(loo_tbl_pre, loo_tbl_i)
}

# select the best model per species, prioritizing Ueq over GP if no difference
loo_tbl = loo_tbl_pre %>% 
  mutate(keep = 1)  %>% 
  group_by(species) %>% 
  mutate(priority = case_when(model_type == "GP" ~ 2,
                              model_type == "Ueq" ~ 1,
                              .default = NA)) %>%
  filter(priority == min(priority))

# confirm selected model diags are good
model_grid_best = mod_sel_grid %>% 
  mutate(priority = case_when(model_type == "GP" ~ 2,
                              model_type == "Ueq" ~ 1,
                              .default = NA)) %>%
  left_join(loo_tbl) %>%
  group_by(species) %>% 
  filter(keep == 1) %>%
  filter(priority == min(priority))

best_pass2_model_fits = pred_data %>% 
  inner_join(model_grid_best, by = c("species","model_type")) 

#==== Plot 
# constrain the ballooning errors for plotting purposes, find the distance (# of years) to the nearest data value
best_pass2_model_fits_dist = best_pass2_model_fits %>% 
  group_by(ts) %>%
  mutate(distance = NA)

ts_fits_distance = NULL
ts_vec_loop = unique(best_pass2_model_fits_dist$ts)
for (i in 1:length(ts_vec_loop)){
  ts_tmp = ts_vec_loop[i]
  ts_fits_tmp = best_pass2_model_fits_dist %>% filter(ts == ts_tmp)
  
  a <- which(!is.na(ts_fits_tmp$data)) 
  b <- which(is.na(ts_fits_tmp$data))
  
  for (i in 1:length(b)){
    dist_tmp = min(abs(a - b[i]),na.rm=T)
    ts_fits_tmp$distance[b[i]] = dist_tmp
  }
  
  ts_fits_distance = rbind(ts_fits_distance,ts_fits_tmp)
}

# constrain the ballooning errors past 5 years
best_pass2_model_fits_error_fix = ts_fits_distance %>% ungroup() %>% # 
  mutate(sd_tmp = sd) %>%
  mutate(sd_tmp = case_when(distance > 5 ~ NA,
                             .default = sd_tmp)) %>%
  dplyr::group_by(ts) %>%
  fill(c(sd_tmp), .direction = "downup") #change to up only if necessary


# add unit scaling to model fits for labeling 
best_pass2_model_fits_error_fix_id = best_pass2_model_fits_error_fix %>%
  mutate(ID = as.numeric(stringr::word(ts, 2, sep = "_")))

scaling_tbl = turt_df_long %>% 
  ungroup() %>%
  dplyr::select(CommonName,ID,Scaling,RMU, Units,Indiv.sum) %>% # 
  distinct()

best_pass2_model_fits_error_fix_scale = best_pass2_model_fits_error_fix_id %>% 
  left_join(scaling_tbl) 

best_pass2_model_fits_error_fix_scale_sd = best_pass2_model_fits_error_fix_scale %>% 
  ungroup() %>% 
  # group_by(year) %>%
  mutate(norm_mean = exp((mean)), # just look at mode, if mean then: +0.5*sd_tmp^2
         normal_sd = exp(mean+(sd_tmp^2)/2)*sqrt(exp(sd_tmp^2)-1))

#=== plot the Stan fits for individual population trends
best_pass2_model_fits_error_fix_scale$label <- paste(best_pass2_model_fits_error_fix_scale$ts,"\n","units = ", best_pass2_model_fits_error_fix_scale$Units,sep="")
best_pass2_species = unique(best_pass2_model_fits_error_fix_scale$species)
# pdf("model/plots/turtles_best_species_marss_stan_ts_fits_20250926.pdf", width = 11, height = 8.5)
ncol = 5
nrow=4
for (species_i in best_pass2_species){
  print(species_i)
  species_best_pass2_fits = best_pass2_model_fits_error_fix_scale %>% filter(species == species_i)
  RMU_vec = unique(species_best_pass2_fits$RMU)
  for(RMU_i in RMU_vec){
    species_best_pass2_fits_proc = species_best_pass2_fits %>% 
      filter(RMU == RMU_i)
  page_length = ceiling(length(unique(species_best_pass2_fits_proc$ts))/(ncol*nrow))
  for (page_i in 1:page_length){
    p = species_best_pass2_fits_proc %>% ggplot() +
      geom_ribbon(aes(x=year,ymin = mean-sd_tmp, ymax = mean+sd_tmp), alpha=0.2) +
      geom_line(aes(x=year,y=mean)) +
      geom_point(aes(x=year,y=data),color="red")+
      facet_wrap_paginate(~label,scales="free", ncol = ncol, nrow = nrow, page = page_i) +
      # facet_wrap(~ts,scales="free") + 
      theme_bw() +
      ggtitle(paste(species_i,RMU_i,
                    sep=" - "))
    print(p)
  }}}
dev.off()


#==== INDEX BY SPECIES
library(dplyr)
library(purrr)
library(broom)

pre_species_plot = best_pass2_model_fits_error_fix_scale_sd %>%
  group_by(species) %>%
  mutate(n=n_distinct(ts)) %>%
  mutate(species_label = paste(species,", ","n = ", n,sep="")) %>%
  mutate(species_mean = mean(tot.q50)) %>%
  group_by(ts) %>%
  mutate(ts_mean = median(exp(data),na.rm=T)) 

ymax_tbl <- pre_species_plot %>%
  group_by(species_label) %>%
  summarise(ymax = max(tot.q75*1.2, na.rm = TRUE), .groups = "drop")

smooth_df <- pre_species_plot %>%
  left_join(ymax_tbl, by = "species_label") %>%
  group_by(species_label, ts) %>%
  group_modify(~ {
    fit <- loess(exp(data)/ts_mean*species_mean ~ year, data = .x)
    out <- augment(fit, newdata = .x)  # guarantees same number of rows
    out$smooth <- out$.fitted
    out$smooth_clipped <- pmin(out$smooth, .x$ymax)  # clip to tot.q75
    out
  }) %>%
  ungroup()

species_trend_plot <- pre_species_plot %>%
  dplyr::select(species_label, year, tot.q50, tot.q25, tot.q75) %>%
  distinct() %>%
  left_join(ymax_tbl, by = "species_label") %>%
  ggplot(aes(x = year, y = tot.q50, ymin = tot.q25, ymax = tot.q75)) +
  
  # Use the clipped smooth
  geom_line(data = smooth_df,
            aes(y = smooth_clipped, group = ts),
            color = alpha("steelblue", 0.2), size = 0.5) +
  
  geom_ribbon(alpha = 0.5, fill = "black") +
  geom_line(size = 1.1, color = "black") +
  
  facet_wrap(~species_label, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = function(x) x / 1000,limits = c(0, NA)) +
  # scale_y_continuous()
  theme_classic() +
  labs(y = "Monitored Adult Female Abundance (1000s)", x = "Year")

# ggsave("plots/turtle_species_trends.pdf", species_trend_plot, width=8, height=10)
