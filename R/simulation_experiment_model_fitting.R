library(tidyverse)
library(ggplot2)
library(ggforce)
library(atsar)
library(rstan)
library(loo)
library(magrittr)
library(atsar)

source("R/marss_stan_functions.R")

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

start_date = Sys.Date()
# start_date = "2025-10-25"
N_sims = 10
for (true_state in c("Ueq","GP")){ #"Ueq",
  for (sim_i in 1:10){ #N_sims
    sim_file_name = paste("data/simulation/seal_sim_true_obsV6_",true_state,"_sim",sim_i,".RDS",sep="")
    
    turt_sim_true_obs = readRDS(sim_file_name)
    turt_df_long = turt_sim_true_obs %>%
      filter(Obs_year == T) %>% # test to make sure simulation good
      mutate(Count = exp(Count_log_obs),
             Scaling = 1, Units = "sim", CommonName = Common.Name, RMU = Subspecies)
    
    species_names <- unique(turt_df_long$CommonName)
    N_species <- length(species_names)
    
    # get the distance matrix
    distance_matrix_list = readRDS("data/site_distance/seal_distance_matrix_subsp_2025_1025.RDS")
    
    mod_sel_grid=NULL
    error_params=NULL
    pred_data=NULL
    loo_all = NULL
    gaussian_rel = NULL

    taxa = paste("seal_SIM",true_state,sim_i,sep="")
    model_type_name = "GPcorr_prior"
    mod_sel_file_name = paste("data/simulation/sim_outputs/",taxa,"_mod_sel_grid_",model_type_name,"_",start_date,".RDS",sep="")
    pred_data_file_name = paste("data/simulation/sim_outputs/",taxa,"_pred_data_",model_type_name,"_",start_date,".RDS",sep="")
    error_params_file_name = paste("data/simulation/sim_outputs/",taxa,"_error_params_",model_type_name,"_",start_date,".RDS",sep="")
    loo_file_name = paste("data/simulation/sim_outputs/",
                          taxa,"_loo_",model_type_name,"_",start_date,".RDS",sep="")
    gaussian_file_name = paste("data/simulation/sim_outputs/",
                               taxa,"_gp_",model_type_name,"_",start_date,".RDS",sep="")
    
    options(mc.cores = parallel::detectCores())
    
    for (i in 1:N_species){
      species_i <- species_names[i]
      cat(species_i,"\n")
      species_i_pre_df <- turt_df_long %>% 
        filter(CommonName == species_i) %>% arrange(ID) %>%
        # filter(RMU == rmu_i) %>% # get only this ESU
        mutate(log.spawner = log(Count+1)) %>% # create a column called log.spawner
        group_by(ID) %>% 
        mutate(n=sum(!is.na(log.spawner))) # %>% 
      #filter(n>2) 
      
      species_i_df <- species_i_pre_df %>% ungroup() %>% 
        dplyr::select(ID,Year,Site, log.spawner) %>%  # get just the columns that I need
        # drop_na(log.spawner) %>%
        group_by(ID,Site) %>% 
        complete(Year = min(Year):max(Year),
                 fill = list(log.spawner = NA)) %>%
        ungroup() %>%
        pivot_wider(names_from = c("Site","ID"), values_from = "log.spawner") %>% 
        arrange(Year) %>% 
        column_to_rownames(var = "Year") %>% # make the years rownames
        # filter(if_any(everything(), ~ !is.na(.)))%>%
        as.matrix() %>% # turn into a matrix with year down the rows
        t() # make time across the columns
      
      # Set up data input vectors for states
      # id_vec for seals which has multiple time series mapping onto 1 location
      id_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>%
        filter(Indiv.sum==1) %>% 
        distinct(ID,Indiv.pop.map) %>% pull(Indiv.pop.map)
      # Boolean to add populations or not
      add_up_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>%
        distinct(ID,Indiv.sum,Scaling) %>% mutate(add_up = Indiv.sum*Scaling) %>% pull(add_up)
      # Subspecies vector
      RMU_vec = species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>% filter(Indiv.sum==1) %>% distinct(ID,RMU) %>%
        pull(RMU)
      
      #== Set up for Stan
      # seal
      n_states <- max(as.numeric(fct_inorder(id_vec)))
      id_states = as.numeric(fct_inorder(id_vec))
      
      # all
      r_unequal <- seq(1, nrow(species_i_df))
      r_equal <- rep(1, nrow(species_i_df))
      uq_unequal <- seq(1, n_states)
      uq_equal <- rep(1, n_states)
      stocks <- (as.numeric(factor(RMU_vec))) 
      
      cat("Species",species_i,"\n")
      
      # set up mcmc options
      mcmc_list  = list(n_mcmc = 4000, n_burn = 1000, n_chain = 3, n_thin = 1,step_size=0.4,adapt_delta=0.9)
      
      method_vec_in = r_equal
      method_vec = "sigma_obs"
      
      #= MARSS specifications
      marss_list = list(states = id_states, #site_vec
                        obsVariances = method_vec_in, #as.numeric(factor(units_vec)), OR r_equal
                        proVariances = uq_equal,
                        trends = uq_equal,
                        stocks = uq_unequal) # stocks = stocks, this was 
      
      # set up data_list
      data_list_tmp <- setup_data(y = species_i_df,
                                  est_nu = FALSE,
                                  est_trend = FALSE,
                                  family = "gaussian", 
                                  # mcmc_list = mcmc_list,
                                  marss = marss_list)
      
      #== subspecies adding
      data_list_tmp$data$n_RMU = length(unique(RMU_vec))
      data_list_tmp$data$RMU = (as.numeric(factor(RMU_vec)))
      
      #== Spatial covariance matrix
      data_list_tmp$data$Dist = distance_matrix_list[[species_i]]
      data_list_tmp$data$add_up = add_up_vec
      
      if (nrow(species_i_df) == 1) {data_list_tmp$data$Dist = matrix(0,1,1) }
      
      model_vec = c("Ueq","GP")
      
      for (model_i in model_vec){
        cat("true_state",true_state,"\n")
        cat("Model",model_i,"\n")
        model_type = model_i
        
        if (model_type == "Ueq") {cmd_file_name = "marss_cmd_prior.stan"}
        if (model_type == "GP") {cmd_file_name = "marss_gp_cmd_prior.stan"}

        
        # Priors for process and observation error params
        # Inverse gamma priors on sigmas: process error (alpha, beta), observation error (alpha, beta)
        # Multiple priors for sensitivity, can be modified as needed
        param1 <- list(c(0.05, 0.2), c(0.1, 0.1))
        param2 <- list(c(0.05, 0.05), c(0.01, 0.01))
        param3 <- list(c(0, 80000), c(0, 160000))
        param4 <- list(c(1, 2), c(1, 1))
        
        # Generate all combinations of the paired values
        prior_combo <- expand.grid(param1, param2, param3, param4)
        
        # Flatten the list columns into separate columns
        prior_combo <- data.frame(
          Param1_1 = sapply(prior_combo$Var1, "[", 1),
          Param1_2 = sapply(prior_combo$Var1, "[", 2),
          Param2_1 = sapply(prior_combo$Var2, "[", 1),
          Param2_2 = sapply(prior_combo$Var2, "[", 2),
          Param3_1 = sapply(prior_combo$Var3, "[", 1),
          Param3_2 = sapply(prior_combo$Var3, "[", 2),
          Param4_1 = sapply(prior_combo$Var4, "[", 1),
          Param4_2 = sapply(prior_combo$Var4, "[", 2)
        )
        
        # Only use the first set of priors for simulation experiment
        prior_combo = prior_combo[1,]
        
        for (prior_i in 1:nrow(prior_combo)){
          data_list_tmp$data$sigma_obs_priors = c(prior_combo[prior_i,1],prior_combo[prior_i,2])
          data_list_tmp$data$sigma_proc_priors = c(prior_combo[prior_i,3],prior_combo[prior_i,4])
          data_list_tmp$data$rho_gp_priors = c(prior_combo[prior_i,5],prior_combo[prior_i,6])
          data_list_tmp$data$alpha_gp_priors = c(prior_combo[prior_i,7],prior_combo[prior_i,8])
          
          
          
          library(cmdstanr)
          file <- file.path(cmdstan_path(), "pinnipeds", cmd_file_name)
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

          fit2_sum = fit2$summary()
          # save model diags
          diag_summ = fit2$diagnostic_summary()
          divergent_perc = sum(diag_summ$num_divergent)/(mcmc_list$n_chain*mcmc_list$n_mcmc)
          
          eff_sample_size = fit2_sum$ess_tail/(mcmc_list$n_chain*mcmc_list$n_mcmc)
          n_eff0.01_tmp = sum(eff_sample_size<0.01,na.rm=T)/length(eff_sample_size)
          
          rhat1.1_tmp = sum(fit2_sum$rhat>1.01,na.rm=T)/length(fit2_sum$rhat)
          
          loo_tmp = loo(fit2$draws("log_lik"))
          
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
          
          
          if (model_type == "GP") {
            gaussian_rel_out = fit2$draws("gaussian_rel",format = "matrix")
            gaussian_rel_tmp = data.frame(tot.q10 = apply(gaussian_rel_out,2,quantile, probs = c(0.10)),
                                          tot.q50 = apply(gaussian_rel_out,2,quantile, probs = c(0.50)),
                                          tot.q90 = apply(gaussian_rel_out,2,quantile, probs = c(0.90)),
                                          var="gaussian_rel",species = species_i) %>%
              mutate(x = 1:1000)
            gaussian_rel = rbind(gaussian_rel,gaussian_rel_tmp) 
            saveRDS(gaussian_rel,gaussian_file_name)
          }
          # ggplot(data=gaussian_rel_tmp,aes(x=x,y=tot.q50,ymin=tot.q10,ymax=tot.q90)) + geom_ribbon(alpha=0.2) + geom_line() +
          #      xlim(c(0,1500))
          
          
          
          loo_list_name = paste(species_i,model_type,sep="-")
          loo_all[[loo_list_name]] = loo_tmp
          
          saveRDS(loo_all,loo_file_name)
          
          model_diags = data.frame(species = species_i, #rmu = rmu_i,proc_state = Q, 
                                   rhat1.1=rhat1.1_tmp,n_eff0.01=n_eff0.01_tmp,div_p=divergent_perc,
                                   loo = loo_tmp$estimates["looic","Estimate"],
                                   loo_se = loo_tmp$estimates["looic","SE"]) %>% 
            mutate(model_type = model_type,prior = prior_i)
          mod_sel_grid = rbind(mod_sel_grid,model_diags) 
          
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
          
          # join raw data to plot 
          raw_data = species_i_df %>% data.frame() %>% 
            set_colnames(colnames(species_i_df)) %>%
            rownames_to_column("ts") %>% 
            pivot_longer(!ts,
                         names_to="year", 
                         values_to = "data") %>% 
            mutate(year=as.numeric(year)) 
          
          
          # extract PRED  from commandstanr
          preds = fit2$summary(variables = "pred")
          preds_matrix = cmdstanr_preds_sd(preds,species_i_df)
          pred_data_tmp = preds_matrix %>% left_join(raw_data) %>% mutate(species = species_i, model_type = model_type,prior = prior_i) %>%
            left_join(tot_sum)
          
          p = pred_data_tmp %>% ggplot() + 
            geom_ribbon(aes(x=year,ymin = mean-1.96*sd, ymax = mean+1.96*sd), alpha=0.2) +
            # geom_line(aes(x=year,y=state_fit)) +
            geom_line(aes(x=year,y=(mean))) +
            geom_point(aes(x=year,y=(data)),color="red")+
            facet_wrap(~ts,scales="free") + theme_bw()
          print(p)
          
          pred_data=rbind(pred_data,pred_data_tmp)
          
          saveRDS(mod_sel_grid,mod_sel_file_name)
          saveRDS(error_params,error_params_file_name)
          saveRDS(pred_data,pred_data_file_name)
          
          
          
        }
      }
    }
  }}

