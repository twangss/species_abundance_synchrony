# # set up state hypotheses
# stock_vec <- species_i_pre_df %>% subset(Common.Name == species_i) %>% distinct(ID,Stock)%>% pull(Stock) 
# custom_vec <- species_i_pre_df %>% subset(Common.Name == species_i) %>% distinct(ID,Custom) %>% pull(Custom) 
# id_vec <- species_i_pre_df %>% subset(Common.Name == species_i) %>% distinct(ID)%>% pull(ID) 
# onestate_vec <- rep("one",times=length(id_vec))
library(tidyverse)
library(atsar)

marss_stan = function(species_i_df_in, state_hyp_in, Q_in){
  
  state_vec_tmp = eval(parse(text=state_hyp_in))
  n_states <- max(as.numeric(factor(state_vec_tmp)))

  r_unequal <- seq(1, nrow(species_i_df))
  r_equal <- rep(1, nrow(species_i_df))
  p_uneq <- seq(1, n_states)
  p_eq <- rep(1, n_states)
  
  pro_tmp = eval(parse(text=Q_in))
  
f_tmp <- fit_stan(y = species_i_df,
                  est_nu = F,
                  est_trend = F,
                  mcmc_list = list(n_mcmc = 3000, n_burn = 500, n_chain = 1, n_thin = 1),
                  marss = list(states = as.numeric(factor(state_vec_tmp)), 
                               obsVariances = r_equal, 
                               proVariances = pro_tmp,
                               trends = p_eq))
return(f_tmp)
}

get_preds = function(species_i_df_in, pars_in,species_name){
  species_i_df =  species_i_df_in
  pars = pars_in
  # extract outputs into long format
  state_fits = apply(pars$pred, c(3,2), mean) %>% data.frame() %>% 
    set_colnames(colnames(species_i_df)) %>%
    magrittr::set_rownames(rownames(species_i_df)) %>%  
    rownames_to_column("ts") %>%
    pivot_longer(!ts,
                 names_to="year", 
                 values_to = "state_fit")
  
  # clunky code to get out means and CI intervals
  state_lo = apply(pars$pred, c(3,2), quantile, 0.025)%>% data.frame() %>% 
    set_colnames(colnames(species_i_df)) %>%
    set_rownames(rownames(species_i_df)) %>%  
    rownames_to_column("ts") %>%
    pivot_longer(!ts,
                 names_to="year", 
                 values_to = "state_lo")
  
  state_hi = apply(pars$pred, c(3,2), quantile, 0.975)%>% data.frame() %>% 
    set_colnames(colnames(species_i_df)) %>%
    set_rownames(rownames(species_i_df)) %>%  
    rownames_to_column("ts") %>%
    pivot_longer(!ts,
                 names_to="year", 
                 values_to = "state_hi")
  
  # join raw data to plot later
  raw_data = species_i_df %>% data.frame() %>% 
    set_colnames(colnames(species_i_df)) %>%
    rownames_to_column("ts") %>% pivot_longer(!ts,
                                              names_to="year", 
                                              values_to = "data")
  
  state_df = state_fits %>% left_join(state_lo) %>% left_join(state_hi)  %>% 
    left_join(raw_data) %>% mutate(year=as.numeric(year),
                                   species = species_name)  
  
  return(state_df)
}


setup_data <- function(y, x = NA,
                       mcmc_list = list(n_mcmc = 1000, n_burn = 500, n_chain = 3, n_thin = 1),
                       family = "gaussian",
                       est_nu = FALSE,
                       est_trend = FALSE,
                       # est_sigma_process_prior = FALSE,
                       marss = list(states = NULL, obsVariances = NULL, proVariances = NULL, trends = NULL,
                                    stocks = NULL,
                                    # sigma_process_prior=NULL
                       ),
                       map_estimation = FALSE,
                       hessian = FALSE, ...) {
  dist <- c("gaussian", "binomial", "poisson", "gamma", "lognormal")
  family <- which(dist == family)
  
  # process potential NAs in data
  if (!is.matrix(y)) {
    N <- length(y)
    pos_indx <- which(!is.na(y))
    y <- y[pos_indx]
    n_pos <- length(pos_indx)
    # catch case where -- needs to be 2 elements for stan vec
    if (length(pos_indx) == 0) {
      pos_indx <- rep(0, 2)
    } else {
      pos_indx <- c(pos_indx, 0, 0)
    }
    
  }
  
  data <- NA
  
  # object <- stanmodels$marss
  if (is.null(marss$states)) marss$states <- rep(1, nrow(y))
  if(length(marss$states) != nrow(y)) stop("Error: state vector must be same length as number of time series in y")
  if (is.null(marss$obsVariances)) marss$obsVariances <- rep(1, nrow(y))
  if(length(marss$obsVariances) != nrow(y)) stop("Error: vector of observation error variances must be same length as number of time series in y")
  if (is.null(marss$proVariances)) marss$proVariances <- rep(1, max(marss$states))
  if(length(marss$proVariances) < max(marss$states)) stop("Error: vector of process error variances is fewer than the number of states")
  if(length(marss$proVariances) > max(marss$states)) stop("Error: vector of process error variances is larger than the number of states")
  if (is.null(marss$trends)) marss$trends <- rep(1, max(marss$states))
  if(length(marss$trends) < max(marss$states)) stop("Error: vector of trends is fewer than the number of states")
  if(length(marss$trends) > max(marss$states)) stop("Error: vector of trends is larger than the number of states")
  
  proVariances <- c(marss$proVariances, 0) # to keep types in stan constant
  trends <- c(marss$trends, 0) # to keep types in stan constant
  N <- ncol(y)
  M <- nrow(y)
  row_indx_pos <- matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
  col_indx_pos <- matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
  n_pos <- length(row_indx_pos)
  y <- y[which(!is.na(y))]
  
  est_A <- rep(1, M)
  for(i in 1:max(marss$states)) {
    indx <- which(marss$states==i)
    est_A[indx[1]] <- 0
  }
  est_A <- which(est_A > 0)
  est_A <- c(est_A, 0, 0)
  n_A <- length(est_A) - 2
  
  data = list("N"=N,"M"=M, "y"=y,
              "states"=marss$states, "S" = max(marss$states), "obsVariances"=marss$obsVariances,
              "n_stocks" = max(marss$stocks),
              "stocks" = marss$stocks,
              "n_obsvar" = max(marss$obsVariances), "proVariances" = proVariances,
              "n_provar" = max(proVariances),
              "trends"=trends, "n_trends" = max(trends),
              "n_pos" = n_pos,
              "col_indx_pos" = col_indx_pos,
              "row_indx_pos" = row_indx_pos,
              "est_A" = est_A,
              "n_A" = n_A,
              "est_nu" = est_nu,
              "est_trend" = est_trend,
              # "est_sigma_process_prior" = est_sigma_process_prior,
              # "sigma_process_prior" = marss$sigma_process_prior,
              "family"=1)
  
  #pars = c("pred", "log_lik","sigma_process","sigma_obs","x0")
  pars = c("pred", "sigma_process","sigma_obs","x0", "log_lik")
  #if(marss$est_B) pars = c(pars, "B")
  if(est_trend) pars = c(pars, "U")
  if(n_A > 0) pars = c(pars,"A")
  if(est_nu) pars = c(pars,"nu")
  
  return(list(data = data, pars = pars))
  #return(list(model = out, data = data, pars = pars))
}

cmdstanr_preds_sd = function(preds_in,species_i_df_in){
  preds_matrix = matrix(preds_in$mean,ncol=nrow(species_i_df_in),nrow=ncol(species_i_df_in)) %>% 
    t() %>%
    set_colnames(colnames(species_i_df_in)) %>%
    set_rownames(rownames(species_i_df_in)) %>%
    data.frame() %>%
    rownames_to_column("ts") %>%
    pivot_longer(!ts,
                 names_to="year", 
                 names_prefix="X",
                 values_to = "mean") %>%
    mutate(year=as.numeric(year))
  
  sd_matrix = matrix(preds_in$sd,ncol=nrow(species_i_df_in),nrow=ncol(species_i_df_in)) %>% 
    t() %>%
    set_colnames(colnames(species_i_df_in)) %>%
    set_rownames(rownames(species_i_df_in)) %>%
    data.frame() %>%
    rownames_to_column("ts") %>%
    pivot_longer(!ts,
                 names_to="year", 
                 names_prefix="X",
                 values_to = "sd") %>%
    mutate(year=as.numeric(year))
  
  preds_sd_long = preds_matrix %>% left_join(sd_matrix, by = c("year","ts"))
  return(preds_sd_long)
}