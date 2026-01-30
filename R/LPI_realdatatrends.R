library(devtools)
library(rlpi)
library(tidyverse)

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

lpi_tbl_all=NULL
    
species_loop = unique(turt_df_long$Subspecies)

  for(species_i in species_loop){
      LPI_D = turt_df_long %>% filter(Subspecies == species_i) %>% 
        # filter(Obs_year ==T) %>%
        select(Site, ID, Year, Count) %>% 
        mutate(Site = str_replace_all(Site, "[^[:alnum:]]", "")) %>%
        arrange(ID,Year)
      
      colnames(LPI_D) = c("Bionomial","ID","year","popvalue")
      
      ID_loop = unique(LPI_D$ID)
      data_file_names = paste("data/LPI_files/","lpi_test_",ID_loop,".txt",sep="")
      
      for (j in 1:length(ID_loop)){
        ID_i = ID_loop[j]
        LPI_D_tmp = LPI_D %>% filter(ID == ID_i)
        
        data_file_save_i = data_file_names[j]
        write.table(LPI_D_tmp,data_file_save_i,sep="\t",row.names=FALSE)
      }
      
      LPI_D_wts = LPI_D %>% group_by(ID) %>% 
        summarise(mean_ind = mean(popvalue,na.rm=T)) %>%
        ungroup() %>%
        mutate(mean_ind_wt = mean_ind/sum(mean_ind,na.rm=T)) %>%
        arrange(ID) %>% pull(mean_ind_wt)
      
      for(weights_opt in c(0,1)){
        if (weights_opt == 0){weights_vec = 1} else
        {weights_vec = LPI_D_wts} 
        
        lpi_infile = data.frame(FileName = data_file_names, Group = 1, Weighting = weights_vec)
        
        lpi_infile_path = paste("data/LPI_files/","lpi_infile_test",".txt",sep="")
        write.table(lpi_infile,lpi_infile_path,sep="\t",row.names=FALSE)
        
        seal_lpi_tmp <- LPIMain("data/LPI_files/lpi_infile_test.txt", use_weightings = weights_opt, VERBOSE=F,plot_lpi=1,REF_YEAR=(min(LPI_D$year+0)))
        
        seal_lpi_tmp_attr = seal_lpi_tmp %>% 
          mutate(species = species_i,year=rownames(seal_lpi_tmp),weights = weights_opt,
                 sim = sim_i , true_state = model_i)
        lpi_tbl_all = rbind(lpi_tbl_all,seal_lpi_tmp_attr)
        
      }
    }

lpi_tbl_all %>% 
  filter(weights == 1) %>%
  mutate(year=as.numeric(year),weights = as.factor(weights)) %>%
  filter(LPI_final > -99) %>%
  ggplot() + geom_line(aes(x=year,y=LPI_final, color = as.factor(weights))) + facet_wrap(~species,scales="free")

saveRDS(lpi_tbl_all,"model/outputs/lpi_tbl_realdata_2026_0115.RDS")


# compare with MARSS trends
subspp = readRDS("model/outputs/seal_realreal0_supspp_GPcorr_prior_2026-01-06.RDS")
IUCN_cutoff = read.csv("data/IUCN_seal_names_popsize.csv") %>%
  select(CommonName, RMU, CommonNamespp,plot_main)

subspp_iucn = subspp %>% 
  inner_join(IUCN_cutoff, by=c("species" = "CommonName", "RMU" = "RMU")) %>%
  filter(plot_main == 1) %>%
  group_by(species,CommonNamespp) %>%
  mutate(meanscale = `50%`/mean(`50%`)) %>%
  ungroup() %>%
  select(CommonNamespp,year,meanscale,model_type) 

lpi_tbl_all = readRDS("model/outputs/lpi_tbl_realdata_2026_0115.RDS")
lpi_to_join_with_MARSS = lpi_tbl_all %>%
  filter(LPI_final > 0) %>%
  filter(weights == 1) %>%
  group_by(species) %>%
  mutate(meanscale = LPI_final/mean(LPI_final)) %>%
  ungroup() %>%
  mutate(model_type = "LPI",year=as.numeric(year)) %>%
  inner_join(IUCN_cutoff, by=c("species" = "RMU")) %>%
  filter(plot_main == 1) %>%
  select(CommonNamespp,year,meanscale,model_type)

MARSS_LPI_species_trend = rbind(subspp_iucn,lpi_to_join_with_MARSS) %>%
  mutate(model_type = case_when(model_type == "GP" ~ "Synchrony Model",
                                model_type == "Ueq" ~ "Independent Model",
                                model_type == "LPI" ~ "LPI"))


cols <- c("Synchrony Model" = "#186689", "LPI" = "forestgreen", 
          "Independent Model" = "#e23b0d")


p = MARSS_LPI_species_trend %>%
  ggplot(aes(color = model_type, linetype = model_type)) +
  # scaled units
  # geom_point(aes(x=year,y=data_norm_trans),alpha=0.2)+
  # geom_line(aes(x=year,y=data_norm_trans,group=ID),alpha=0.1)+
  # geom_ribbon(aes(x=year,ymin = pmax((x.tot.norm-year_species_sd)/y_axis_divide,0), ymax = (x.tot.norm+year_species_sd)/y_axis_divide), alpha=0.4) +
  geom_line(aes(x=year,y=(meanscale)),lwd=0.9) +
  scale_x_continuous(breaks=scales::pretty_breaks())+
  scale_y_continuous(breaks=scales::pretty_breaks())+
  labs(x="Year",y=paste("Mean scaled abundance",sep="")) +
  scale_color_manual(values = cols) + 
  # geom_hline(yintercept=0, linetype = "dashed", alpha=0.3) +
  facet_wrap(~CommonNamespp, scales = "free_y", ncol = 5) +
  labs(color = "Model Type", linetype = "Model Type") +
  expand_limits(y = c(0)) +
  xlim(c(1965,2022)) +
  theme_classic() +
  theme(legend.position = "bottom")
p

pdf("model/plots/seal_compare_model_species_trends.pdf", width = 8.5, height = 11)
print(p)
dev.off()
