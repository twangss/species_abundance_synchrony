library(tidyverse)
library(ggplot2)

subspp = readRDS("model/outputs/seal_realreal0_supspp_GPcorr_prior_2026-01-06.RDS")
IUCN_cutoff = read.csv("data/IUCN_seal_names_popsize.csv") %>%
  select(CommonName, RMU, CommonNamespp, GL_days,IUCN.cat, IUCN.pop,IUCN.pop.year,plot_main)
model_grid_best = readRDS("model/outputs/model_grid_best_20260106.RDS")

# plot species trends
RMU_species_number_ts = turt_df_long %>% ungroup()%>%
  # dplyr::select(CommonName,ID,Units_clean) %>% 
  # distinct(CommonName,ID,Units_clean) %>%
  group_by(CommonName,RMU) %>% 
  mutate(species_n = n_distinct(ID)) %>% 
  # filter(Indiv.sum == 1) %>%
  filter(Scaling==1) %>% # to get the common denominator untis
  dplyr::select(CommonName,RMU,species_n,Units_clean) %>% distinct()

subspp_best = subspp %>% 
  inner_join(model_grid_best, by = c("species","model_type")) %>%
  inner_join(IUCN_cutoff, by=c("species" = "CommonName", "RMU" = "RMU"))

coverage = subspp_best %>% 
  group_by(species,RMU) %>%
  filter(year == IUCN.pop.year) %>%
  mutate(coverage = min(round(`50%`/IUCN.pop*100),100)) %>%
  select(species,RMU,coverage)


subspp_label = subspp_best %>% 
  left_join(coverage) %>%
  left_join(RMU_species_number_ts, by = c("species" = "CommonName","RMU"="RMU")) %>%
  mutate(Units_clean = case_when(Units_clean == "All Individuals" ~ "I",
                                 Units_clean == "Pups" ~ "P",
                                 Units_clean == "Adults" ~ "A",
                                 Units_clean == "Non-pups" ~ "NP"))

subspp_label$species_label <- paste(subspp_label$CommonNamespp,"<br>",
                                    "n<sub>p</sub> = ",subspp_label$species_n,", ",
                                    "cov = ",subspp_label$coverage,"%, ",
                                    subspp_label$Units_clean,
                                    sep="")  

subspp_label$species_label <- iconv(subspp_label$species_label, to = "UTF-8", sub = "")


# pdf("model/plots/seals_species_trends.pdf", width = 10, height = 8)
iucn_order <- c("LC", "NT", "VU", "EN", "CR", "DD")
subspp_label$IUCN.cat <- factor(subspp_label$IUCN.cat, levels = iucn_order)
y_axis_divide= 1000
p=subspp_label %>% filter(plot_main ==1) %>%
  ggplot(aes(color=IUCN.cat,fill=IUCN.cat)) + 
  geom_ribbon(color=NA,aes(x=year,ymin = `25%`/y_axis_divide, ymax = `75%`/y_axis_divide), alpha=0.4) +
  geom_line(aes(x=year,y=(`50%`/y_axis_divide)),lwd=0.9) +
  # scale_color_manual(values = c( "#610B0B","#0E5917","#4D4F70"))+
  # scale_fill_manual(values = c("#610B0B","#0E5917","#4D4F70"))+
  scale_x_continuous(breaks=scales::pretty_breaks(n=4))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  labs(x="Year",y="Monitored Abundance (1000s units)",
       color = "IUCN category", 
       fill = "IUCN category") +
  # geom_hline(yintercept=0, linetype = "dashed", alpha=0.3) +
  facet_wrap(~species_label, scales = "free_y", ncol = 5) +
  # expand_limits(y = c(0,data_in$ymax_val)) + 
  coord_cartesian(y=c(0,NA),x=c(1970,2024)) + 
  # xlim(c(1965,2022)) +
  # theme_bw() +
  scale_fill_manual(values=c("VU"="#CD9A02",
                             "EN"="#cd6630",
                             "CR"="#D62828",
                             "NT"="#016766",
                             "LC"="#558B2F",
                             "DD"= "grey30")) +
  scale_color_manual(values=c("VU"="#CD9A02",
                              "EN"="#cd6630",
                              "CR"="#D62828",
                              "NT"="#016766",
                              "LC"="#558B2F",
                              "DD"= "grey30")) +
  theme(strip.background = element_rect(colour=NA, fill=NA),
        strip.text = ggfacet::element_textbox_highlight(
          size = 9, face = "plain",
          fill = NA, box.color = NA, color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)
        ),
        panel.border = element_rect(colour = "gray10", fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "gray10"),
        axis.text.y = element_text(color = "gray10"),
        axis.title.x = element_text(color = "gray10"),
        axis.title.y = element_text(color = "gray10"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "bottom",legend.direction="horizontal",
        # plot.title = element_text(hjust = .01, vjust=-7),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))
print(p)
# dev.off()

file_list <- list.files(path = "model/outputs", pattern = "^posterior_r_.*\\.RDS$",full.names = TRUE)

posterior_data <- file_list %>%
  setNames(., .) %>%  # FIX: Sets the names of the list to the filenames themselves
  map_dfr(readRDS) %>% 
  as_tibble()

cutoff = posterior_data %>%
  inner_join(model_grid_best, by = c("species","model_type")) %>%
  inner_join(IUCN_cutoff, by=c("species" = "CommonName", "RMU" = "RMU")) %>%
  filter(plot_main == 1) %>%
  group_by(CommonNamespp) %>%
  mutate(
    year_range = max(year) - min(year),
    max_year_ref = max(year),
    IUCN_cutoff_yr = floor(max(year) - (3*GL_days)/365)) %>%
  filter(year_range > 10) %>%
  filter(year == IUCN_cutoff_yr | (!any(year == IUCN_cutoff_yr) & year == min(year))) %>%
  mutate(duration_years = max_year_ref - year,
         lambda = (exp(r / duration_years) - 1) * 100)

lambda_median = cutoff %>%
  group_by(CommonNamespp,IUCN.cat) %>%
  summarise(
    min_year = min(year),
    median_lambda = median(lambda, na.rm = TRUE))


pdf("model/plots/lambda_IUCN_status.pdf", width = 3.5, height = 3)
iucn_order <- c("LC", "NT", "VU", "EN", "CR", "DD")
lambda_median$IUCN.cat <- factor(lambda_median$IUCN.cat, levels = iucn_order)
p=ggplot(data=lambda_median,aes(x=median_lambda,y=IUCN.cat, color=IUCN.cat)) + 
  geom_vline(xintercept = 0, color = "black",alpha=0.8, linetype = "dashed") + 
  geom_boxplot(alpha=1,outlier.shape = NA) + 
  # geom_jitter(aes(size=year_tot),alpha=0.8) + 
  geom_jitter(size=3,alpha=0.8) +
  theme_classic(base_size = 8) +  
  scale_color_manual(values=c("VU"="#CD9A02",
                              "EN"="#cd6630",
                              "CR"="#D62828",
                              "NT"="#016766",
                              "LC"="#558B2F",
                              "DD"= "grey30")) +
  theme(
    legend.position = c(1,1), # top left position
    legend.justification = c(1.2, 1.1)) + 
  labs(x = "Total species abundance growth rate (%)",y = "IUCN Category",
       color = "IUCN Category", 
       fill = "IUCN Category")
print(p)
dev.off()


pred_data_file_name = "model/outputs/seal_realreal0_pred_data_GPcorr_prior_2026-01-06.RDS"
pred_data = readRDS(pred_data_file_name) %>%
  mutate(ID = as.numeric(str_extract(ts, "\\d+$"))) %>%
  left_join(turt_df_long %>% select(ID, CommonName, RMU) %>% distinct(), by="ID")


lambda_pops = pred_data %>%
  inner_join(model_grid_best, by = c("species","model_type")) %>%
  group_by(ts) %>%
  mutate(
    min_data_year = min(year[!is.na(data)]),
    max_data_year = max(year[!is.na(data)])) %>%
  filter(year >= min_data_year & year <= max_data_year) %>%
  arrange(ts,year) %>%
  mutate(growth_rate = (exp(mean)-exp(lag(mean)))/exp(lag(mean))) %>%
  group_by(species,RMU,ts) %>%
  summarise(mean_growth_rate = mean(growth_rate, na.rm=T)) %>%
  inner_join(IUCN_cutoff, by=c("species" = "CommonName", "RMU" = "RMU")) %>%
  filter(plot_main == 1)

lambda_pops_plot = ggplot(data =lambda_pops,aes(x=100*mean_growth_rate,y=reorder(CommonNamespp, mean_growth_rate, median))) + 
  geom_boxplot(outliers = F) + 
  geom_jitter(size=1,alpha=0.5) +
  # geom_point() + 
  geom_vline(xintercept = 0 , alpha=0.4, linetype= "dashed") + 
  theme_bw() +
  labs(x = "Mean abundance growth rate (%) by population",y = "Species") +
  coord_cartesian(xlim = c(-5, 20))


lambda_species_plot = subspp_label %>%
  filter(plot_main == 1) %>%
  group_by(CommonNamespp) %>%
  mutate(n_year = n_distinct(year)) %>%
  filter(n_year > 9) %>%
  arrange(year) %>%
  mutate(lambda = (`50%`-lag(`50%`))/lag(`50%`)) %>%
  summarise(mean_growth_rate = mean(lambda, na.rm=T)) %>%
  ggplot(aes(x=100*mean_growth_rate,y="")) + 
  geom_boxplot(outliers = F) + 
  geom_jitter(width = 0,size=2,alpha=0.5) +
  # geom_point() + 
  geom_vline(xintercept = 0 , alpha=0.4, linetype= "dashed") + 
  theme_bw() +
  labs(x = "Mean abundance growth rate (%) by species",y = "All Species")


lambda_plots = cowplot::plot_grid(lambda_species_plot, lambda_pops_plot, ncol = 1, 
                                  rel_widths = c(1, 1),
                                  rel_heights = c(0.5, 1),
                                  labels=c("a","b"))  
pdf("model/plots/species_population_lambda.pdf", width = 8, height = 7)
print(lambda_plots)
dev.off()


#=== plot lambda trends across years
# join to generation length data, IUCN table is the middle step to join
subspp_lambda_year = subspp_best %>%
  filter(plot_main==1)%>%
  group_by(RMU) %>%
  arrange(year) %>%
  mutate(lambda = (`50%`-lag(`50%`))/lag(`50%`))


pdf("model/plots/lambda_years.pdf", width = 3.5, height = 3)
lambda_year_boxplot_n = subspp_lambda_year %>% group_by(year) %>%
  dplyr::mutate(n_species = n()) %>% filter(lambda<0.3 & lambda > -0.3)

p=ggplot(lambda_year_boxplot_n, aes(y = lambda*100, x = year)) +
  coord_cartesian(ylim=c(-5,9),xlim=c(1950,2020))+ # using coord cartesian to avoid changing boxplots
  xlim(c(1950,2024))+
  geom_hline(yintercept=0,color="grey50") + 
  geom_boxplot(aes(group=year,fill=n_species),outlier.shape=NA, coef = 0,linewidth=0.2,color="white") +
  # geom_point(aes(group=year,color=n_species),alpha=0.3) +
  viridis::scale_fill_viridis(option="mako",direction=-1) +
  stat_summary(aes(group=year),fun.y=mean, geom="point", shape=20, size=1.5, color="#E9724C") +
  geom_smooth(color="#C5283D",fill="#C5283D",se=T,level = 0.95,lwd=1.5,alpha=0.2,
              method = loess, method.args = list(family = "symmetric")) +
  theme_classic(base_size = 8) +
  labs(y="Annual growth rate (%)",x="Year",fill="Number of species") +
  # guides(guide_legend("Number of species"))+
  theme(
    legend.position = c(0,1), # top left position
    legend.justification = c(-.1, 1.1), # top left justification
    legend.key.size = unit(0.3, 'cm')
    # legend.box.margin = margin(5, l = 5, unit = "mm") # small margin
  )
print(p)
dev.off()


# Observation error plots
error_params_file_name = "model/outputs/seal_realreal0_error_params_GPcorr_prior_2026-01-06.RDS"
error_params = readRDS(error_params_file_name)

pdf("model/plots/observation_error_survey_method.pdf", width = 3.5, height = 3)
survey_plot = error_params %>% 
  inner_join(model_grid_best, by = c("species","model_type")) %>%
  filter(var == "sigma_obs") %>% filter(units != "Unknown") %>%
  ggplot(aes(x=q50,y=reorder(units, q50, median))) +
  geom_boxplot(alpha=1,outlier.shape = NA) + 
  geom_jitter(alpha=0.5,size=2) + theme_classic(base_size = 8) +
  labs(x="Log-space observation error variance", y= "Survey method")
print(survey_plot)
dev.off()

