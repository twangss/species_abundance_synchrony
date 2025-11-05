library(tidyverse)
library(ggplot2)
library(sf)

# Spatiotemporal representation plots
LPD2022_public_long_turtle = read.csv("/Users/twangs/Documents/Github/LPI/LPD2022_public/LPD2022_public.csv") %>% 
  filter(Common_name == "Green turtle" | Common_name == "Hawksbill turtle" | Common_name == "Leatherback turtle" | Common_name == "Loggerhead turtle" | Common_name == "Olive ridley turtle") %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "year",
    names_prefix = "X",
    values_to = "popvalue",
    values_drop_na = TRUE) %>%
  mutate(popvalue=as.numeric(popvalue),Year=as.numeric(year)) %>%
  drop_na(popvalue) %>%
  group_by(Year) %>%
  summarise(n=n()) %>%
  mutate(Source = "LPI")

data_coverage_plot = turt_df_long %>% # same turt_df_long from "turtle model fitting and selection.R"
  group_by(Year) %>%
  summarise(n=n()) %>%
  mutate(Source = "Our Data") %>%
  bind_rows(LPD2022_public_long_turtle) %>%
  ggplot(aes(x=Year,y=n,color=Source)) +
  geom_line(size=1.2) +
  theme_classic() +
  labs(y="Data Points",x="Year") +
  scale_color_manual(values=c("Our Data"="steelblue","LPI"="firebrick")) +
  theme(legend.position = c(0.1,0.8), legend.background = element_rect(fill = "white", color = "grey80"))

turt_df_sf = turt_df_long %>%
  ungroup() %>%
  dplyr::select(CommonName, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude","Latitude"))%>% 
  st_set_crs(4326) #sets the coordinates to WGS 84
  
library(ochRe)
# get world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent != "Antarctica")

world_map = ggplot() +
  geom_sf(data = world, fill = "grey70", color = NA) + # World map
  geom_sf(data = turt_df_sf, aes(fill = CommonName,shape = CommonName),color="gray20", alpha = 0.2,size=4) +   
  theme_classic() +
  theme(
    # plot.margin=unit(c(0,0,0,0), "pt"),
    legend.position = c(0.15,0.2),
    legend.text = element_text(size = 8),
    # legend.box = "vertical", # Set legend box arrangement
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    axis.ticks = element_blank(),axis.text = element_blank(),axis.line = element_blank()
  )+
  scale_fill_ochre(palette = "winmar",name="Species")+ 
  scale_shape_manual(values=c(21, 21, 21,21,23,23,23))+
  labs(fill = "Species", shape = "Species")+
  guides(
    fill  = guide_legend(override.aes = list(alpha = 1),ncol = 2),
    shape = guide_legend(override.aes = list(alpha = 1),ncol = 2)
  ) +
  coord_sf(xlim = c(-165, 170),ylim = c(-55, 80))



spatiotemporal_plot <- cowplot::plot_grid(data_coverage_plot, world_map, 
                           labels = c("a", "b"),
                           ncol = 1)

ggsave("plots/turtle_spatiotemporal_coverage.pdf", spatiotemporal_plot, width=8, height=7)

