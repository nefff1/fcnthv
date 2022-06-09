# Estimation of mass loss to herbivory on main tree species ####################.

# restricted to (close-to) monocultures
# based on LAI, SLA and herbivory rates

# setup ########################################################################
################################################################################.

library(tidyverse); options(dplyr.summarise.inform = FALSE)
library(car)

# read data ####################################################################
################################################################################.

# tree layer cover
# from tree and regeneration inventory
d.tree <- read.table("d.tree.txt", header = T)

# plant traits 
# condensed from BExIS ID 24807
d_p_traits_f <- read.table("d_p_traits_f.txt", header = T)

# leaf area index (LAI)
# BExIS ID 24886
d_lai <- read.table("d_lai.txt", header = T)

# plant-species-level data (herbivory rates)
d_spec <- read.table("Data/Dat_specieslevel.txt", header = T)

# massloss estimation ##########################################################
################################################################################.

# determine (close-to, i.e. > 90% one tree species) stands
monoculture_ids <- d.tree %>% 
  group_by(PlotID) %>% 
  mutate(Cover_rel = Cover / sum(Cover)) %>% 
  filter(Cover_rel > .90) %>% 
  select(PlotID, Host_Name_std) %>% 
  filter(Host_Name_std %in% c("Fagus sylvatica", "Picea abies", "Pinus sylvestris")) %>% 
  mutate(monoculture = T)

# aggregate SLA
d_p_traits_f_m <- d_p_traits_f %>% 
  group_by(plotid_withzero, Host_Name_std) %>% 
  summarise(sla = mean(sla))

# caluclate massloss
d_massloss <- d_spec %>% 
  filter(system == "forest") %>% 
  left_join(monoculture_ids, by = c(plot = "PlotID", "Host_Name_std")) %>% 
  filter(monoculture) %>% 
  left_join(d_p_traits_f_m, by = c(plot = "plotid_withzero", "Host_Name_std")) %>% 
  left_join(d_lai, by = c(plot = "PlotID")) %>% 
  mutate(massloss = LAI * 10^6 / (sla * 100) * a.herb2.prop) # unit masslos: g/m2

# plotting #####################################################################
################################################################################.

p1 <- d_massloss %>% 
  ggplot(aes(x = Host_Name_std, y = massloss)) +
  geom_boxplot() +
  ylab(expression(Estimated~mass~loss~(g/m^2))) +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(NA, 50))

p2 <- d_massloss %>% 
  ggplot(aes(x = Host_Name_std, y = massloss)) +
  geom_point() +
  theme(axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(119.5, 124), breaks = 120)

plot_grid(p2, p1, ncol = 1, rel_heights = c(1, 7), align = "v")

# stats ########################################################################
################################################################################.

mod_massloss <- lm(massloss ~ Host_Name_std, d_massloss)
Anova(mod_massloss, type = "III")