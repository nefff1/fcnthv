################################################################################.
################################################################################.
# SETUP ########################################################################
################################################################################.
################################################################################.

library(tidyverse); options(dplyr.summarise.inform = FALSE)
library(FD)
library(furrr)
library(hillR)
library(cluster)
library(TPD)
library(ape)
library(parallel)
library(vegan)
library(ecodist)

# Source code for functional diversity from Chao et al. (2018)
# function FunD
code_lines <- scan("https://raw.githubusercontent.com/AnneChao/FunD/master/FunD_Rcode.txt", 
                   what = "character", 
                   skip = 38, nlines = 164, sep = "\n")
code_lines <- paste(code_lines, collapse = "\n")
source(textConnection(code_lines))

# function FD_MLE
code_lines <- scan("https://raw.githubusercontent.com/AnneChao/FunD/master/FunD_Rcode.txt", 
                   what = "character", 
                   skip = 13, nlines = 23, sep = "\n")
code_lines <- paste(code_lines, collapse = "\n")
source(textConnection(code_lines))
rm(code_lines)
closeAllConnections()


################################################################################.
# define functions #############################################################
################################################################################.

# extract part of int db for a plant species -----------------------------------.

f_int_extract_plant <- function(plant_i, poly_all = T, poly_sub = T){
  
  int_target <- int.db.leaves %>% 
    filter(Host_Name_std == plant_i)
  
  # go down
  
  children <- host.std %>% 
    filter(Host_IsChildTaxonOf == plant_i)
  
  if (nrow(children) > 0){
    goon <- T
    while (goon){
      int_target <- int.db.leaves %>% 
        filter(Host_Name_std %in% children$Host_Name_std) %>% 
        bind_rows(int_target, .)
      
      children <- host.std %>% 
        filter(Host_IsChildTaxonOf %in% children$Host_Name_std)
      
      if (nrow(children) == 0){
        goon <- F
      }
    }
    
  }
  
  # go up
  
  parent <- host.std %>% 
    filter(Host_Name_std == plant_i) %>% 
    select(Host_IsChildTaxonOf) %>% 
    distinct() %>% 
    deframe()
  
  goon <- T
  while (goon){
    int_target <- int.db.leaves %>% 
      filter(Host_Name_std == parent) %>% 
      bind_rows(int_target, .)
    
    parent <- host.std %>% 
      filter(Host_Name_std == parent) %>% 
      select(Host_IsChildTaxonOf) %>% 
      distinct() %>% 
      deframe()
    
    if (parent == "ROOT"){
      goon <- F
    }
  }
  
  # polyphagous
  
  if (poly_sub){
    gr_target <- pl.groups[pl.groups$Host_Name_std == plant_i, 
                           c(grepl("Host_Group", names(pl.groups)))] %>% 
      unlist()
    
    gr_target <- gr_target[!is.na(gr_target)]
    
    gr_target <- paste("polyphagous -", gr_target)
    
    int_target <- int.db.leaves %>% 
      filter(Host_Name_std %in% gr_target) %>% 
      bind_rows(int_target, .)
    
  }
  
  if (poly_all){
    int_target <- int.db.leaves %>% 
      filter(Host_Name_std == "polyphagous") %>% 
      bind_rows(int_target, .)
    
  }
  
  int_target
}

################################################################################.
# read data ####################################################################
################################################################################.

options(stringsAsFactors = FALSE)

# herbivory rates aggregated at community level + land use ---------------------.
# condensed combination of BExIS ID 20826 (RW, HW), 24646, 24806, LUI

herb.dat.comm <- read.table("Data/herb.dat.comm.txt", header = T)

# herbivory aggregated at species level + land use -----------------------------.
# condensed combination of BExIS ID 20826 (RW, HW), 24646, 24806, LUI, 21426, 23686, 24247

herb.dat.agg <- read.table("Data/herb.dat.agg.txt", header = T)

# plant cover ------------------------------------------------------------------.
# condensed combination of BExIS ID 21426, 23686, 24247

plant.cover <- read.table("Data/plant.cover.txt", header = T)

# plant biomass grasslands -----------------------------------------------------.
# BExIS ID 24166

biomass_grasslands <- read.table("Data/biomass_grasslands.txt", header = T)

# plant groups for interactions data base polyphags ----------------------------.

pl.groups <- read.table("Data/pl.groups.txt", header = T)

# standardized host plant taxonomy ---------------------------------------------.
# based on GermanSL v. 1.5

host.std <- read.table("Data/host.std.txt", header = T)

# insect abundance grasslands (species level) ----------------------------------.
# condensed version BExIS ID 21969

swnet_agg <- read.table("Data/swnet_agg.txt", header = T)

# insect abundance forests (species level) -------------------------------------.
# condensed version BExIS ID 22008

wintr_agg <- read.table("Data/wintr_agg.txt", header = T)

# insect abundance  (higher taxon level) ---------------------------------------.

herbivore.abundance <- read.table("Data/herbivore.abundance.txt", header = T)

# interaction data base --------------------------------------------------------.
# restrict to interactions on (potentially) leaves, condensed & filtered version of BExIS ID 26926
int.db.leaves <- read.table("Data/int.db.leaves.txt", header = T)

# insect traits ----------------------------------------------------------------.
# from Gossner et al. (2015), ScientificData, DOI 10.1038/sdata.2015.13

traits_lh <- read.table("Data/traits_lh.txt", header = T)

# Insect length / weight allometric relation -----------------------------------.
# data from Sohlström et al. (2018), Dryad Digital Repository, DOI 10.5061/dryad.vk24fr1
# and from Sohlström et al. (2018), Ecology and Evolution, DOI 10.1002/ece3.4702

regs_mass_length <- read.table("Data/regs_mass_length.txt", header = T)

# SLA and LDMC -----------------------------------------------------------------.
# condensed from BExIS ID 24807

# grassland data set:
d_p_traits_g <- read.table("Data/d_p_traits_g.txt", header = T)
 
# forest data set:
d_p_traits_f <- read.table("Data/d_p_traits_f.txt", header = T)

# Nutrient contents ------------------------------------------------------------. 
# condensed from BExIS ID 23367 and 26526

d_nutrient_g <- read.table("Data/d_nutrient_g.txt", header = T)

# plant height -----------------------------------------------------------------.
# data from Jäger et al. (2017): Rothmaler – Exkursionsflora von Deutschland. Gefäßpflanzen: Atlasband

d_p_height <- read.table("Data/d_p_height.txt", header = T)

# leaf area index LAI ----------------------------------------------------------.
# BExIS ID 24886

d_lai <- read.table("Data/d_lai.txt", header = T)

# forest plant biomass ---------------------------------------------------------.
# estimations based on BExIS ID 21426, 23686, tree regeneration data

d_bm_f <- read.table("Data/d_bm_f.txt", header = T)

# .................................................................. #################################################
# .................................................................. #################################################
# .................................................................. #################################################

################################################################################.
################################################################################.
# ///////// ANALYSIS GRASSLAND \\\\\\\\\\ ######################################
################################################################################.
################################################################################.

################################################################################.
#  COMMUNITY LEVEL #############################################################
################################################################################.

d_g_path_plotlevel <- data.frame(PlotID = unique(herb.dat.comm$plot),
                                 stringsAsFactors = F) %>% 
  filter(grepl("G", PlotID))
     
################################################################################.
# ........ PLANTS ------------------------------- ##############################
################################################################################.

# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# dissimilarity based on plant species composition

dissim.bc.p.plot <-
  plant.cover %>% 
  filter(grepl("G", PlotID)) %>% 
  select(PlotID, Host_Name_std, cover) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  vegdist()

# perform PCoA
pcoa.bc.p.plot <- pco(dissim.bc.p.plot, negvals = "zero")

d_g_path_plotlevel <- data.frame(PlotID = labels(dissim.bc.p.plot), 
                                 pcoa.bc.p.plot$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(p.tax.comp.1 = X1,
         p.tax.comp.2 = X2,
         p.tax.comp.3 = X3) %>%
  left_join(d_g_path_plotlevel, ., by = "PlotID")


# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# functional groups ------------------------------------------------------------.
d_g_path_plotlevel <- herb.dat.agg %>% 
  filter(system == "grassland") %>% 
  group_by(plot) %>% 
  summarise(p.fgr.comp.grass  = sum(prop.cover * (pl.group == "grass")),
            p.fgr.comp.legume = sum(prop.cover * (pl.group == "legume"))) %>% 
  ungroup() %>%
  left_join(d_g_path_plotlevel, ., by = c("PlotID" = "plot"))


# SLA / LDMC -------------------------------------------------------------------.

# calculate CWMs
plant_cover_g_ct <- herb.dat.agg %>% 
  filter(system == "grassland") %>% 
  mutate(ID = paste(plot, Host_Name_std)) %>% 
  select(ID, cover, plot) %>% 
  pivot_wider(names_from = ID, 
              values_from = cover, 
              values_fill = list(cover = 0)) %>% 
  column_to_rownames("plot")

d_p_traits_g_mean <- d_p_traits_g %>% 
  mutate(ID = paste(plotid_withzero, Host_Name_std)) %>% 
  group_by(ID) %>% 
  summarise_at(vars(ldmc, sla), ~mean(., na.rm = T)) %>% 
  ungroup() %>% 
  filter(ID %in% names(plant_cover_g_ct)) %>% 
  column_to_rownames("ID")

CWM_plant_g <- functcomp(d_p_traits_g_mean, as.matrix(plant_cover_g_ct[, rownames(d_p_traits_g_mean)])) %>% 
  rownames_to_column("PlotID")

d_g_path_plotlevel <- CWM_plant_g %>%
  rename(p.fun.comp.ldmc = ldmc,
         p.fun.comp.sla = sla) %>%
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# nutrient concentrations ------------------------------------------------------.

d_g_path_plotlevel <- d_nutrient_g %>% 
  mutate(prim.fiber = NDF - ADL) %>% 
  select(PlotID, P, N, Ca, Mg, ADL, prim.fiber) %>% 
  rename(p.fun.comp.P = P, 
         p.fun.comp.N = N,
         p.fun.comp.Ca = Ca,
         p.fun.comp.Mg = Mg,
         p.fun.comp.lignin = ADL,
         p.fun.comp.prim.fiber = prim.fiber) %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_plotlevel <- plant.cover %>% 
  filter(grepl("G", PlotID)) %>% 
  group_by(PlotID) %>% 
  summarise(p.tax.abund = sum(cover)) %>% 
  ungroup() %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_plotlevel <- biomass_grasslands %>%
  rename(p.fun.abund = biomass) %>% 
  select(PlotID, p.fun.abund) %>% 
  left_join(d_g_path_plotlevel, ., by = c("PlotID"))


# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_plotlevel <- plant.cover %>% 
  filter(grepl("G", PlotID)) %>% 
  as.data.frame %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(p.tax.div.0 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

d_g_path_plotlevel <- plant.cover %>% 
  filter(grepl("G", PlotID)) %>% 
  spread(Host_Name_std, cover, fill = 0) %>%  
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(p.tax.div.1 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

d_g_path_plotlevel <- plant.cover %>% 
  filter(grepl("G", PlotID)) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(p.tax.div.2 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")


# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

m_p_traits_sample <- d_p_traits_g %>% 
  group_by(plotid_withzero, Host_Name_std) %>% 
  summarise(ldmc = mean(ldmc, na.rm = T),
            sla = mean(sla, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(ID = paste(plotid_withzero, Host_Name_std)) %>% 
  column_to_rownames("ID") %>% 
  select(-c(plotid_withzero, Host_Name_std))

dij_p <- m_p_traits_sample %>% 
  daisy(metric = "euclidean", stand = T) %>% 
  as.matrix()

g.p.fund <-
  herb.dat.agg %>%
  filter(system == "grassland") %>% 
  mutate(ID = paste(plot, Host_Name_std)) %>% 
  select(plot, ID, cover) %>% 
  filter(ID %in% labels(dij_p)[[1]]) %>% # exclude few NAs
  spread(plot, cover, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  {FunD(data = .,
        dij = dij_p[rownames(.), rownames(.)],
        tau = max(dij_p),
        q = seq(0, 2, by = 1), 
        boot = 0, 
        datatype = "abundance")$fortau}


d_g_path_plotlevel <- g.p.fund %>% 
  rename(PlotID = site) %>% 
  select(estimate, q, PlotID) %>% 
  spread(q, estimate) %>% 
  rename(p.fun.div.0 = `0`,
         p.fun.div.1 = `1`,
         p.fun.div.2 = `2`) %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

################################################################################.
# ........ INSECTS ------------------------------ ##############################
################################################################################.

# restrict to leaf-feeding herbivores ------------------------------------------.
swnet_herb_agg <- swnet_agg %>% 
  filter(Herb_Rank == "SPECIES" & Herb_Species %in% int.db.leaves$Herb_Species |
           Herb_Rank == "GENUS" & Herb_Genus %in% int.db.leaves$Herb_Genus) 

# prepare cross table ----------------------------------------------------------.

m.swnet.ct <- swnet_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID")

# exclude plots, where only one species was recorded (no effect)
m.swnet.ct <- m.swnet.ct[rowSums(m.swnet.ct > 0) > 1, ] 
m.swnet.ct <- m.swnet.ct[, colSums(m.swnet.ct) > 0]


# prepare trait data -----------------------------------------------------------.

sel.traits <- c("Body_Size", 
                "Dispersal_ability", 
                "Feeding_mode",
                "Feeding_specialization", 
                "Stratum_use_short")

m.traits <- swnet_herb_agg %>%
  select(Herb_Name_std) %>%
  distinct() %>%
  left_join(traits_lh %>%
              select(Herb_Name_std, !! sel.traits), by = "Herb_Name_std") %>%
  mutate(Dispersal_ability = as.numeric(Dispersal_ability)) %>% 
  mutate_at(vars(Feeding_mode,
                 Stratum_use_short),
            ~ as.factor(.)) %>% 
  mutate_at(vars(Feeding_specialization,
              Dispersal_ability),
         ~ as.ordered(.)) %>% 
  column_to_rownames("Herb_Name_std")


# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

dissim.bc.i.plot <- 
  swnet_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  vegdist()

# perform PCoA
pcoa.bc.plot <- pco(dissim.bc.i.plot, negvals = "zero")

d_g_path_plotlevel <- data.frame(PlotID = labels(dissim.bc.i.plot), 
                                 pcoa.bc.plot$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(i.tax.comp.1 = X1,
         i.tax.comp.2 = X2,
         i.tax.comp.3 = X3) %>%
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# perform PCA to reduce trait dimensions
ord.traits <- ade4::dudi.mix(m.traits, scannf = F, nf = 4)

# extract components for species:
ord.traits <- ord.traits$li %>%
  mutate(Herb_Name_std = rownames(m.traits))

tpd.means <- ord.traits %>%
  mutate_at(vars(-Herb_Name_std), scale) %>% # scale 
  column_to_rownames("Herb_Name_std")

tps.sds <- matrix(0.5, nrow = nrow(tpd.means), ncol = ncol(tpd.means),
                  byrow = T) # use fixed SD values of 0.5

set.seed(9)
TPDs <- TPDsMean(species = rownames(tpd.means),
                 means = tpd.means,
                 sds = tps.sds) # calculate the TPDs

TPDc <- TPDc(TPDs = TPDs , sampUnit = m.swnet.ct) # calculate the target TPDc

dissim.tpd.plot  <- dissim(TPDc)$communities$dissimilarity # calculate community dissimilarity

rm(TPDs); rm(TPDc)

# perform PCoA
pcoa.tpd <- pco(as.dist(dissim.tpd.plot), negvals = "zero")


d_g_path_plotlevel <- data.frame(PlotID = rownames(dissim.tpd.plot), 
                                 pcoa.tpd$vectors[, 1:3], 
                                 stringsAsFactors = F) %>%
  rename(i.fun.comp.1 = X1,
         i.fun.comp.2 = X2,
         i.fun.comp.3 = X3) %>%
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# what do the composition axes mean?
# cwms.i.g.plot <- functcomp(m.traits, as.matrix(m.swnet.ct)[, rownames(m.traits)])
# cwms.i.g.plot %>%
#   rownames_to_column("PlotID") %>%
#   left_join(d_g_path_plotlevel %>%
#               select(i.fun.comp.1, i.fun.comp.2, i.fun.comp.3, PlotID),
#             by = "PlotID") %>%
#   select(-PlotID) %>%
#   pairs()
# 
# cwms.i.g.plot %>%
#   rownames_to_column("PlotID") %>%
#   left_join(d_g_path_plotlevel, by = "PlotID") %>% 
#   ggplot(aes(x = Feeding_mode, y = i.fun.comp.1)) +
#   geom_boxplot()

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_plotlevel <- herbivore.abundance %>% 
  rename(i.tax.abund.sort = herbivore.abundance) %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

d_g_path_plotlevel <- swnet_herb_agg %>% 
  group_by(PlotID) %>% 
  summarise(i.tax.abund.det = sum(NumberAdults)) %>%
  ungroup() %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# ....................functional -----------------------------------------------
# ------------------------------------------------------------------------------.

# extract allometric regression coefficient for herbivores:
f_lw <- function(herb_i, data){
  genus_i <- as.character(unique(swnet_herb_agg$Herb_Genus[swnet_herb_agg$Herb_Name_std == herb_i]))
  family_i <- as.character(unique(swnet_herb_agg$Herb_Family[swnet_herb_agg$Herb_Name_std == herb_i]))
  suborder_i <- as.character(unique(swnet_herb_agg$Herb_Suborder[swnet_herb_agg$Herb_Name_std == herb_i]))
  order_i <- as.character(unique(swnet_herb_agg$Herb_Order[swnet_herb_agg$Herb_Name_std == herb_i]))
  
  if (genus_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == genus_i],
                      a = data$a[data$Herb_Name_std == genus_i],
                      stringsAsFactors = F)
    out
  } else if (family_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == family_i],
                      a = data$a[data$Herb_Name_std == family_i],
                      stringsAsFactors = F)
    out
  } else if (suborder_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == suborder_i],
                      a = data$a[data$Herb_Name_std == suborder_i],
                      stringsAsFactors = F)
    out
  } else if (order_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == order_i],
                      a = data$a[data$Herb_Name_std == order_i],
                      stringsAsFactors = F)
    out
  }
}

# apply function
regs_swnet <- unique(swnet_herb_agg$Herb_Name_std) %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

regs_swnet <- regs_swnet %>% 
  left_join(select(traits_lh, Herb_Name_std, Body_Size), by = "Herb_Name_std") 

regs_swnet <- regs_swnet %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)


d_g_path_plotlevel <- swnet_herb_agg %>% 
  left_join(regs_swnet, by = "Herb_Name_std") %>% 
  mutate(metabolic_rel_sum = metabolic_rel * NumberAdults) %>% 
  group_by(PlotID) %>% 
  summarise(i.fun.abund = sum(metabolic_rel_sum)) %>% 
  ungroup() %>%
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# q = 0
d_g_path_plotlevel <- swnet_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(i.tax.div.0 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# q = 1
d_g_path_plotlevel <- swnet_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(i.tax.div.1 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# q = 2
d_g_path_plotlevel <- swnet_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(i.tax.div.2 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

dij <- m.traits %>% 
  daisy(metric = "gower") %>% 
  as.matrix()

g.i.fund <- FunD(data = as.data.frame(t(m.swnet.ct)), 
     dij = dij[names(m.swnet.ct), names(m.swnet.ct)], 
     tau = max(dij),
     q = seq(0, 2, by = 1), 
     boot = 1, 
     datatype = "abundance")$fortau

d_g_path_plotlevel <- g.i.fund %>% 
  rename(PlotID = site) %>% 
  select(estimate, q, PlotID) %>% 
  spread(q, estimate) %>% 
  rename(i.fun.div.0 = `0`,
         i.fun.div.1 = `1`,
         i.fun.div.2 = `2`) %>% 
  left_join(d_g_path_plotlevel, ., by = "PlotID")




  
################################################################################.
# ........ FINALIZE FOR ANALYSIS ------------------------------ ################
################################################################################.

# standardize data
d_g_path_plotlevel_z <- d_g_path_plotlevel %>%
  mutate_if(is.numeric, ~ (. - mean(., na.rm = T)) / (2 * sd(.,  na.rm = T))) # scale

# .................................................................. #################################################
# .................................................................. #################################################

################################################################################.
# SPECIES LEVEL ################################################################
################################################################################.

# initiate data frames ---------------------------------------------------------.
d_g_path_samplelevel <- herb.dat.agg %>% 
  filter(system == "grassland") %>% 
  select(plot, Host_Name_std) %>% 
  mutate(ID = paste(plot, Host_Name_std))


################################################################################.
# ........ PLANTS ------------------------------- ##############################
################################################################################.

# ............ composition #####################################################
################################################################################.

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# functional groups ------------------------------------------------------------.

d_g_path_samplelevel <- d_g_path_samplelevel %>% 
  left_join(herb.dat.agg %>% 
              filter(system == "grassland") %>% 
              select(plot, Host_Name_std, pl.group),
            by = c("plot", "Host_Name_std")) %>%
  mutate(p.fgr.comp.grass = ifelse(pl.group == "grass", 1, 0),
         p.fgr.comp.legume = ifelse(pl.group == "legume", 1, 0)) %>% 
  select(-pl.group) %>% 
  mutate_at(vars(contains("p.fgr.comp")), ~as.factor(.)) 

# SLA / LDMC -------------------------------------------------------------------.

d_g_path_samplelevel <- d_p_traits_g %>% 
  group_by(plotid_withzero, Host_Name_std) %>% 
  summarise(p.fun.comp.ldmc = mean(ldmc, na.rm = T),
            p.fun.comp.sla = mean(sla, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(d_g_path_samplelevel, ., by = c(plot = "plotid_withzero", "Host_Name_std")) 

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_samplelevel <- d_g_path_samplelevel %>% 
  left_join(herb.dat.agg %>% 
              filter(system == "grassland") %>% 
              select(plot, Host_Name_std, cover),
            by = c("plot", "Host_Name_std")) %>% 
  rename(p.tax.abund = cover)

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# add plant height to exploratory species
d_p_height_g <-
  d_g_path_samplelevel %>% 
  select(Host_Name_std) %>% 
  distinct() %>% 
  left_join(d_p_height, by = "Host_Name_std") %>% 
  group_by(Host_Name_std) %>% 
  summarise(pl.height = mean(pl.height, na.rm = T)) %>%
  ungroup()

# filter missing cases (mainly Genus-level cases)
(missing <- d_p_height_g %>% 
  filter(is.na(pl.height)) %>% 
  select(Host_Name_std))

d_p_height_g <- d_p_height_g %>% 
  filter(!Host_Name_std %in% missing$Host_Name_std)

# deal with missing cases (use children entries of these):
missing_cases <- data.frame()
for (missing_i in missing$Host_Name_std){
  
  
  goon <- T
  children <- data.frame()
  while(goon){
    addon <- host.std %>% 
      filter(Host_IsChildTaxonOf == missing_i |
               Host_IsChildTaxonOf %in% children$Host_Name_std) %>% 
      filter(!Host_Name_original %in% children$Host_Name_original)
    
    children <- bind_rows(children, addon)
    
    if (nrow(addon) == 0) goon <- F
  }
  
  
  out <- d_p_height %>% 
    filter(Host_Name_std %in% children$Host_Name_std) %>% 
    summarise(pl.height = mean(pl.height))
  
  if (is.na(out$pl.height)) print(missing_i)
  
  missing_cases <- data.frame(Host_Name_std = missing_i,
                              pl.height = out$pl.height,
                              stringsAsFactors = F) %>% 
    bind_rows(missing_cases, .)
  
}

# add extracted data
d_p_height_g <- d_p_height_g %>% 
  bind_rows(missing_cases) %>% 
  arrange(Host_Name_std)

# calculate biomass from cover * height
d_g_path_samplelevel <- herb.dat.agg %>% 
  filter(system == "grassland") %>% 
  left_join(d_p_height_g, by = "Host_Name_std") %>% 
  mutate(p.fun.abund = cover * pl.height) %>% 
  select(plot, Host_Name_std, p.fun.abund) %>% 
  left_join(d_g_path_samplelevel, ., by = c("plot", "Host_Name_std"))


################################################################################.
# ........ INSECTS ------------------------------ ##############################
################################################################################.

# list of all plant species recorded
plant_list <- herb.dat.agg %>%  
  filter(system == "grassland") %>% 
  select(Host_Name_std) %>% 
  distinct() %>% 
  deframe() %>% 
  as.character()

# list of all insect species recorded plus rank
herb_list <- swnet_herb_agg %>% 
  select(Herb_Name_std, Herb_Rank, Herb_Genus, Herb_Species) %>% 
  distinct() %>% 
  mutate_all(~as.character(.))

herb_list <- herb_list %>% 
  mutate(Herb_Rank = paste0(substr(Herb_Rank, 1, 1), 
                            substr(tolower(Herb_Rank), 2, nchar(Herb_Rank))),
         Herb_Rank = paste0("Herb_", Herb_Rank))

# function returning for each plant species a list of herbivore species that were recorded in the explos
f_intlist <- function(pl_i, poly_all = T, poly_sub = T){
  int_target_i <- f_int_extract_plant(pl_i, poly_all = poly_all, poly_sub = poly_sub)
  
  out <- c()
  
  for (hr_i in 1:nrow(herb_list)){
    rank_target <- herb_list$Herb_Rank[hr_i]
    
    tmp <- int_target_i %>% 
      filter(!! sym(rank_target) == herb_list[hr_i, rank_target])
    
    if (nrow(tmp) > 0) out <- c(out, herb_list$Herb_Name_std[hr_i])
    
  }
  out
}

# apply the function:
l_pl_herbivores_g <- parallel::mclapply(plant_list, f_intlist)
names(l_pl_herbivores_g) <- plant_list


# create data frame with herbivore communities per plant and plot
swnet_agg_plant <- data.frame()
for (i in 1:nrow(d_g_path_samplelevel)){
  plt_i <- d_g_path_samplelevel$plot[i]
  pl_i <- d_g_path_samplelevel$Host_Name_std[i]
  
  swnet_agg_plant <- swnet_herb_agg %>% 
    filter(PlotID == plt_i,
           Herb_Name_std %in% l_pl_herbivores_g[[pl_i]]) %>% 
    select(PlotID, Herb_Name_std, NumberAdults) %>% 
    mutate(Host_Name_std = pl_i) %>% 
    bind_rows(swnet_agg_plant, .)
  
}

swnet_agg_plant <- swnet_agg_plant %>% 
  mutate(ID = paste(PlotID, Host_Name_std))

# prepare functional data ------------------------------------------------------.

sel.traits <- c("Body_Size", 
                "Dispersal_ability", 
                "Feeding_mode",  
                "Feeding_specialization", 
                "Stratum_use_short")

m.traits <- swnet_agg_plant %>% 
  select(Herb_Name_std) %>% 
  distinct() %>%
  left_join(traits_lh %>% 
              select(Herb_Name_std, !! sel.traits), by = "Herb_Name_std") %>%
  mutate(Dispersal_ability = as.numeric(Dispersal_ability)) %>% 
  mutate_at(vars(Feeding_mode,
                 Stratum_use_short),
            ~ as.factor(.)) %>% 
  mutate_at(vars(Feeding_specialization,
                 Dispersal_ability),
            ~ as.ordered(.)) %>% 
  column_to_rownames("Herb_Name_std")


swnet_ct_plant <- swnet_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID")

# exclude plots, where only one species was recorded (for functional diversity)
sum(rowSums(swnet_ct_plant > 0) <= 1)
swnet_ct_plant_no1 <- swnet_ct_plant[rowSums(swnet_ct_plant > 0) > 1, ] 
swnet_ct_plant_no1 <- swnet_ct_plant_no1[, colSums(swnet_ct_plant_no1) > 0]


# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

dissim.bc <- 
  swnet_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  vegdist()

# perform PCoA
pcoa.bc <- pco(dissim.bc, negvals = "zero")


d_g_path_samplelevel <- data.frame(ID = labels(dissim.bc), pcoa.bc$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(i.tax.comp.1 = X1,
         i.tax.comp.2 = X2,
         i.tax.comp.3 = X3) %>%
  left_join(d_g_path_samplelevel, ., by = "ID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# perform PCA to reduce trait dimensions
ord.traits <- ade4::dudi.mix(m.traits, scannf = F, nf = 4)

# extract components for species:
ord.traits <- ord.traits$li %>%
  mutate(Herb_Name_std = rownames(m.traits))

tpd.means <- ord.traits %>%
  mutate_at(vars(-Herb_Name_std), scale) %>% # scale 
  column_to_rownames("Herb_Name_std")

tps.sds <- matrix(0.5, nrow = nrow(tpd.means), ncol = ncol(tpd.means),
                  byrow = T) # use fixed SD values of 0.5


# following code takes a while
set.seed(213)
TPDs <- TPDsMean(species = rownames(tpd.means),
                 means = tpd.means,
                 sds = tps.sds) # calculate the TPDs

TPDc <- TPDc(TPDs = TPDs , sampUnit = swnet_ct_plant) # calculate the target TPDc

# following line takes several hours!
dissim.tpd.g  <- dissim(TPDc)$communities$dissimilarity # calculate community dissimilarity
rm(TPDs); rm(TPDc)

# perform PCoA
pcoa.tpd <- pco(as.dist(dissim.tpd.g), negvals = "zero")

d_g_path_samplelevel <- data.frame(ID = rownames(dissim.tpd.g), pcoa.tpd$vectors[, 1:3], 
                                   stringsAsFactors = F) %>%
  rename(i.fun.comp.1 = X1,
         i.fun.comp.2 = X2,
         i.fun.comp.3 = X3) %>%
  left_join(d_g_path_samplelevel, ., by = "ID")


# what do the traits mean?
cwms.i.g.sample <- functcomp(m.traits, as.matrix(swnet_ct_plant)[, rownames(m.traits)])
cwms.i.g.sample %>%
  rownames_to_column("ID") %>%
  left_join(d_g_path_samplelevel %>%
              select(i.fun.comp.1, i.fun.comp.2, i.fun.comp.3, ID),
            by = "ID") %>%
  select(-ID) %>%
  pairs()


# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_g_path_samplelevel <- swnet_agg_plant %>% 
  group_by(ID) %>% 
  summarise(i.tax.abund = sum(NumberAdults)) %>% 
  ungroup() %>% 
  left_join(d_g_path_samplelevel, ., by = "ID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# extract allometric regression coefficient for herbivores:
f_lw <- function(herb_i, data){
  genus_i <- as.character(unique(swnet_herb_agg$Herb_Genus[swnet_herb_agg$Herb_Name_std == herb_i]))
  family_i <- as.character(unique(swnet_herb_agg$Herb_Family[swnet_herb_agg$Herb_Name_std == herb_i]))
  suborder_i <- as.character(unique(swnet_herb_agg$Herb_Suborder[swnet_herb_agg$Herb_Name_std == herb_i]))
  order_i <- as.character(unique(swnet_herb_agg$Herb_Order[swnet_herb_agg$Herb_Name_std == herb_i]))
  
  
  if (genus_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == genus_i],
                      a = data$a[data$Herb_Name_std == genus_i],
                      stringsAsFactors = F)
    out
  } else if (family_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == family_i],
                      a = data$a[data$Herb_Name_std == family_i],
                      stringsAsFactors = F)
    out
  } else if (suborder_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == suborder_i],
                      a = data$a[data$Herb_Name_std == suborder_i],
                      stringsAsFactors = F)
    out
  } else if (order_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == order_i],
                      a = data$a[data$Herb_Name_std == order_i],
                      stringsAsFactors = F)
    out
  }
}

# apply function
regs_swnet <- unique(swnet_agg_plant$Herb_Name_std) %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

regs_swnet <- regs_swnet %>% 
  left_join(select(traits_lh, Herb_Name_std, Body_Size), by = "Herb_Name_std") 

regs_swnet <- regs_swnet %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)

d_g_path_samplelevel <- swnet_agg_plant %>% 
  left_join(regs_swnet, by = "Herb_Name_std") %>% 
  mutate(metabolic_rel_sum = metabolic_rel * NumberAdults) %>% 
  group_by(ID) %>% 
  summarise(i.fun.abund = sum(metabolic_rel_sum, na.rm = T)) %>% 
  ungroup() %>%
  left_join(d_g_path_samplelevel, ., by = "ID")

# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# q = 0
d_g_path_samplelevel <- swnet_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(i.tax.div.0 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_g_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.0 = ifelse(is.na(i.tax.div.0), 0, i.tax.div.0))

# q = 1
d_g_path_samplelevel <- swnet_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(i.tax.div.1 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_g_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.1 = ifelse(is.na(i.tax.div.1), 0, i.tax.div.1))

# q = 2
d_g_path_samplelevel <- swnet_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(i.tax.div.2 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_g_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.2 = ifelse(is.na(i.tax.div.2), 0, i.tax.div.2))


# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

dij <- m.traits %>% 
  daisy(metric = "gower") %>% 
  as.matrix()

g.i.fund.plant <- FunD(data = as.data.frame(t(swnet_ct_plant_no1)), 
                       dij = dij[names(swnet_ct_plant_no1), names(swnet_ct_plant_no1)], 
                       tau = max(dij),
                       q = seq(0, 2, by = 1), 
                       boot = 1, 
                       datatype = "abundance")$fortau

d_g_path_samplelevel <- g.i.fund.plant %>% 
  rename(ID = site) %>% 
  select(estimate, q, ID) %>% 
  spread(q, estimate) %>% 
  rename(i.fun.div.0 = `0`,
         i.fun.div.1 = `1`,
         i.fun.div.2 = `2`) %>% 
  left_join(d_g_path_samplelevel, ., by = "ID")





################################################################################.
# ........ FINALIZE FOR ANALYSIS ------------------------------ ################
################################################################################.

# corrections of NA values
d_g_path_samplelevel <-
  d_g_path_samplelevel %>% 
  mutate(i.tax.abund = ifelse(is.na(i.tax.abund), 0, i.tax.abund),
         i.fun.abund = ifelse(is.na(i.fun.abund), 0, i.fun.abund),
         # set functional diversity values for 0 or 1 species to 1, which is the minimum possible value
         i.fun.div.0 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.0),
         i.fun.div.1 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.1),
         i.fun.div.2 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.2),
         # set composition to zero and add extra factor for case of no herbivores 
         i.tax.comp.1 = ifelse(i.tax.abund == 0, 0, i.tax.comp.1),
         i.tax.comp.2 = ifelse(i.tax.abund == 0, 0, i.tax.comp.2),
         i.tax.comp.3 = ifelse(i.tax.abund == 0, 0, i.tax.comp.3),
         i.fun.comp.1 = ifelse(i.fun.abund == 0, 0, i.fun.comp.1),
         i.fun.comp.2 = ifelse(i.fun.abund == 0, 0, i.fun.comp.2),
         i.fun.comp.3 = ifelse(i.fun.abund == 0, 0, i.fun.comp.3),
         i.tax.comp.iszero = factor(ifelse(i.tax.abund == 0, "yes", "no"))) 

# scale variables, add herbivory and region data and scaled LUI data
d_g_path_samplelevel_z <- d_g_path_samplelevel %>% 
  mutate_if(is.numeric, ~ (. - mean(., na.rm = T)) / (2 * sd(.,  na.rm = T)))

# remove rows containing NAs
d_g_path_samplelevel_z <- d_g_path_samplelevel_z %>% 
  filter_all(~!is.na(.)) 
  
  
# .................................................................. #################################################
# .................................................................. #################################################
# .................................................................. #################################################
# .................................................................. #################################################
# .................................................................. #################################################
# .................................................................. #################################################


################################################################################.
################################################################################.
# ////////// ANALYSIS FOREST \\\\\\\\\\\\ ######################################
################################################################################.
################################################################################.


################################################################################.
#  COMMUNITY LEVEL #############################################################
################################################################################.

d_f_path_plotlevel <- data.frame(PlotID = unique(herb.dat.comm$plot),
                                 stringsAsFactors = F) %>% 
  filter(grepl("W", PlotID))

################################################################################.
# ........ PLANTS ------------------------------- ##############################
################################################################################.

# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# dissimilarity based on plant species composition

dissim.bc.p.plot <-
  plant.cover %>% 
  filter(grepl("W", PlotID)) %>% 
  select(PlotID, Host_Name_std, cover) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  vegdist()

# perform PCoA
pcoa.bc.p.plot <- pco(dissim.bc.p.plot, negvals = "zero")

d_f_path_plotlevel <- data.frame(PlotID = labels(dissim.bc.p.plot), 
                                 pcoa.bc.p.plot$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(p.tax.comp.1 = X1,
         p.tax.comp.2 = X2,
         p.tax.comp.3 = X3) %>%
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# functional groups ------------------------------------------------------------.

d_f_path_plotlevel <- herb.dat.agg %>% 
  filter(system == "forest") %>% 
  group_by(plot) %>% 
  summarise(p.fgr.comp.grass  = sum(prop.cover * (pl.group == "grass")),
            p.fgr.comp.forb = sum(prop.cover * (pl.group == "forb")),
            p.fgr.comp.tree = sum(prop.cover * (pl.group == "tree")),
            p.fgr.comp.geophyt = sum(prop.cover * (pl.group == "geophyt"))) %>% 
  ungroup() %>%
  left_join(d_f_path_plotlevel, ., by = c("PlotID" = "plot"))


# SLA / LDMC -------------------------------------------------------------------.

plant_cover_f_ct <- herb.dat.agg %>% 
  filter(system == "forest") %>% 
  mutate(ID = paste(plot, Host_Name_std)) %>% 
  select(ID, cover, plot) %>% 
  pivot_wider(names_from = ID, 
              values_from = cover, 
              values_fill = list(cover = 0)) %>% 
  column_to_rownames("plot")

d_p_traits_f_mean <- d_p_traits_f %>% 
  mutate(ID = paste(plotid_withzero, Host_Name_std)) %>% 
  group_by(ID) %>% 
  summarise_at(vars(ldmc, sla), ~mean(., na.rm = T)) %>% 
  ungroup() %>% 
  filter(ID %in% names(plant_cover_f_ct)) %>% 
  column_to_rownames("ID")

CWM_plant_f <- functcomp(d_p_traits_f_mean, as.matrix(plant_cover_f_ct[, rownames(d_p_traits_f_mean)])) %>% 
  rownames_to_column("PlotID")

d_f_path_plotlevel <- CWM_plant_f %>%
  rename(p.fun.comp.ldmc = ldmc,
         p.fun.comp.sla = sla) %>%
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_plotlevel <- plant.cover %>% 
  filter(grepl("W", PlotID)) %>% 
  group_by(PlotID) %>% 
  summarise(p.tax.abund = sum(cover)) %>% 
  ungroup() %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# approach one: use LAI

d_f_path_plotlevel <- d_lai %>% 
  rename(p.fun.abund.lai = LAI) %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# approach two: use estimated biomass

d_f_path_plotlevel <- d_bm_f %>% 
  filter(PlotID != "HEW51") %>% # incomplete data (no trees!) for this plot
  group_by(PlotID) %>% 
  summarise(p.fun.abund.bm = sum(biomass)) %>% 
  ungroup() %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_plotlevel <- plant.cover %>% 
  filter(grepl("W", PlotID)) %>% 
  select(PlotID, Host_Name_std, cover) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(p.tax.div.0 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

d_f_path_plotlevel <- plant.cover %>% 
  filter(grepl("W", PlotID)) %>% 
  select(PlotID, Host_Name_std, cover) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(p.tax.div.1 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

d_f_path_plotlevel <- plant.cover %>% 
  filter(grepl("W", PlotID)) %>% 
  select(PlotID, Host_Name_std, cover) %>% 
  spread(Host_Name_std, cover, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(p.tax.div.2 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")


# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

m_p_traits_sample <- d_p_traits_f %>% 
  group_by(plotid_withzero, Host_Name_std) %>% 
  summarise(ldmc = mean(ldmc, na.rm = T),
            sla = mean(sla, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(ID = paste(plotid_withzero, Host_Name_std)) %>% 
  column_to_rownames("ID") %>% 
  select(-c(plotid_withzero, Host_Name_std))

dij_p <- m_p_traits_sample %>% 
  daisy(metric = "euclidean", stand = T) %>% 
  as.matrix()

f.p.fund <-
  herb.dat.agg %>%
  filter(system == "forest") %>% 
  mutate(ID = paste(plot, Host_Name_std)) %>% 
  select(plot, ID, cover) %>% 
  filter(ID %in% labels(dij_p)[[1]]) %>% # exclude one NA
  spread(plot, cover, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  {FunD(data = .,
        dij = dij_p[rownames(.), rownames(.)],
        tau = max(dij_p),
        q = seq(0, 2, by = 1), 
        boot = 0, 
        datatype = "abundance")$fortau}


d_f_path_plotlevel <- f.p.fund %>% 
  rename(PlotID = site) %>% 
  select(estimate, q, PlotID) %>% 
  spread(q, estimate) %>% 
  rename(p.fun.div.0 = `0`,
         p.fun.div.1 = `1`,
         p.fun.div.2 = `2`) %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")




################################################################################.
# ........ INSECTS ------------------------------ ##############################
################################################################################.

# restrict to leaf-feeding herbivores
wintr_herb_agg <- wintr_agg %>% 
  filter(Herb_Rank == "SPECIES" & Herb_Species %in% int.db.leaves$Herb_Species |
           Herb_Rank == "GENUS" & Herb_Genus %in% int.db.leaves$Herb_Genus) 
  
# make cross table -------------------------------------------------------------.

m.wintr.ct <- wintr_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID")

# exclude plots, where only one species was recorded (no effect)
m.wintr.ct <- m.wintr.ct[rowSums(m.wintr.ct > 0) > 1, ]
m.wintr.ct <- m.wintr.ct[, colSums(m.wintr.ct) > 0]

# prepare trait data -----------------------------------------------------------.

sel.traits <- c("Body_Size", 
                "Dispersal_ability", 
                "Feeding_mode",
                "Feeding_specialization", 
                "Stratum_use_short")


m.traits <-
  wintr_herb_agg %>% 
  select(Herb_Name_std) %>% 
  distinct() %>%
  left_join(traits_lh %>%
              select(Herb_Name_std, !! sel.traits), by = "Herb_Name_std") %>%
  mutate(Dispersal_ability = as.numeric(Dispersal_ability)) %>% 
  mutate_at(vars(Feeding_mode,
                 Stratum_use_short),
            ~ as.factor(.)) %>% 
  mutate_at(vars(Feeding_specialization,
                 Dispersal_ability),
            ~ as.ordered(.)) %>% 
  column_to_rownames("Herb_Name_std")


# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

dissim.bc.i.plot <- 
  wintr_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  vegdist()

# perform PCoA
pcoa.bc.plot <- pco(dissim.bc.i.plot, negvals = "zero")


d_f_path_plotlevel <- data.frame(PlotID = labels(dissim.bc.i.plot), 
                                 pcoa.bc.plot$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(i.tax.comp.1 = X1,
         i.tax.comp.2 = X2,
         i.tax.comp.3 = X3) %>%
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# perform PCA to reduce trait dimensions
ord.traits <- ade4::dudi.mix(m.traits, scannf = F, nf = 4)

# extract components for species:
ord.traits <- ord.traits$li %>%
  mutate(Herb_Name_std = rownames(m.traits))

tpd.means <- ord.traits %>%
  mutate_at(vars(-Herb_Name_std), scale) %>% # scale 
  column_to_rownames("Herb_Name_std")

tps.sds <- matrix(0.5, nrow = nrow(tpd.means), ncol = ncol(tpd.means),
                  byrow = T) # use fixed SD values of 0.5

set.seed(61)
TPDs <- TPDsMean(species = rownames(tpd.means),
                 means = tpd.means,
                 sds = tps.sds) # calculate the TPDs

TPDc <- TPDc(TPDs = TPDs , sampUnit = m.wintr.ct) # calculate the target TPDc

dissim.tpd.plot  <- dissim(TPDc)$communities$dissimilarity # calculate community dissimilarity

rm(TPDs); rm(TPDc)

# perform PCoA
pcoa.tpd <- pco(as.dist(dissim.tpd.plot), negvals = "zero")


d_f_path_plotlevel <- data.frame(PlotID = rownames(dissim.tpd.plot), 
                                 pcoa.tpd$vectors[, 1:3], 
                                 stringsAsFactors = F) %>%
  rename(i.fun.comp.1 = X1,
         i.fun.comp.2 = X2,
         i.fun.comp.3 = X3) %>%
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_plotlevel <- herbivore.abundance %>% 
  rename(i.tax.abund.sort = herbivore.abundance) %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

d_f_path_plotlevel <- wintr_herb_agg %>% 
  group_by(PlotID) %>% 
  summarise(i.tax.abund.det = sum(NumberAdults)) %>%
  ungroup() %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# ....................functional -----------------------------------------------
# ------------------------------------------------------------------------------.

# extract allometric regression coefficient for herbivores:
f_lw <- function(herb_i, data){
  genus_i <- as.character(unique(wintr_herb_agg$Herb_Genus[wintr_herb_agg$Herb_Name_std == herb_i]))
  family_i <- as.character(unique(wintr_herb_agg$Herb_Family[wintr_herb_agg$Herb_Name_std == herb_i]))
  suborder_i <- as.character(unique(wintr_herb_agg$Herb_Suborder[wintr_herb_agg$Herb_Name_std == herb_i]))
  order_i <- as.character(unique(wintr_herb_agg$Herb_Order[wintr_herb_agg$Herb_Name_std == herb_i]))
  
  if (genus_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == genus_i],
                      a = data$a[data$Herb_Name_std == genus_i],
                      stringsAsFactors = F)
    out
  } else if (family_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == family_i],
                      a = data$a[data$Herb_Name_std == family_i],
                      stringsAsFactors = F)
    out
  } else if (suborder_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == suborder_i],
                      a = data$a[data$Herb_Name_std == suborder_i],
                      stringsAsFactors = F)
    out
  } else if (order_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == order_i],
                      a = data$a[data$Herb_Name_std == order_i],
                      stringsAsFactors = F)
    out
  }
}

# apply function
regs_wintr <- unique(wintr_herb_agg$Herb_Name_std) %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

regs_wintr <- regs_wintr %>% 
  left_join(select(traits_lh, Herb_Name_std, Body_Size), by = "Herb_Name_std") 


regs_wintr <- regs_wintr %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)

d_f_path_plotlevel <- wintr_herb_agg %>% 
  left_join(regs_wintr, by = "Herb_Name_std") %>% 
  mutate(metabolic_rel_sum = metabolic_rel * NumberAdults) %>% 
  group_by(PlotID) %>% 
  summarise(i.fun.abund = sum(metabolic_rel_sum, na.rm = T)) %>%
  ungroup() %>%
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# q = 0
d_f_path_plotlevel <- wintr_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(i.tax.div.0 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# q = 1
d_f_path_plotlevel <- wintr_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(i.tax.div.1 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# q = 2
d_f_path_plotlevel <- wintr_herb_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("PlotID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(i.tax.div.2 = .) %>% 
  rownames_to_column("PlotID") %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

dij <- m.traits %>% 
  daisy(metric = "gower") %>% 
  as.matrix()

f.i.fund <- FunD(data = as.data.frame(t(m.wintr.ct)), 
                 dij = dij[names(m.wintr.ct), names(m.wintr.ct)], 
                 tau = max(dij),
                 q = seq(0, 2, by = 1), 
                 boot = 1, 
                 datatype = "abundance")$fortau

d_f_path_plotlevel <- f.i.fund %>% 
  rename(PlotID = site) %>% 
  select(estimate, q, PlotID) %>% 
  spread(q, estimate) %>% 
  rename(i.fun.div.0 = `0`,
         i.fun.div.1 = `1`,
         i.fun.div.2 = `2`) %>% 
  left_join(d_f_path_plotlevel, ., by = "PlotID")



################################################################################.
# ........ FINALIZE FOR ANALYSIS ------------------------------ ################
################################################################################.

# add land-use data, standardize data
d_f_path_plotlevel_z <- d_f_path_plotlevel %>%
  filter_all(~!is.na(.)) %>%  # remove all entries with NAs in any of the variables
  mutate_if(is.numeric, ~ (. - mean(., na.rm = T)) / (2 * sd(.,  na.rm = T))) # scale

# .................................................................. #################################################
# .................................................................. #################################################

################################################################################.
# SPECIES LEVEL ################################################################
################################################################################.

# initiate data frames ---------------------------------------------------------.
d_f_path_samplelevel <- herb.dat.agg %>% 
  filter(system == "forest") %>% 
  select(plot, Host_Name_std) %>% 
  mutate(ID = paste(plot, Host_Name_std))


################################################################################.
# ........ PLANTS ------------------------------- ##############################
################################################################################.

# ............ composition #####################################################
################################################################################.

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# functional groups ------------------------------------------------------------.

d_f_path_samplelevel <-
  d_f_path_samplelevel %>% 
  left_join(herb.dat.agg %>%
              filter(system == "forest") %>% 
              select(plot, Host_Name_std, pl.group),
            by = c("plot", "Host_Name_std")) %>%
  mutate(p.fgr.comp.grass = ifelse(pl.group == "grass", 1, 0),
         p.fgr.comp.forb = ifelse(pl.group == "forb", 1, 0),
         p.fgr.comp.tree = ifelse(pl.group == "tree", 1, 0),
         p.fgr.comp.geophyt = ifelse(pl.group == "geophyt", 1, 0)) %>%
  select(-pl.group) %>% 
  mutate_at(vars(contains("p.fgr.comp")), ~as.factor(.)) 

# SLA / LDMC -------------------------------------------------------------------.

d_f_path_samplelevel <-
  d_p_traits_f %>% 
  group_by(plotid_withzero, Host_Name_std) %>% 
  summarise(p.fun.comp.ldmc = mean(ldmc, na.rm = T),
            p.fun.comp.sla = mean(sla, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(d_f_path_samplelevel, ., by = c(plot = "plotid_withzero", "Host_Name_std"))


# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_samplelevel <- d_f_path_samplelevel %>%  
  left_join(herb.dat.agg %>% 
              filter(system == "forest") %>% 
              select(plot, Host_Name_std, cover),
            by = c("plot", "Host_Name_std")) %>% 
  rename(p.tax.abund = cover)

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_samplelevel <- d_bm_f %>% 
  group_by(PlotID, Name_std) %>% 
  summarise(p.fun.abund = sum(biomass)) %>% 
  ungroup() %>% 
  left_join(d_f_path_samplelevel, ., by = c(plot = "PlotID", Host_Name_std = "Name_std")) %>% 
  mutate(p.fun.abund = ifelse(p.fgr.comp.tree == "1" & plot == "HEW51", NA, p.fun.abund),
         p.fun.abund = ifelse(is.na(p.fun.abund) & plot != "HEW51", 0, p.fun.abund)) # treelayer missing in HEW51, thus leaves NAs there


################################################################################.
# ........ INSECTS ------------------------------ ##############################
################################################################################.

# list of all plant species recorded
plant_list <- herb.dat.agg %>%  
  filter(system == "forest") %>% 
  select(Host_Name_std) %>% 
  distinct() %>% 
  deframe()

# list of all insect species recorded plus rank
herb_list <- wintr_herb_agg %>% 
  select(Herb_Name_std, Herb_Rank, Herb_Genus, Herb_Species) %>% 
  distinct()

herb_list <- herb_list %>% 
  mutate(Herb_Rank = paste0(substr(Herb_Rank, 1, 1), 
                            substr(tolower(Herb_Rank), 2, nchar(Herb_Rank))),
         Herb_Rank = paste0("Herb_", Herb_Rank))

# function returning for each plant species a list of herbivore species that were recorded in the explos
f_intlist <- function(pl_i, poly_all = T, poly_sub = T){
  int_target_i <- f_int_extract_plant(pl_i, poly_all = poly_all, poly_sub = poly_sub)

  out <- c()
  
  for (hr_i in 1:nrow(herb_list)){
    rank_target <- herb_list$Herb_Rank[hr_i]
    
    tmp <- int_target_i %>% 
      filter(!! sym(rank_target) == herb_list[hr_i, rank_target])
    
    if (nrow(tmp) > 0) out <- c(out, herb_list$Herb_Name_std[hr_i])
    
  }
  out
}


# apply the function:
l_pl_herbivores_f <- parallel::mclapply(plant_list, f_intlist)
names(l_pl_herbivores_f) <- plant_list


# create data frame with herbivore communities per plant and plot
wintr_agg_plant <- data.frame()
for (i in 1:nrow(d_f_path_samplelevel)){
  plt_i <- d_f_path_samplelevel$plot[i]
  pl_i <- d_f_path_samplelevel$Host_Name_std[i]
  
  wintr_agg_plant <- wintr_herb_agg %>% 
    filter(PlotID == plt_i,
           Herb_Name_std %in% l_pl_herbivores_f[[pl_i]]) %>% 
    select(PlotID, Herb_Name_std, NumberAdults) %>% 
    mutate(Host_Name_std = pl_i) %>% 
    bind_rows(wintr_agg_plant, .)
  
}

wintr_agg_plant <- wintr_agg_plant %>% 
  mutate(ID = paste(PlotID, Host_Name_std))

# prepare functional data ------------------------------------------------------.

sel.traits <- c("Body_Size", 
                "Dispersal_ability", 
                "Feeding_mode",  
                "Feeding_specialization", 
                "Stratum_use_short")

m.traits <- wintr_agg_plant %>% 
  select(Herb_Name_std) %>% 
  distinct() %>%
  left_join(traits_lh %>%
              select(Herb_Name_std, !! sel.traits), by = "Herb_Name_std") %>%
  mutate(Dispersal_ability = as.numeric(Dispersal_ability)) %>% 
  mutate_at(vars(Feeding_mode,
                 Stratum_use_short),
            ~ as.factor(.)) %>% 
  mutate_at(vars(Feeding_specialization,
                 Dispersal_ability),
            ~ as.ordered(.)) %>% 
  column_to_rownames("Herb_Name_std")


wintr_ct_plant <- wintr_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>%
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID")

# exclude plots, where only one species was recorded (for functional diversity)
sum(rowSums(wintr_ct_plant > 0) <= 1)
wintr_ct_plantt_no1 <- wintr_ct_plant[rowSums(wintr_ct_plant > 0) > 1, ] 
wintr_ct_plantt_no1 <- wintr_ct_plantt_no1[, colSums(wintr_ct_plantt_no1) > 0]


# ............ composition #####################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

dissim.bc <- 
  wintr_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  vegdist()

# perform PCoA
pcoa.bc <- pco(dissim.bc, negvals = "zero")

d_f_path_samplelevel <- data.frame(ID = labels(dissim.bc), pcoa.bc$vectors[, 1:3], stringsAsFactors = F) %>%
  rename(i.tax.comp.1 = X1,
         i.tax.comp.2 = X2,
         i.tax.comp.3 = X3) %>%
  left_join(d_f_path_samplelevel, ., by = "ID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# perform PCA to reduce trait dimensions
ord.traits <- ade4::dudi.mix(m.traits, scannf = F, nf = 4)

# extract components for species:
ord.traits <- ord.traits$li %>%
  mutate(Herb_Name_std = rownames(m.traits))

tpd.means <- ord.traits %>%
  mutate_at(vars(-Herb_Name_std), scale) %>% # scale 
  column_to_rownames("Herb_Name_std")

tps.sds <- matrix(0.5, nrow = nrow(tpd.means), ncol = ncol(tpd.means),
                  byrow = T) # use fixed SD values of 0.5

set.seed(24)
TPDs <- TPDsMean(species = rownames(tpd.means),
                 means = tpd.means,
                 sds = tps.sds) # calculate the TPDs

TPDc <- TPDc(TPDs = TPDs , sampUnit = wintr_ct_plant) # calculate the target TPDc

# following line takes several hours!
dissim.tpd.f  <- dissim(TPDc)$communities$dissimilarity # calculate community dissimilarity
rm(TPDs); rm(TPDc)

# perform PCoA
pcoa.tpd <- pco(as.dist(dissim.tpd.f), negvals = "zero")

d_f_path_samplelevel <- data.frame(ID = rownames(dissim.tpd.f), pcoa.tpd$vectors[, 1:3], 
                                   stringsAsFactors = F) %>%
  rename(i.fun.comp.1 = X1,
         i.fun.comp.2 = X2,
         i.fun.comp.3 = X3) %>%
  left_join(d_f_path_samplelevel, ., by = "ID")

# ............ abundance #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

d_f_path_samplelevel <- wintr_agg_plant %>% 
  group_by(ID) %>% 
  summarise(i.tax.abund = sum(NumberAdults)) %>% 
  ungroup() %>% 
  left_join(d_f_path_samplelevel, ., by = "ID")

# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

# extract allometric regression coefficient for herbivores:
f_lw <- function(herb_i, data){
  genus_i <- as.character(unique(wintr_herb_agg$Herb_Genus[wintr_herb_agg$Herb_Name_std == herb_i]))
  family_i <- as.character(unique(wintr_herb_agg$Herb_Family[wintr_herb_agg$Herb_Name_std == herb_i]))
  suborder_i <- as.character(unique(wintr_herb_agg$Herb_Suborder[wintr_herb_agg$Herb_Name_std == herb_i]))
  order_i <- as.character(unique(wintr_herb_agg$Herb_Order[wintr_herb_agg$Herb_Name_std == herb_i]))
  
  if (genus_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == genus_i],
                      a = data$a[data$Herb_Name_std == genus_i],
                      stringsAsFactors = F)
    out
  } else if (family_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == family_i],
                      a = data$a[data$Herb_Name_std == family_i],
                      stringsAsFactors = F)
    out
  } else if (suborder_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == suborder_i],
                      a = data$a[data$Herb_Name_std == suborder_i],
                      stringsAsFactors = F)
    out
  } else if (order_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == order_i],
                      a = data$a[data$Herb_Name_std == order_i],
                      stringsAsFactors = F)
    out
  }
}

# apply function
regs_wintr <- unique(wintr_agg_plant$Herb_Name_std) %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

regs_wintr <- regs_wintr %>% 
  left_join(select(traits_lh, Herb_Name_std, Body_Size), by = "Herb_Name_std") 

regs_wintr <- regs_wintr %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)

d_f_path_samplelevel <- wintr_agg_plant %>% 
  left_join(regs_wintr, by = "Herb_Name_std") %>% 
  mutate(metabolic_rel_sum = metabolic_rel * NumberAdults) %>% 
  group_by(ID) %>% 
  summarise(i.fun.abund = sum(metabolic_rel_sum)) %>% 
  ungroup() %>%
  left_join(d_f_path_samplelevel, ., by = "ID")

# ............ diversity #######################################################
################################################################################.

# .................... taxonomic -----------------------------------------------
# ------------------------------------------------------------------------------.

# q = 0
d_f_path_samplelevel <- wintr_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 0) %>% 
  data.frame(i.tax.div.0 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_f_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.0 = ifelse(is.na(i.tax.div.0), 0, i.tax.div.0))

# q = 1
d_f_path_samplelevel <- wintr_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 1) %>% 
  data.frame(i.tax.div.1 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_f_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.1 = ifelse(is.na(i.tax.div.1), 0, i.tax.div.1))

# q = 2
d_f_path_samplelevel <- wintr_agg_plant %>% 
  select(ID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  column_to_rownames("ID") %>% 
  hill_taxa(q = 2) %>% 
  data.frame(i.tax.div.2 = .) %>% 
  rownames_to_column("ID") %>% 
  left_join(d_f_path_samplelevel, ., by = "ID") %>% 
  mutate(i.tax.div.2 = ifelse(is.na(i.tax.div.2), 0, i.tax.div.2))


# .................... functional ----------------------------------------------
# ------------------------------------------------------------------------------.

dij <- m.traits %>% 
  daisy(metric = "gower") %>% 
  as.matrix()

g.i.fund.plant <- FunD(data = as.data.frame(t(wintr_ct_plantt_no1)), 
                       dij = dij[names(wintr_ct_plantt_no1), names(wintr_ct_plantt_no1)], 
                       tau = max(dij),
                       q = seq(0, 2, by = 1), 
                       boot = 1, 
                       datatype = "abundance")$fortau

d_f_path_samplelevel <- g.i.fund.plant %>% 
  rename(ID = site) %>% 
  select(estimate, q, ID) %>% 
  spread(q, estimate) %>% 
  rename(i.fun.div.0 = `0`,
         i.fun.div.1 = `1`,
         i.fun.div.2 = `2`) %>% 
  left_join(d_f_path_samplelevel, ., by = "ID")




################################################################################.
# ........ FINALIZE FOR ANALYSIS ------------------------------ ################
################################################################################.

# corrections of NA values
d_f_path_samplelevel <-
  d_f_path_samplelevel %>% 
  mutate(i.tax.abund = ifelse(is.na(i.tax.abund), 0, i.tax.abund),
         i.fun.abund = ifelse(is.na(i.fun.abund), 0, i.fun.abund),
         # set functional diversity values for 0 or 1 species to 1, which is the minimum possible value
         i.fun.div.0 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.0),
         i.fun.div.1 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.1),
         i.fun.div.2 = ifelse(i.tax.div.0 %in% c(0, 1), 1, i.fun.div.2),
         # set composition to zero and add extra factor for case of no herbivores 
         i.tax.comp.1 = ifelse(i.tax.abund == 0, 0, i.tax.comp.1),
         i.tax.comp.2 = ifelse(i.tax.abund == 0, 0, i.tax.comp.2),
         i.tax.comp.3 = ifelse(i.tax.abund == 0, 0, i.tax.comp.3),
         i.fun.comp.1 = ifelse(i.fun.abund == 0, 0, i.fun.comp.1),
         i.fun.comp.2 = ifelse(i.fun.abund == 0, 0, i.fun.comp.2),
         i.fun.comp.3 = ifelse(i.fun.abund == 0, 0, i.fun.comp.3),
         i.tax.comp.iszero = factor(ifelse(i.tax.abund == 0, "yes", "no"))) 


# scale variables, add herbivory, region data and scaled ForMI data
d_f_path_samplelevel_z <- d_f_path_samplelevel %>% 
  mutate_if(is.numeric, ~ (. - mean(., na.rm = T)) / (2 * sd(.,  na.rm = T))) 

# remove rows with NAs
d_f_path_samplelevel_z <- d_f_path_samplelevel_z %>% 
  filter_all(~!is.na(.)) # remove all entries with NAs in any of the variables


################################################################################.
################################################################################.
# ///////// EXPORT COMBINED DATA \\\\\\\\\\ ####################################
################################################################################.
################################################################################.

herb.dat.comm %>% 
  full_join(d_f_path_plotlevel_z %>% 
              bind_rows(d_g_path_plotlevel_z) %>% 
            by = c("plot" = "PlotID")) %>% 
  select(system, region, plot, starts_with("a."), 
         LUI, G_std, M_std, F_std, 
         ForMI, Inonat, Iharv, Idwcut, 
         RW, HW,
         p.tax.comp.1, p.tax.comp.2, p.tax.comp.3,
         p.fgr.comp.grass, p.fgr.comp.forb, p.fgr.comp.tree,
         p.fgr.comp.geophyt, p.fgr.comp.legume, 
         p.fun.comp.ldmc, p.fun.comp.sla,
         p.fun.comp.lignin, p.fun.comp.prim.fiber,
         p.fun.comp.P, p.fun.comp.N, p.fun.comp.Ca, p.fun.comp.Mg,
         p.tax.abund, p.fun.abund,
         p.fun.abund.lai, p.fun.abund.bm, 
         p.tax.div.0, p.tax.div.1, p.tax.div.2, 
         p.fun.div.0, p.fun.div.1, p.fun.div.2, 
         i.tax.comp.1, i.tax.comp.2, i.tax.comp.3,
         i.fun.comp.1, i.fun.comp.2, i.fun.comp.3,
         i.tax.abund.sort, i.tax.abund.det,
         everything())  %>% 
  write.table("Data/Dat_communitylevel.txt")

herb.dat.agg %>% 
  full_join(d_f_path_samplelevel_z %>% 
              bind_rows(d_g_path_samplelevel_z),
            by = c("plot" = "PlotID", "Host_Name_std")) %>% 
  select(system, region, ID, plot, Host_Name_std, pl.species.abbr, 
         starts_with("a."),
         pl.group, cover, prop.cover,
         LUI, G_std, M_std, F_std, 
         ForMI, Inonat, Iharv, Idwcut, 
         RW, HW,
         p.fgr.comp.grass, p.fgr.comp.forb, p.fgr.comp.tree,
         p.fgr.comp.geophyt, p.fgr.comp.legume, 
         p.fun.comp.ldmc, p.fun.comp.sla,
         p.tax.abund, p.fun.abund,
         i.tax.comp.iszero,
         everything()) %>% 
  write.table("Data/Dat_specieslevel.txt")