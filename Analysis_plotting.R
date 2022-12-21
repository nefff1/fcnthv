################################################################################.
# SET UP -------------- ########################################################
################################################################################.

# load packages ----------------------------------------------------------------

library(tidyverse); options(dplyr.summarise.inform = FALSE); theme_set(theme_classic())
library(glmmTMB)
library(sp)
library(gstat)
library(piecewiseSEM)
library(parallel)
library(cowplot)

# standardisation function
f_std <- function(x) {
  if (is.numeric(x)){
    (x - mean(x, na.rm = T)) / (2 * sd(x,  na.rm = T))
  } else {
    x
  }
}

# define metric names for renaming ---------------------------------------------

# community level
metrics_comm <- c(system = "System",
                    region = "Region",
                    plot = "Plot",
                    a.herb2.prop.comm = "Herb_rate_total",
                    a.chewing.prop.comm = "Herb_rate_chewing",
                    a.scraping.prop.comm = "Herb_rate_scraping",
                    a.sucking.prop.comm = "Herb_rate_sucking",
                    a.mines.prop.comm = "Herb_rate_mines",
                    a.gallstot.prop.comm = "Herb_rate_galls",
                    ForMI = "ForMI",
                    Inonat = "Inonat",
                    Idwcut = "Idwcut",
                    Iharv = "Iharv",
                    LUI = "LUI",
                    F_std = "Fertilization",
                    M_std = "Mowing",
                    G_std = "Grazing",
                    p.tax.comp.1 = "P_Com_Tax_axis1",
                    p.tax.comp.2 = "P_Com_Tax_axis2",
                    p.tax.comp.3 = "P_Com_Tax_axis3",
                    p.fgr.comp.grass = "P_Com_Fun_grass",
                    p.fgr.comp.legume = "P_Com_Fun_legume",
                    p.fgr.comp.forb = "P_Com_Fun_forb",
                    p.fgr.comp.tree = "P_Com_Fun_tree",
                    p.fgr.comp.geophyt = "P_Com_Fun_geophyte",
                    p.fun.comp.ldmc = "P_Com_Fun_LDMC",
                    p.fun.comp.sla = "P_Com_Fun_SLA",
                    p.fun.comp.N = "P_Com_Fun_N",
                    p.fun.comp.P = "P_Com_Fun_P",
                    p.fun.comp.Ca = "P_Com_Fun_Ca",
                    p.fun.comp.Mg = "P_Com_Fun_Mg",
                    p.fun.comp.lignin = "P_Com_Fun_lignin",
                    p.fun.comp.prim.fiber = "P_Com_Fun_primfiber",
                    p.tax.abund = "P_Abu_Tax_cover",
                    p.fun.abund.bm = "P_Abu_Fun_biomass",
                    p.fun.abund.lai = "P_Abu_Fun_LAI",
                    p.tax.div.0 = "P_Div_Tax_q0",
                    p.tax.div.1 = "P_Div_Tax_q1",
                    p.tax.div.2 = "P_Div_Tax_q2",
                    p.fun.div.0 = "P_Div_Fun_q0",
                    p.fun.div.1 = "P_Div_Fun_q1",
                    p.fun.div.2 = "P_Div_Fun_q2",
                    i.tax.comp.1 = "H_Com_Tax_axis1",
                    i.tax.comp.2 = "H_Com_Tax_axis2",
                    i.tax.comp.3 = "H_Com_Tax_axis3",
                    i.fun.comp.1 =  "H_Com_Fun_axis1",
                    i.fun.comp.2 =  "H_Com_Fun_axis2",
                    i.fun.comp.3 = "H_Com_Fun_axis3",
                    i.tax.abund.sort = "H_Abu_Tax_all",
                    i.tax.abund.det = "H_Abu_Tax_id",
                    i.fun.abund = "H_Abu_Fun_met",
                    i.tax.div.0 = "H_Div_Tax_q0",
                    i.tax.div.1 = "H_Div_Tax_q1",
                    i.tax.div.2 = "H_Div_Tax_q2",
                    i.fun.div.0 = "H_Div_Fun_q0",
                    i.fun.div.1 = "H_Div_Fun_q1",
                    i.fun.div.2 = "H_Div_Fun_q2")

# species level
metrics_spec <- c(system = "System",
                    region = "Region",
                    plot = "Plot",
                    Host_Name_std = "Pl_spec",
                    a.herb2.prop = "Herb_rate_total",
                    a.chewing.prop = "Herb_rate_chewing",
                    a.scraping.prop = "Herb_rate_scraping",
                    a.sucking.prop = "Herb_rate_sucking",
                    a.mines.prop = "Herb_rate_mines",
                    a.gallstot.prop = "Herb_rate_galls",
                    ForMI = "ForMI_std",
                    Inonat = "Inonat_std",
                    Idwcut = "Idwcut_std",
                    Iharv = "Iharv_std",
                    LUI = "LUI_std",
                    F_std = "Fertilization_std",
                    M_std = "Mowing_std",
                    G_std = "Grazing_std",
                    p.fgr.comp.grass = "P_Com_Fun_grass",
                    p.fgr.comp.legume = "P_Com_Fun_legume",
                    p.fgr.comp.forb = "P_Com_Fun_forb",
                    p.fgr.comp.tree = "P_Com_Fun_tree",
                    p.fgr.comp.geophyt = "P_Com_Fun_geophyte",
                    p.fun.comp.ldmc = "P_Com_Fun_LDMC",
                    p.fun.comp.sla = "P_Com_Fun_SLA",
                    p.tax.abund = "P_Abu_Tax_cover",
                    p.fun.abund = "P_Abu_Fun_biomass",
                    i.tax.comp.iszero = "H_Com_0",
                    i.tax.comp.1 = "H_Com_Tax_axis1",
                    i.tax.comp.2 = "H_Com_Tax_axis2",
                    i.tax.comp.3 = "H_Com_Tax_axis3",
                    i.fun.comp.1 =  "H_Com_Fun_axis1",
                    i.fun.comp.2 =  "H_Com_Fun_axis2",
                    i.fun.comp.3 = "H_Com_Fun_axis3",
                    i.tax.abund = "H_Abu_Tax_id",
                    i.fun.abund = "H_Abu_Fun_met",
                    i.tax.div.0 = "H_Div_Tax_q0",
                    i.tax.div.1 = "H_Div_Tax_q1",
                    i.tax.div.2 = "H_Div_Tax_q2",
                    i.fun.div.0 = "H_Div_Fun_q0",
                    i.fun.div.1 = "H_Div_Fun_q1",
                    i.fun.div.2 = "H_Div_Fun_q2")

# load data --------------------------------------------------------------------

# check Community_metrics.R for calculation of community metrics in the following data sets

# Community-level data (BEXIS ID 31412, https://www.bexis.uni-jena.de/)
d_comm <- read.table("Data/Dat_communitylevel.txt", header = T)

d_comm <- d_comm %>% 
  mutate_all(~ifelse(. == "nd", NA, .)) %>% 
  mutate_at(vars(-c(System, Region, Plot)), ~ as.numeric(.)) %>% 
  rename(metrics_comm) %>% 
  group_by(system) %>% 
  mutate_at(vars(-starts_with("a.")),
            ~ f_std(.)) %>% 
  ungroup() %>% 
  mutate(a.herb2.prop.comm_log = log(a.herb2.prop.comm))

# Species-level data (BEXIS ID 31413, https://www.bexis.uni-jena.de/)
d_spec <- read.table("Data/Dat_specieslevel.txt", header = T)

d_spec <- d_spec %>% 
  mutate_all(~ifelse(. == "nd", NA, .)) %>% 
  mutate_at(vars(-c(System, Region, Plot, Pl_spec, H_Com_0)), ~ as.numeric(.)) %>% 
  rename(metrics_spec) %>% 
  mutate(a.herb2.prop_log = log(a.herb2.prop + 0.001))

d_spec <- d_spec %>% 
  mutate_at(vars(p.fgr.comp.grass, p.fgr.comp.forb, p.fgr.comp.tree, 
                 p.fgr.comp.geophyt, p.fgr.comp.legume, i.tax.comp.iszero),
            ~ as.factor(.))

# global parameters ------------------------------------------------------------

# text size of figures
s.axis.text <- 8
s.axis.title <- 10

# colours
col.plgroups <- c(fern = "#FB6542", grass = "#3F681C", forb = "#FFBB00", 
                  tree = "#375E97", geophyt = "#C00000", legume = "#9437FF")

# labels
labs.plgroups <- c(fern = "Ferns", grass = "Grasses", forb = "Forbs", 
                   tree = "Trees", geophyt = "Geophytes", legume = "Legumes")

varlabs <- c("p.tax.comp.1" = expression("Com"["Tax" ]~"axis1"),
             "p.tax.comp.2" = expression("Com"["Tax" ]~"axis2"),
             "p.tax.comp.3" = expression("Com"["Tax" ]~"axis3"),
             "p.fgr.comp.grass" = expression("Com"["Fun" ]~"grass"),
             "p.fgr.comp.grass1" = expression("Com"["Fun" ]~"grass"),
             "p.fgr.comp.forb" = expression("Com"["Fun" ]~"forb"),
             "p.fgr.comp.tree" = expression("Com"["Fun" ]~"tree"),
             "p.fgr.comp.geophyt" = expression("Com"["Fun" ]~"geophyte"),
             "p.fgr.comp.forb1" = expression("Com"["Fun" ]~"forb"),
             "p.fgr.comp.tree1" = expression("Com"["Fun" ]~"tree"),
             "p.fgr.comp.geophyt1" = expression("Com"["Fun" ]~"geophyte"),
             "p.fgr.comp.legume" = expression("Com"["Fun" ]~"legume"),
             "p.fgr.comp.legume1" = expression("Com"["Fun" ]~"legume"),
             "p.fun.comp.ldmc" = expression("Com"["Fun" ]~"LDMC"),
             "p.fun.comp.sla" = expression("Com"["Fun" ]~"SLA"),
             "p.tax.abund" = expression("Abu"["Tax" ]~"cover"),
             "p.fun.abund" = expression("Abu"["Fun" ]~"biomass"),
             "p.fun.abund.bm" = expression("Abu"["Fun" ]~"biomass"),
             "p.fun.abund.lai" = expression("Abu"["Fun" ]~"LAI"),
             "p.tax.div.0" = expression("Div"["Tax" ]~"q=0"),
             "p.tax.div.1" = expression("Div"["Tax" ]~"q=1"),
             "p.tax.div.2" = expression("Div"["Tax" ]~"q=2"),
             "p.fun.div.0" = expression("Div"["Fun" ]~"q=0"),
             "p.fun.div.1" = expression("Div"["Fun" ]~"q=1"),
             "p.fun.div.2" = expression("Div"["Fun" ]~"q=2"),
             "i.tax.comp.iszeroyes" = "Com = 0",
             "i.tax.comp.iszero" = "Com = 0",
             "i.tax.comp.1" =  expression("Com"["Tax" ]~"axis1"),
             "i.tax.comp.2" = expression("Com"["Tax" ]~"axis2"),
             "i.tax.comp.3" =  expression("Com"["Tax" ]~"axis3"),
             "i.fun.comp.1" =  expression("Com"["Fun" ]~"axis1"),
             "i.fun.comp.2" =  expression("Com"["Fun" ]~"axis2"),
             "i.fun.comp.3" = expression("Com"["Fun" ]~"axis3"),
             "i.tax.abund.sort" = expression("Abu"["Tax" ]~"all"),
             "i.tax.abund.det" = expression("Abu"["Tax" ]~"id."),
             "i.tax.abund" = expression("Abu"["Tax" ]~"id."),
             "i.fun.abund" = expression("Abu"["Fun" ]~"met."),
             "i.tax.div.0" = expression("Div"["Tax" ]~"q=0"),
             "i.tax.div.1" = expression("Div"["Tax" ]~"q=1"),
             "i.tax.div.2" = expression("Div"["Tax" ]~"q=2"),
             "i.fun.div.0" = expression("Div"["Fun" ]~"q=0"),
             "i.fun.div.1" = expression("Div"["Fun" ]~"q=1"),
             "i.fun.div.2" = expression("Div"["Fun" ]~"q=2"),
             "Inonat" = "Inonat",
             "Idwcut" = "Idwcut",
             "Iharv" = "Iharv",
             "F_std" = "Fertilisation",
             "M_std" = "Mowing",
             "G_std" = "Grazing",
             "regionSCH" = "region SCH",
             "regionHAI" = "region HAI",
             "a.herb2.prop.comm_log" = "Herbivory rate",
             "a.herb2.prop_log" = "Herbivory rate",
             "direct" = "direct",
             "indirect" = "indirect",
             "net_effect" = "net effect",
             "R2" = expression(R^2),
             "R2m" = expression(marg.~R^2),
             "R2c" = expression(cond.~R^2)) 

# functions --------------------------------------------------------------------

# R2 for glmmTMB models from github --------------------------------------------.

source("https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/glmmTMB/inst/misc/rsqglmm.R")

# extract coefficients and confidence intervals from models --------------------.

f.p.coef.ci <- function(mod){
  
  if (class(mod) == "glmmTMB"){
    p <- summary(mod)$coefficients$cond  %>% 
      as.data.frame() %>% 
      rownames_to_column("var") %>% 
      select(var, `Pr(>|z|)`) %>% 
      rename(p_value = `Pr(>|z|)`)
    
    coef.ci <- confint(mod, parm = "beta_", method = "profile") %>% 
      cbind(Estimate = as.data.frame(confint(mod, parm = "beta_"))$Estimate) %>% 
      as.data.frame() %>% 
      rownames_to_column("var")
    
    p %>% 
      full_join(coef.ci, by = "var")
  } else if (class(mod) == "lm"){
    p <- summary(mod)$coefficients %>% 
      as.data.frame() %>% 
      rownames_to_column("var") %>% 
      select(var, `Pr(>|t|)`, Estimate) %>%
      rename(p_value = `Pr(>|t|)`)
    
    coef.ci <- confint(mod) %>% 
      as.data.frame() %>% 
      rownames_to_column("var")
    
    p %>% 
      full_join(coef.ci, by = "var")
  }
  
}

# Model selection adjusted  ----------------------------------------------------.
# (deals partly with collinearity, also works for glmmTMB)

f.step.adj <- function(mod, fixed, fixcomb = F){
  
  resp <- formula(mod)[2] %>% as.character
  
  vars <- attr(terms(mod), "term.labels")
  vars <- vars[!vars %in% fixed]
  
  mod_up <- mod
  vars_up <- vars
  goon <- T
  while(goon){ # cycle through updating models
    d_AIC <- data.frame(var1 = NA, var2 = NA, AIC = AIC(mod_up)) # initiate data frame
    
    # function 1: remove only 1 predictor variable at a time
    f.var1 <- function(var_i){
      mod_i <- update(mod_up, as.formula(paste(". ~ . -", var_i)))
      
      data.frame(var1 = var_i, var2 = NA, AIC = AIC(mod_i),
                 stringsAsFactors = F)
    }
    out1 <- mclapply(vars_up, f.var1) # apply the function to all variables
    out1 <- do.call(rbind, out1)
    
    d_AIC <- bind_rows(d_AIC, out1)
    
    # first, select variables that are highly collinear with another
    mcor <- performance::check_collinearity(mod_up)
    mcor <- mcor %>% 
      filter(VIF > 10, # threshold from the function guide
             Parameter %in% vars)
    
    if (nrow(mcor) > 1){
      # combinations of 2 variables out of these collinear set (cannot see which are colinear with which, so just try all combinations)
      doublecombs <- combn(mcor$Parameter, 2)
      
      f.var2 <- function(vars_i){
        mod_i <- update(mod_up, as.formula(paste(". ~ . -", doublecombs[1, vars_i],
                                                 "-", doublecombs[2, vars_i])))
        data.frame(var1 = doublecombs[1, vars_i], 
                   var2 = doublecombs[2, vars_i], 
                   AIC = AIC(mod_i),
                   stringsAsFactors = F)
      }
      out2 <- mclapply(1:ncol(doublecombs), f.var2)
      out2 <- do.call(rbind, out2)
      
      d_AIC <- bind_rows(d_AIC, out2)
    }
    
    # keep factor denoting that there are no species sampled in the model
    # is excluded below if no compositional variables are left
    if (!isFALSE(fixcomb)){
      d_AIC <-  d_AIC %>%
        filter((is.na(var1) | var1 != fixcomb) |
                 (!is.na(var2) & var2 != fixcomb))
    }
    
    sel <- d_AIC %>% 
      filter(AIC == min(AIC))
    print(sel)
    if (is.na(sel$var1)){
      goon <- F
    } else if(is.na(sel$var2)){
      mod_up <- update(mod_up, as.formula(paste(". ~ . -", sel$var1)))
      vars_up <- vars_up[!vars_up %in% c(sel$var1)]
    } else{
      mod_up <- update(mod_up, as.formula(paste(". ~ . -", sel$var1, "-", sel$var2)))
      vars_up <- vars_up[!vars_up %in% c(sel$var1, sel$var2)]
    }
    
    
    # if all compositional variables are excluded, also exclude the factor denoting whether there are any herbivores of the plant species recorded
    if (!isFALSE(fixcomb)){
      ncomp <- sum(grepl("i.tax.comp", vars_up)) + sum(grepl("i.fun.comp", vars_up))
      if (ncomp == 0){
        mod_up <- update(mod_up, as.formula(paste(". ~ . -", fixcomb)))
        vars_up <- vars_up[!vars_up %in% c(fixcomb)]
      }
    }
    # go on
  }
  mod_up
}

# stepwise SEM improvement based on missing pathways ---------------------------.

# assumptions
# 1. bottom-up control (plants --> insects)
# 2. comp --> abund --> div (comp first necessary for plants!)
# 3. within category, tax --> fun (also necessary for plants)
# 4. if all categories equal, <-->

f.SEM.step <- function(modlist_i, data, fixcomb = NULL, ...){
  missing_paths <- sem.missing.paths(modlist_i, data, ...)[, 1:6]

  missing_paths <- missing_paths %>% 
    filter(p.value <= 0.05) %>% 
    mutate(missing.path = as.character(missing.path))
  
  modlist_updated <- modlist_i
  
  if (!is.null(fixcomb) & nrow(missing_paths) > 0) {
    missing_paths <-  missing_paths %>% 
      rowwise() %>% 
      filter(!grepl(fixcomb, strsplit(missing.path, "~")[[1]][2])) %>% 
      ungroup()
  }
  
  while (nrow(missing_paths) > 0){
    d_to_add <- data.frame(response = character(0), predictor = character(0))
    
    
    for (i in 1:nrow(missing_paths)){
      vars_i <- as.character(missing_paths$missing.path[i])
      
      var1_i <- substr(vars_i, 1, regexpr(" ~ ", vars_i) - 1)
      var2_i <- substr(vars_i, regexpr(" ~ ", vars_i) + 3, regexpr(" \\+ ", vars_i) - 1)
      
      if (substr(var1_i, 1, 1) == "p" & substr(var2_i, 1, 1) == "i"){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
      } else if (substr(var2_i, 1, 1) == "p" & substr(var1_i, 1, 1) == "i"){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
      } else if (grepl("comp", var1_i) & !grepl("comp", var2_i)){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
        if (!is.null(fixcomb) & substr(var1_i, 1, 1) == "i") {
          d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = fixcomb)
        }
      } else if (grepl("comp", var2_i) & !grepl("comp", var1_i)){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
        if (!is.null(fixcomb) & substr(var2_i, 1, 1) == "i") {
          d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = fixcomb)
        }
      } else if (grepl("abu", var1_i) & grepl("div", var2_i)){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
      } else if (grepl("abu", var2_i) & grepl("div", var1_i)){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
      } else if (grepl("tax", var1_i) & grepl("fgr", var2_i)){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
      } else if (grepl("tax", var1_i) & grepl("fun", var2_i)){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
      } else if (grepl("fgr", var1_i) & grepl("fun", var2_i)){
        d_to_add <- d_to_add %>% add_row(response = var2_i, predictor = var1_i)
      } else if (grepl("tax", var2_i) & grepl("fgr", var1_i)){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
      } else if (grepl("tax", var2_i) & grepl("fun", var1_i)){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
      } else if (grepl("fgr", var2_i) & grepl("fun", var1_i)){
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i)
      } else {
        d_to_add <- d_to_add %>% add_row(response = var1_i, predictor = var2_i) %>% 
          add_row(response = var2_i, predictor = var1_i)
      }
    }
    
    mod_responses <- sapply(modlist_i,
                            function(x) as.character(formula(x)[2]))
    
    for (i in 1:nrow(d_to_add)){
      res_i <- d_to_add$response[i]
      pred_i <- d_to_add$predictor[i]
      
      mod_nr_i <- which(mod_responses == res_i)
      
      if (!grepl(pred_i, formula(modlist_updated[[mod_nr_i]])[3])){
        modlist_updated[[mod_nr_i]] <- update(modlist_updated[[mod_nr_i]], paste(". ~ . +", pred_i))
      }
      
    }
    
    missing_paths <- sem.missing.paths(modlist_updated, data, ...)[, 1:6]

    missing_paths <- missing_paths %>% 
      filter(p.value <= 0.05) %>% 
      mutate(missing.path = as.character(missing.path))
    
    if (!is.null(fixcomb) & nrow(missing_paths) > 0) {
      missing_paths <- missing_paths %>% 
        rowwise() %>% 
        filter(!grepl(fixcomb, strsplit(missing.path, "~")[[1]][2])) %>% 
        ungroup()
    }
  }
  
  modlist_updated  
}

# function that adjust the model coefficients of logistic regression models ----.
# so that they represent changes in probability

f_adjust_semcoefs <- function(d_semcoefs, modlist){
  
  binomials_ind <- which(lapply(modlist, function(x) family(x)$family) == "binomial")
  binomials <- unlist(lapply(modlist, function(x) as.character(formula(x))[2])[binomials_ind])
  
  
  for (i in 1:length(binomials_ind)){
    mod <- modlist[[binomials_ind[i]]]
    
    # get model data
    dat_i <- mod$frame
    
    dat_i <- dat_i %>% 
      mutate_if(is.character, as.factor)
    
    # get response & predictors
    resp_i <- binomials[i]
    pred_i <- names(dat_i)[names(dat_i) != resp_i]
    pred_i <- pred_i[!pred_i %in% c("plot", "Host_Name_std")]
    
    # create empty template to fill with according data
    newdat_i <- dat_i %>% 
      mutate_if(is.numeric, ~0)
    newdat_i[, names(newdat_i) %in% c("plot", "Host_Name_std")] <- NA
    newdat_i <- newdat_i[, names(newdat_i) != resp_i]
    newdat_i <- newdat_i %>% distinct()
    
    out_i <- data.frame()
    for (j in pred_i){
      
      # numeric variables
      if (is.numeric(dat_i[, j])){
        # "low" state data
        newdat_j_low <- newdat_i %>% 
          mutate(!! sym(j) := -.5)
        # "high" state data
        newdat_j_high <- newdat_i %>% 
          mutate(!! sym(j) := .5)
        
        estimate_j <- mean(predict(mod, newdata = newdat_j_high, type = "response")) -
          mean(predict(mod, newdata = newdat_j_low, type = "response"))
        
        out_i <- data.frame(predictor = j,
                            estimate = estimate_j) %>% 
          bind_rows(out_i, .)
        
        # factor variables
      } else if (is.factor(dat_i[, j])){
        for (level_m in levels(dat_i[, j])[-1]){
          # "low" level data
          newdat_j_low <- newdat_i %>% 
            filter(!!sym(j) == levels(dat_i[, j])[1])
          # "high" level data
          newdat_j_high <- newdat_i %>% 
            filter(!!sym(j) == level_m)
          
          # if there are several p.tax predictors, all should be 0 at low level (reference group!)
          if (grepl("p.tax", j)){
            newdat_j_low <- newdat_j_low %>% 
              mutate_at(vars(contains("p.tax")), ~ "0") %>% 
              distinct()
          }
          
          estimate_j <- mean(predict(mod, newdata = newdat_j_high, type = "response")) -
            mean(predict(mod, newdata = newdat_j_low, type = "response"))
          
          out_i <- data.frame(predictor = paste0(j, level_m),
                              estimate = estimate_j) %>% 
            bind_rows(out_i, .)
        }
      } else {
        warning("something wrong") # neither numeric nor factor variable
      }
    }
    
    out_i <- out_i %>% 
      mutate(response = resp_i,
             std.error = "dummy")
    
    # make changes to coefficients data frame
    d_semcoefs <- d_semcoefs %>% 
      left_join(out_i, by = c("predictor", "response"), suffix = c("", ".adapted")) %>% 
      mutate(estimate = ifelse(!is.na(estimate.adapted), estimate.adapted, estimate),
             std.error = ifelse(!is.na(std.error.adapted), NA, std.error)) %>% 
      select(-c(estimate.adapted, std.error.adapted))
  }
  
  d_semcoefs
}

# Labeller for ggplots  --------------------------------------------------------.

ggplot_labeller <- function(label, angle, fontsize, fill, frame = F, hjust = "middle", parse = F){
  p <- ggplot(data.frame(x = 0, y = 0))
  
  if(frame != F) p <- p + geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = fill, color = frame)
  
  p <- p +
    annotate(geom = "text", label = label, x = 1, y = 0, size = fontsize/ggplot2:::.pt, angle = angle, hjust = hjust, parse = parse) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = fill),
          plot.margin = unit(c(0,0,0,0),"cm"))
  
  return(p)
}

# Variograms for regional subsets ----------------------------------------------.

f_vario <- function(mod, data, subset_region){
  d_res <-data[names(resid(mod)), ] %>% 
    mutate(residuals = resid(mod)) %>% 
    filter(region == subset_region) %>% 
    select(residuals, RW, HW)
  
  coordinates(d_res) <- c("RW", "HW")
  
  variogram(residuals ~ 1, d_res) %>% 
    mutate(region = subset_region)
}

# show random effects for each plot on coordinate grid -------------------------.

f_ranefplots <- function(mod){
  ranef_plots <- list()
  for (reg_i in c("Alb", "Hai", "Sch")){
    ranef_plots[[reg_i]] <- ranef(mod)$cond[["plot"]] %>% 
      rename(ranef = `(Intercept)`) %>% 
      rownames_to_column("PlotID") %>% 
      mutate(sign = factor(sign(ranef), label = c("neg", "pos"))) %>% 
      left_join(d_comm %>% select(plot, region, RW, HW), by = c(PlotID = "plot")) %>%
      mutate(RW = RW / 1000,
             HW = HW / 1000) %>% 
      filter(region == reg_i) %>% 
      ggplot(aes(x = RW, y = HW, size = abs(ranef), col = sign)) +
      geom_point(shape = 1) +
      coord_fixed() +
      facet_wrap(~ region) +
      theme(legend.position = "none")
  }
  plot_grid(plotlist = ranef_plots)
}

################################################################################.
# ANALYSIS ---------- ##########################################################
################################################################################.

# Land-use models ##############################################################
################################################################################.

# i.e., models used to analyze land-use effects on herbivory rates
# not accounting for plant/herbivore communities

# community-level models -------------------------------------------------------

# forests ----------------------------------------------------------------------.

mod.comm.formi <- lm(log(a.herb2.prop.comm) ~ ForMI + region,
                     data = d_comm)

summary(mod.comm.formi)

mod.comm.formi.sep <- lm(log(a.herb2.prop.comm) ~ Inonat + Iharv + Idwcut + region,
                         data = d_comm)

summary(mod.comm.formi.sep)

# grasslands -------------------------------------------------------------------.


mod.comm.lui <- lm(log(a.herb2.prop.comm) ~ LUI + region,
                   data = d_comm)
summary(mod.comm.lui)

mod.comm.lui.sep <- lm(log(a.herb2.prop.comm) ~ G_std + M_std + F_std + region,
                       data = d_comm)

summary(mod.comm.lui.sep)

# extract coefficients and confidence intervals  -------------------------------.

coef.ci.plot <- f.p.coef.ci(mod.comm.formi) %>% 
  mutate(model = "comb",
         system = "Forest") %>% 
  bind_rows(f.p.coef.ci(mod.comm.formi.sep) %>% 
              mutate(model = "sep",
                     system = "Forest")) %>% 
  bind_rows(f.p.coef.ci(mod.comm.lui) %>% 
              mutate(model = "comb",
                     system = "Grassland")) %>% 
  bind_rows(f.p.coef.ci(mod.comm.lui.sep) %>% 
              mutate(model = "sep",
                     system = "Grassland")) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(var != "(Intercept)",
         !grepl("region", var)) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                          "LUI", "G_std", "M_std", "F_std")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                     "LUI", "Grazing", "Mowing", "Fertilisation"))),
         sign = ifelse(p_value < 0.05, "sign", "ns"))

# species-level models ---------------------------------------------------------

# forests ----------------------------------------------------------------------.

mod.sample.formi <- glmmTMB(a.herb2.prop_log ~ ForMI + region +
                              (1 | plot) + (1 | Host_Name_std),
                            data = d_spec)

summary(mod.sample.formi)


mod.sample.formi.sep <- glmmTMB(a.herb2.prop_log ~ Iharv + Inonat + Idwcut + region +
                                  (1 | plot) + (1 | Host_Name_std),
                                data = d_spec)

summary(mod.sample.formi.sep)


# grasslands -------------------------------------------------------------------.

mod.sample.lui <- glmmTMB(a.herb2.prop_log ~ LUI + region +
                            (1 | plot) + (1 | Host_Name_std),
                          data = d_spec)

summary(mod.sample.lui)

mod.sample.lui.sep <- glmmTMB(a.herb2.prop_log ~ G_std + M_std + F_std + region +
                                (1 | plot) + (1 | Host_Name_std),
                              data = d_spec)

summary(mod.sample.lui.sep)

# extract coefficients and confidence intervals  -------------------------------.

coef.ci.sample <- f.p.coef.ci(mod.sample.formi) %>% 
  mutate(model = "comb",
         system = "Forest") %>% 
  bind_rows(f.p.coef.ci(mod.sample.formi.sep) %>% 
              mutate(model = "sep",
                     system = "Forest")) %>% 
  bind_rows(f.p.coef.ci(mod.sample.lui) %>% 
              mutate(model = "comb",
                     system = "Grassland")) %>% 
  bind_rows(f.p.coef.ci(mod.sample.lui.sep) %>% 
              mutate(model = "sep",
                     system = "Grassland")) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(var != "(Intercept)",
         !grepl("region", var)) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                          "LUI", "G_std", "M_std", "F_std")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                     "LUI", "Grazing", "Mowing", "Fertilisation"))),
         sign = ifelse(p_value < 0.05, "sign", "ns"))

# plotting (Fig. 3) ------------------------------------------------------------

coef.ci.plot %>% 
  mutate(level = "Community level") %>% 
  bind_rows(coef.ci.sample %>% 
              mutate(level = "Species level")) %>% 
  ggplot(aes(x = Estimate, y = var)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 3.5, lty = 3) +
  geom_segment(aes(x = lower, xend = upper, yend = var, size = sign, col = sign)) +
  geom_point(col = "white", size = 3) +
  geom_point(aes(col = sign), size = 2) +
  facet_wrap(level ~ system, scales = "free_y") +
  xlab("Standardised model coefficients (± CI)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = s.axis.title),
        axis.text = element_text(size = s.axis.text),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing.y = unit(60, "pt"),
        plot.margin = margin(40, 2, 2, 2, "pt")) +
  scale_size_manual(values = c(ns = .5, sign = 1.5), guide = F) +
  scale_color_manual(values = c(ns = "grey50", sign = 1), guide = F)

# exclusion of some forests plots ----------------------------------------------

# On some forest plots, insect sampling was not conducted. These plots were not
# included in the analysis of pathways. How do overall models change, if they
# are excluded from these analysis?

# community level --------------------------------------------------------------.

d_comm_sub <- d_comm %>% 
  filter(system == "forest",
         # exclude plots on which insects were not sampled
         !plot %in% c("AEW41", "SEW22", "SEW23", "SEW24", "SEW25", 
                      "SEW26", "SEW27", "SEW28")) %>% 
  mutate_at(vars(ForMI, Iharv, Inonat, Idwcut),
            ~as.numeric((. - mean(., na.rm = T)) / (2 * sd(., na.rm = T))))


mod.comm.formi.sub <- lm(log(a.herb2.prop.comm) ~ ForMI + region,
                         data = d_comm_sub)

summary(mod.comm.formi.sub)


mod.comm.formi.sep.sub <- lm(log(a.herb2.prop.comm) ~ Inonat + Iharv + Idwcut + region,
                             data = d_comm_sub)

summary(mod.comm.formi.sep.sub)

# species level ----------------------------------------------------------------.

d_spec_sub <- 
  d_spec %>% 
  filter(system == "forest") %>% 
  # exclude plots on which insects were not sampled
  filter(!plot %in% c("AEW41", "SEW22", "SEW23", "SEW24", "SEW25", 
                      "SEW26", "SEW27", "SEW28")) %>% 
  select(a.herb2.prop_log, a.specialist.prop, a.generalist.prop, 
         a.chewing.prop, a.scraping.prop, a.sucking.prop, a.mines.prop, a.gallstot.prop,
         plot, Host_Name_std, region, RW, HW) %>% 
  left_join(formi %>% mutate_at(vars(ForMI, Iharv, Inonat, Idwcut),
                                ~as.numeric((. - mean(., na.rm = T)) / 
                                              (2 * sd(., na.rm = T)))),
            by = c(plot = "PlotID"))


mod.sample.formi.sub <- glmmTMB(a.herb2.prop_log ~ ForMI + region +
                                  (1 | plot) + (1 | Host_Name_std),
                                data = d_spec_sub)

summary(mod.sample.formi.sub)


mod.sample.formi.sep.sub <- glmmTMB(a.herb2.prop_log ~ Iharv + Inonat + Idwcut + region +
                                      (1 | plot) + (1 | Host_Name_std),
                                    data = d_spec_sub)

summary(mod.sample.formi.sep.sub)

# extract coefficient estimates ------------------------------------------------.

coef.ci.plot.sub <- f.p.coef.ci(mod.comm.formi.sub) %>% 
  mutate(model = "comb") %>% 
  bind_rows(f.p.coef.ci(mod.comm.formi.sep.sub) %>% 
              mutate(model = "sep")) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(var != "(Intercept)",
         !grepl("region", var)) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut"))),
         sign = ifelse(p_value < 0.05, "sign", "ns"))


coef.ci.sample.sub <- f.p.coef.ci(mod.sample.formi.sub) %>% 
  mutate(model = "comb") %>% 
  bind_rows(f.p.coef.ci(mod.sample.formi.sep.sub) %>% 
              mutate(model = "sep")) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(var != "(Intercept)",
         !grepl("region", var)) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut"))),
         sign = ifelse(p_value < 0.05, "sign", "ns"))

# plotting ---------------------------------------------------------------------.

coef.ci.plot.sub %>% 
  mutate(level = "Community level") %>% 
  bind_rows(coef.ci.sample.sub %>% 
              mutate(level = "Sample level")) %>% 
  ggplot(aes(x = Estimate, y = var)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 3.5, lty = 3) +
  geom_segment(aes(x = lower, xend = upper, yend = var, size = sign, col = sign)) +
  geom_point(col = "white", size = 3) +
  geom_point(aes(col = sign), size = 2) +
  facet_wrap(level ~ ., scales = "free_y") +
  xlab("Standardised model coefficients (± CI)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = s.axis.title),
        axis.text = element_text(size = s.axis.text),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing.y = unit(60, "pt"),
        plot.margin = margin(40, 2, 2, 2, "pt")) +
  scale_size_manual(values = c(ns = .5, sign = 1.5), guide = F) +
  scale_color_manual(values = c(ns = "grey50", sign = 1), guide = F)

# separate models per damage type ----------------------------------------------

# community level --------------------------------------------------------------.

damage_types <- c("a.chewing.prop.comm", "a.gallstot.prop.comm", 
                  "a.mines.prop.comm", "a.scraping.prop.comm", 
                  "a.sucking.prop.comm")

mod_damage_types <- list()
for (res_i in damage_types){
  target <- d_comm %>% 
    mutate(response = !! sym(res_i),
           response = log(response + 0.001))
  
  mod_damage_types[[res_i]]$Forest$comb <- lm(response ~ ForMI + region,
                                              data = target)
  
  mod_damage_types[[res_i]]$Forest$split <- lm(response ~ Inonat + Idwcut + Iharv + region,
                                               data = target)
  
  target <- d_comm %>% 
    mutate(response = !! sym(res_i),
           response = log(response + 0.001))
  
  mod_damage_types[[res_i]]$Grassland$comb <- lm(response ~ LUI + region,
                                                 data = target)
  
  mod_damage_types[[res_i]]$Grassland$split <- lm(response ~ G_std + M_std + F_std + region,
                                                  data = target)
  
}

labels_dt <- c(Forest.Chewing = " Chewing", Grassland.Chewing = "Chewing ",
               Forest.Scraping = " Scraping", Grassland.Scraping = "Scraping ",
               Forest.Sucking = " Sucking", Grassland.Sucking = "Sucking ",
               Forest.Mines = " Mines", Grassland.Mines = "Mines ",
               Forest.Galls = " Galls")

mod_damage_types %>% 
  map(function(x) map(x, 
                      function(y) map(y, f.p.coef.ci) %>% 
                        enframe("model") %>% 
                        unnest(cols = everything())) %>% 
        enframe("system") %>% 
        unnest(cols = everything())) %>% 
  enframe("damage_type") %>% 
  unnest(cols = everything()) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(!var %in% c("(Intercept)", "regionHAI", "regionSCH"),
         !(damage_type == "a.gallstot.prop.comm" & system == "Grassland")) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                          "LUI", "G_std", "M_std", "F_std")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                     "LUI", "Grazing", "Mowing", "Fertilisation"))),
         damage_type = factor(damage_type, 
                              levels = c("a.chewing.prop.comm", "a.scraping.prop.comm", "a.sucking.prop.comm",
                                         "a.mines.prop.comm", "a.gallstot.prop.comm"),
                              labels = c("Chewing", "Scraping", "Sucking", "Mines", "Galls")),
         sign = ifelse(p_value < 0.05, "sign", "ns"),
         ID = factor(interaction(system, damage_type), labels = labels_dt)) %>% 
  ggplot(aes(x = Estimate, y = var)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 3.5, lty = 3) +
  geom_segment(aes(x = lower, xend = upper, yend = var, size = sign, col = sign)) +
  geom_point(col = "white", size = 3) +
  geom_point(aes(col = sign), size = 2) +
  facet_wrap(~ ID, scales = "free_y", ncol = 2) +
  xlab("Standardised model coefficients (± CI)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = s.axis.title),
        axis.text = element_text(size = s.axis.text),
        strip.text = element_text(size = s.axis.title,
                                  margin = margin(.5,0,1,0, "mm")),
        plot.margin = margin(40, 2, 2, 2, "pt")) +
  scale_size_manual(values = c(ns = .5, sign = 1.5), guide = F) +
  scale_color_manual(values = c(ns = "grey50", sign = 1), guide = F)

# species level ----------------------------------------------------------------.

damage_types <- c("a.chewing.prop", "a.gallstot.prop", "a.mines.prop", 
                  "a.scraping.prop", "a.sucking.prop")

mod_damage_types_sample <- list()
for (res_i in damage_types){
  target <- d_spec %>% 
    mutate(response = !! sym(res_i),
           response = log(response + 0.001))
  
  mod_damage_types_sample[[res_i]]$Forest$comb <- glmmTMB(response ~ ForMI + region +
                                                            (1 | plot) + (1 | Host_Name_std),
                                                          data = target)
  
  mod_damage_types_sample[[res_i]]$Forest$split <- glmmTMB(response ~ Inonat + Idwcut + Iharv + region +
                                                             (1 | plot) + (1 | Host_Name_std),
                                                           data = target)
  
  if (res_i != "a.gallstot.prop"){
    target <- d_spec %>% 
      mutate(response = !! sym(res_i),
             response = log(response + 0.001))
    
    mod_damage_types_sample[[res_i]]$Grassland$comb <- glmmTMB(response ~ LUI + region +
                                                                 (1 | plot) + (1 | Host_Name_std),
                                                               data = target)
    
    mod_damage_types_sample[[res_i]]$Grassland$split <- glmmTMB(response ~ G_std + M_std + F_std + region +
                                                                  (1 | plot) + (1 | Host_Name_std),
                                                                data = target)
  }
  
}

mod_damage_types_sample %>% 
  map(function(x) map(x, 
                      function(y) map(y, f.p.coef.ci) %>% 
                        enframe("model") %>% 
                        unnest(cols = everything())) %>% 
        enframe("system") %>% 
        unnest(cols = everything())) %>% 
  enframe("damage_type") %>% 
  unnest(cols = everything()) %>% 
  rename(lower = `2.5 %`,
         upper = `97.5 %`) %>% 
  filter(!var %in% c("(Intercept)", "regionHAI", "regionSCH")) %>% 
  mutate(var = factor(var, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                          "LUI", "G_std", "M_std", "F_std")),
                      labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                     "LUI", "Grazing", "Mowing", "Fertilisation"))),
         damage_type = factor(damage_type,
                              levels = c("a.chewing.prop", "a.scraping.prop", "a.sucking.prop",
                                         "a.mines.prop", "a.gallstot.prop"),
                              labels = c("Chewing", "Scraping", "Sucking", "Mines", "Galls")),
         sign = ifelse(p_value < 0.05, "sign", "ns"),
         ID = factor(interaction(system, damage_type), labels = labels_dt)) %>% 
  ggplot(aes(x = Estimate, y = var)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 3.5, lty = 3) +
  geom_segment(aes(x = lower, xend = upper, yend = var, size = sign, col = sign)) +
  geom_point(col = "white", size = 3) +
  geom_point(aes(col = sign), size = 2) +
  facet_wrap(~ ID, scales = "free_y", ncol = 2) +
  xlab("Standardised model coefficients (± CI)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = s.axis.title),
        axis.text = element_text(size = s.axis.text),
        strip.text = element_text(size = s.axis.title,
                                  margin = margin(.5,0,1,0, "mm")),
        plot.margin = margin(40, 2, 2, 2, "pt")) +
  scale_size_manual(values = c(ns = .5, sign = 1.5), guide = F) +
  scale_color_manual(values = c(ns = "grey50", sign = 1), guide = F)


# models for a selection of species  -------------------------------------------

# forests ----------------------------------------------------------------------.

n.imp <- 10 # minimum number of observations per plant
range.prop <- .66 # proportion of ForMI range that should be covered
range.formi <- range.prop * diff(range(d_spec$ForMI, na.rm = T))

important.pl.species.abbr <-
  d_spec %>%
  filter(!is.na(ForMI)) %>%
  group_by(pl.species.abbr) %>%
  filter(n() > n.imp,
         diff(range(ForMI)) >= range.formi) %>%
  ungroup() %>% 
  select(pl.species.abbr) %>%
  distinct()

f.stats.herb.pool <- list()
for (pl.sp in important.pl.species.abbr$pl.species.abbr){
  target <-
    d_spec %>%
    filter(pl.species.abbr == pl.sp, !is.na(ForMI)) %>% 
    droplevels() %>%
    mutate(region = as.factor(region))
  
  # comb 
  if (nlevels(target$region) == 1){ # region not included when there is only one
    f.stats.herb.pool[[pl.sp]][["comb"]] <-
      lm(log(a.herb2.prop + 0.001) ~ ForMI,
         data = target)
  } else {
    
    f.stats.herb.pool[[pl.sp]][["comb"]] <-
      lm(log(a.herb2.prop + 0.001) ~ ForMI + region,
         data = target)
    
  }
  
  # split 
  if (nlevels(target$region) == 1){ # region not included when there is only one
    f.stats.herb.pool[[pl.sp]][["split"]] <-
      lm(log(a.herb2.prop + 0.001) ~ Iharv + Inonat + Idwcut,
         data = target)
  } else {
    
    f.stats.herb.pool[[pl.sp]][["split"]] <-
      lm(log(a.herb2.prop + 0.001) ~ Iharv + Inonat + Idwcut + region,
         data = target)
    
  }
  
}

# function to extract results
f.extract <- function(x) {
  out <- c()
  
  CIs <- confint(x$comb, level = .95)
  
  if (class(x$comb) == "lm") CIs <- cbind(CIs, Estimate = coefficients(x$comb))
  
  out$lower <-  CIs[grepl("ForMI", rownames(CIs)), "2.5 %"]
  out$upper <-  CIs[grepl("ForMI", rownames(CIs)), "97.5 %"]
  out$est <-  CIs[grepl("ForMI", rownames(CIs)), "Estimate"]
  
  ANOVA <- car::Anova(x$comb, type = "III")
  out$p.value <- ANOVA["ForMI", grepl("Pr", names(ANOVA))]
  out$n <- nobs(x$comb)
  out$explanatory <- "ForMI"
  out$combsplit <- "comb"
  
  out.df <- data.frame(out, stringsAsFactors = F)
  
  CIs <- confint(x$split, level = .95, parm = 1:4)
  
  if (class(x$split) == "lm") CIs <- cbind(CIs, Estimate = coefficients(x$split)[1:4])
  
  for (var in c("Iharv", "Inonat", "Idwcut")){
    out <- c()
    
    out$lower <-  CIs[grepl(var, rownames(CIs)), "2.5 %"]
    out$upper <-  CIs[grepl(var, rownames(CIs)), "97.5 %"]
    out$est <-  CIs[grepl(var, rownames(CIs)), "Estimate"]
    
    ANOVA <- car::Anova(x$split, type = "III")
    out$p.value <- ANOVA[var, grepl("Pr", names(ANOVA))]
    out$n <- nobs(x$split)
    out$explanatory <- var
    out$combsplit <- "split"
    
    out.df <- out.df %>%
      bind_rows(data.frame(out, stringsAsFactors = F))
    
  }
  
  out.df
  
}

# apply function
f.cf.df <- f.stats.herb.pool %>%
  map(f.extract) %>%
  map(cbind.data.frame)  %>%
  enframe("pl.species.abbr") %>%
  unnest(cols = everything())

# adjust data
f.cf.df <- f.cf.df %>%
  mutate(significance = factor(ifelse(p.value <= 0.05, "significant", "ns")),
         combsplit = as.factor(combsplit)) %>% 
  left_join(d_spec %>% 
              select(pl.species.abbr, pl.group) %>% distinct(), 
            by = "pl.species.abbr")

# grasslands -------------------------------------------------------------------.

n.imp <- 10 # minimum number of observations per plant
range.prop <- .66 # proportion of ForMI range that should be covered
range.lui <- range.prop * diff(range(d_spec$LUI, na.rm = T))

important.pl.species.abbr <-
  d_spec %>%
  filter(!is.na(LUI)) %>%
  group_by(pl.species.abbr) %>%
  filter(n() > n.imp,
         diff(range(LUI)) >= range.lui) %>%
  ungroup() %>% 
  select(pl.species.abbr) %>%
  distinct()

g.stats.herb.pool <- list()
for (pl.sp in important.pl.species.abbr$pl.species.abbr){
  target <-
    d_spec %>%
    filter(pl.species.abbr == pl.sp, !is.na(LUI)) %>% 
    droplevels()
  
  # comb 
  if (nlevels(target$region) == 1){ # no mixed effects model when there is only one region
    g.stats.herb.pool[[pl.sp]][["comb"]] <-
      lm(log(a.herb2.prop + 0.001) ~ LUI,
         data = target)
  } else {
    g.stats.herb.pool[[pl.sp]][["comb"]] <-
      lm(log(a.herb2.prop + 0.001) ~ LUI + region,
         data = target)
    
  }
  
  # split 
  if (nlevels(target$region) == 1){ # no mixed effects model when there is only one region
    g.stats.herb.pool[[pl.sp]][["split"]] <-
      lm(log(a.herb2.prop + 0.001) ~ M_std + F_std + G_std,
         data = target)
  } else {
    
    g.stats.herb.pool[[pl.sp]][["split"]] <-
      lm(log(a.herb2.prop + 0.001) ~ M_std + F_std + G_std + region,
         data = target)
    
  }
  
}

# function to extract results
f.extract <- function(x) {
  out <- c()
  
  CIs <- confint(x$comb, level = .95)
  
  if (class(x$comb) == "lm") CIs <- cbind(CIs, Estimate = coefficients(x$comb))
  
  
  out$lower <-  CIs[grepl("LUI", rownames(CIs)), "2.5 %"]
  out$upper <-  CIs[grepl("LUI", rownames(CIs)), "97.5 %"]
  out$est <-  CIs[grepl("LUI", rownames(CIs)), "Estimate"]
  
  ANOVA <- car::Anova(x$comb, type = "III")
  out$p.value <- ANOVA["LUI", grepl("Pr", names(ANOVA))]
  out$n <- nobs(x$comb)
  out$explanatory <- "LUI"
  out$combsplit <- "comb"
  
  out.df <- data.frame(out, stringsAsFactors = F)
  
  CIs <- confint(x$split, level = .95, parm = 1:4)
  
  if (class(x$split) == "lm") CIs <- cbind(CIs, Estimate = coefficients(x$split)[1:4])
  
  
  for (var in c("M_std", "F_std", "G_std")){
    out <- c()
    
    out$lower <-  CIs[grepl(var, rownames(CIs)), "2.5 %"]
    out$upper <-  CIs[grepl(var, rownames(CIs)), "97.5 %"]
    out$est <-  CIs[grepl(var, rownames(CIs)), "Estimate"]
    
    ANOVA <- car::Anova(x$split, type = "III")
    out$p.value <- ANOVA[var, grepl("Pr", names(ANOVA))]
    out$n <- nobs(x$split)
    out$explanatory <- var
    out$combsplit <- "split"
    
    out.df <- out.df %>% 
      bind_rows(data.frame(out, stringsAsFactors = F))
    
  }
  
  out.df
  
}

# apply function
g.cf.df <- g.stats.herb.pool %>% 
  map(f.extract) %>% 
  map(cbind.data.frame)  %>%
  enframe("pl.species.abbr") %>% 
  unnest(cols = everything())

# adjust data
g.cf.df <- g.cf.df %>% 
  mutate(significance = factor(ifelse(p.value <= 0.05, "significant", "ns")),
         combsplit = as.factor(combsplit)) %>% 
  left_join(d_spec %>% 
              select(pl.species.abbr, pl.group) %>% distinct(), 
            by = "pl.species.abbr")

# plotting ---------------------------------------------------------------------

f.cf.df %>%
  filter(combsplit == "comb") %>%
  arrange(-est) %>%
  mutate(pl.species.abbr = factor(pl.species.abbr, levels = unique(pl.species.abbr))) %>%
  {ggplot(., aes(y = pl.species.abbr, x = est, size = significance, col = pl.group)) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 2) +
      geom_hline(yintercept = max(which(.$est > 0)) + 0.5, lty = 3) +
      geom_segment(aes(x = lower, xend = upper, yend = pl.species.abbr)) +
      geom_point(col = "white", size = 3) +
      geom_point(size = 2) +
      scale_size_manual(values = c(ns = 0.5, significant =  1.5), guide = F) +
      scale_color_manual(values = col.plgroups,
                         name = "Plant group",
                         labels = labs.plgroups) +
      geom_label(aes(x = min(lower) - (max(upper) - min(lower))/30, label = n),
                 size = s.axis.text/ggplot2:::.pt, col = 1) +
      xlab(expression("Standardised" ~ italic("ForMI")~ "model coefficient (\u00B1 CI)")) +
      ylab("Plant species") +
      theme(axis.text = element_text(size = s.axis.text),
            axis.title = element_text(size = s.axis.title),
            legend.title = element_text(size = s.axis.title),
            legend.text = element_text(size = s.axis.text),
            axis.text.y = element_text(face = "italic"))}


xlim <- 2 # threshold for plotting region

f.cf.df %>%
  arrange(combsplit, -est) %>%
  mutate(pl.species.abbr = factor(pl.species.abbr, levels = unique(pl.species.abbr))) %>%
  filter(combsplit == "split") %>%
  mutate(explanatory = gsub(".z", "", explanatory),
         explanatory = factor(explanatory, levels = c("Inonat", "Iharv", "Idwcut")),
         outlow = ifelse(lower < -xlim, TRUE, FALSE),
         outhigh = ifelse(upper > xlim, TRUE, FALSE),
         outlow.e = ifelse(est < -xlim, TRUE, FALSE),
         outhigh.e = ifelse(est > xlim, TRUE, FALSE),
         lower = ifelse(lower < -xlim, -xlim, lower),
         upper = ifelse(upper > xlim, xlim, upper)) %>% 
  {ggplot(., aes(y = pl.species.abbr, x = est, size = significance, col = pl.group)) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 2) +
      geom_segment(aes(x = lower, xend = upper, yend = pl.species.abbr)) +
      geom_segment(aes(x = upper, xend = -xlim, yend = pl.species.abbr, alpha = outlow),
                   arrow = arrow(length = unit(5, "pt"))) +
      geom_segment(aes(x = lower, xend = xlim, yend = pl.species.abbr, alpha = outhigh),
                   arrow = arrow(length = unit(5, "pt"))) +
      geom_point(col = "white", size = 3) +
      geom_point(size = 2) +
      geom_text(data = filter(., outhigh.e), 
                aes(x = xlim , y = as.numeric(pl.species.abbr) + .65, 
                    label = round(est, 1)),
                size = s.axis.text/ggplot2:::.pt) +
      scale_size_manual(values = c(ns = 0.5, significant =  1.5), guide = F) +
      scale_color_manual(values = col.plgroups,
                         name = "Plant group",
                         labels = labs.plgroups) +
      geom_label(data = filter(., explanatory == "Inonat"), 
                 aes(x = -xlim - .4 , label = n),
                 size = s.axis.text/ggplot2:::.pt, col = 1) +
      xlab("Standardised model coefficient (\u00B1 CI)") +
      ylab("Plant species") +
      theme(axis.text = element_text(size = s.axis.text),
            axis.title = element_text(size = s.axis.title),
            legend.title = element_text(size = s.axis.title),
            legend.text = element_text(size = s.axis.text),
            axis.text.y = element_text(face = "italic")) +
      scale_x_continuous(limits = c(-(xlim + .5), xlim + .5)) +
      scale_alpha_manual(values = c(0, 1), guide = F) +
      facet_grid(. ~ explanatory)}


g.cf.df %>%
  filter(combsplit == "comb") %>%
  arrange(-est) %>%
  mutate(pl.species.abbr = factor(pl.species.abbr, levels = unique(pl.species.abbr))) %>%
  {ggplot(., aes(y = pl.species.abbr, x = est, size = significance, col = pl.group)) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 2) +
      geom_hline(yintercept = ifelse(any(.$est > 0), max(which(.$est > 0)) + 0.5, 2000), lty = 3) +
      geom_segment(aes(x = lower, xend = upper, yend = pl.species.abbr)) +
      geom_point(col = "white", size = 3) +
      geom_point(size = 2) +
      scale_size_manual(values = c(ns = 0.5, significant =  1.5), guide = F) +
      scale_color_manual(values = col.plgroups,
                         name = "Plant group",
                         labels = labs.plgroups) +
      geom_label(aes(x = min(lower) - (max(upper) - min(lower))/30, label = n),
                 size = s.axis.text/ggplot2:::.pt, col = 1) +
      xlab(expression("Standardised" ~ italic("LUI")~ "model coefficient (\u00B1 CI)")) +
      ylab("Plant species") +
      theme(axis.text = element_text(size = s.axis.text),
            axis.title = element_text(size = s.axis.title),
            legend.title = element_text(size = s.axis.title),
            legend.text = element_text(size = s.axis.text),
            axis.text.y = element_text(face = "italic"))}

xlim <- 3.2 # threshold for plotting region

g.cf.df %>%
  arrange(combsplit, -est) %>%
  mutate(pl.species.abbr = factor(pl.species.abbr, levels = unique(pl.species.abbr))) %>%
  filter(combsplit == "split") %>%
  mutate(explanatory = gsub(".z", "", explanatory),
         explanatory = factor(explanatory, levels = c("G_std", "M_std", "F_std"),
                              labels = c("Grazing", "Mowing", "Fertilisation")),
         outlow = ifelse(lower < -xlim, TRUE, FALSE),
         outhigh = ifelse(upper > xlim, TRUE, FALSE),
         outlow.e = ifelse(est < -xlim, TRUE, FALSE),
         outhigh.e = ifelse(est > xlim, TRUE, FALSE),
         lower = ifelse(lower < -xlim, -xlim, lower),
         upper = ifelse(upper > xlim, xlim, upper)) %>% 
  {ggplot(., aes(y = pl.species.abbr, x = est, size = significance, col = pl.group)) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 2) +
      geom_segment(aes(x = lower, xend = upper, yend = pl.species.abbr)) +
      geom_segment(aes(x = upper, xend = -xlim, yend = pl.species.abbr, alpha = outlow),
                   arrow = arrow(length = unit(5, "pt"))) +
      geom_segment(aes(x = lower, xend = xlim, yend = pl.species.abbr, alpha = outhigh),
                   arrow = arrow(length = unit(5, "pt"))) +
      geom_point(col = "white", size = 3) +
      geom_point(size = 2) +
      scale_size_manual(values = c(ns = 0.5, significant =  1.5), guide = F) +
      scale_color_manual(values = col.plgroups,
                         name = "Plant group",
                         labels = labs.plgroups) +
      geom_label(data = filter(., explanatory == "Grazing"), 
                 aes(x = -xlim - .4 , label = n),
                 size = s.axis.text/ggplot2:::.pt, col = 1) +
      xlab("Standardised model coefficient (\u00B1 CI)") +
      ylab("Plant species") +
      theme(axis.text = element_text(size = s.axis.text),
            axis.title = element_text(size = s.axis.title),
            legend.title = element_text(size = s.axis.title),
            legend.text = element_text(size = s.axis.text),
            axis.text.y = element_text(face = "italic")) +
      scale_x_continuous(limits = c(-(xlim + .5), xlim + .5)) +
      scale_alpha_manual(values = c(0, 1), guide = F) +
      facet_grid(. ~ explanatory)}


# Structural equation modelling ################################################
################################################################################.

# forests community level ------------------------------------------------------

# prepare scaled data set
d_f_path_plotlevel_z <- d_comm %>% 
  select(plot, region, RW, HW, 
         a.herb2.prop.comm_log, ForMI, Inonat, Iharv, Idwcut, 
         p.tax.comp.1, p.tax.comp.2, p.tax.comp.3, 
         p.fgr.comp.grass, p.fgr.comp.forb, p.fgr.comp.tree, p.fgr.comp.geophyt, 
         p.fun.comp.sla, p.fun.comp.ldmc,
         p.tax.abund,
         p.fun.abund.lai, p.fun.abund.bm,
         p.tax.div.0, p.tax.div.1, p.tax.div.2,
         p.fun.div.0, p.fun.div.1, p.fun.div.2,
         i.tax.comp.1, i.tax.comp.2, i.tax.comp.3,
         i.fun.comp.1, i.fun.comp.2, i.fun.comp.3,
         i.tax.abund.sort, i.tax.abund.det,
         i.fun.abund,
         i.tax.div.0, i.tax.div.1, i.tax.div.2,
         i.fun.div.0, i.fun.div.1, i.fun.div.2) %>% 
  filter_all(~!is.na(.)) %>%  # remove all entries with NAs in any of the variables
  mutate_at(vars(-c(plot, region, RW, HW, a.herb2.prop.comm_log)), 
            ~ f_std(.)) 

# ForMI combined ---------------------------------------------------------------.

# define global model
mod_f_path_plot <- lm(a.herb2.prop.comm_log ~ ForMI + region +
                        p.tax.comp.1 + p.tax.comp.2 + p.tax.comp.3 + 
                        p.fgr.comp.grass + p.fgr.comp.forb + p.fgr.comp.tree + 
                        p.fgr.comp.geophyt + 
                        p.fun.comp.sla + p.fun.comp.ldmc +
                        p.tax.abund +
                        p.fun.abund.lai + p.fun.abund.bm +
                        p.tax.div.0 + p.tax.div.1 + p.tax.div.2 +
                        p.fun.div.0 + p.fun.div.1 + p.fun.div.2 +
                        i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                        i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                        i.tax.abund.sort + i.tax.abund.det +
                        i.fun.abund +
                        i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                        i.fun.div.0 + i.fun.div.1 + i.fun.div.2,
                      data = d_f_path_plotlevel_z)

# perform model selection
mod_f_path_plot_sel <- f.step.adj(mod_f_path_plot, c("ForMI", "region"))

# manually check and adjust for collinear variables 
performance::check_collinearity(mod_f_path_plot_sel)

update(mod_f_path_plot_sel, . ~ . - p.fun.div.0) %>% AIC
update(mod_f_path_plot_sel, . ~ . - p.fun.div.1) %>% AIC

mod_f_path_plot_sel <- update(mod_f_path_plot_sel, . ~ . - p.fun.div.0)
performance::check_collinearity(mod_f_path_plot_sel)

update(mod_f_path_plot_sel, . ~ . - i.tax.abund.det) %>% AIC
update(mod_f_path_plot_sel, . ~ . - i.tax.abund.sort) %>% AIC

mod_f_path_plot_sel <- update(mod_f_path_plot_sel, . ~ . - i.tax.abund.sort)
performance::check_collinearity(mod_f_path_plot_sel)

update(mod_f_path_plot_sel, . ~ . - p.tax.div.0) %>% AIC
update(mod_f_path_plot_sel, . ~ . - p.tax.div.2) %>% AIC

mod_f_path_plot_sel <- update(mod_f_path_plot_sel, . ~ . - p.tax.div.0)
performance::check_collinearity(mod_f_path_plot_sel)

update(mod_f_path_plot_sel, . ~ . - p.fun.comp.sla) %>% AIC
update(mod_f_path_plot_sel, . ~ . - p.fun.comp.ldmc) %>% AIC

mod_f_path_plot_sel <- update(mod_f_path_plot_sel, . ~ . - p.fun.comp.ldmc)
performance::check_collinearity(mod_f_path_plot_sel)

summary(mod_f_path_plot_sel)

# SEM

l_modlist_f_plot_1 <- list(
  lm(p.fgr.comp.forb ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fun.comp.sla ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(p.tax.div.2 ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fun.div.1 ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(i.tax.abund.det ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(i.fun.abund ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(i.tax.div.0 ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  lm(i.fun.div.2 ~ ForMI + region, 
     data = d_f_path_plotlevel_z),
  mod_f_path_plot_sel
)

# perform stepwise SEM adaption
l_modlist_f_plot_sel <- f.SEM.step(l_modlist_f_plot_1, d_f_path_plotlevel_z)


# summarise pathways

d_semcoefs_f_plot <- sem.coefs(l_modlist_f_plot_sel, d_f_path_plotlevel_z,
                               standardize = "none")

names(d_semcoefs_f_plot)[6] <- "sign.code"

d_semcoefs_f_plot <- d_semcoefs_f_plot %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate_at(vars(response, predictor), ~gsub("~~ ", "", .))

# add significance
d_semcoefs_f_plot <- d_semcoefs_f_plot %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_f_plot <- data.frame()

vars <- unique(d_semcoefs_f_plot$response)
vars <- vars[!vars %in% c("ForMI")]

for (var_i in c("ForMI", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_f_plot %>% 
    filter(predictor == var_i) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop.comm_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop.comm_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_f_plot,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop.comm_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop.comm_log")
  }
  
  
  pathways_f_plot <- pathways_f_plot %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
}


pathways_f_plot_sum <- pathways_f_plot %>% 
  filter(predictor %in% c("ForMI")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()


# ForMI split ------------------------------------------------------------------.

# define global model
mod_f_path_plot_sep <- lm(a.herb2.prop.comm_log ~ Inonat + Iharv + Idwcut + region +
                            p.tax.comp.1 + p.tax.comp.2 + p.tax.comp.3 + 
                            p.fgr.comp.grass +  p.fgr.comp.forb + 
                            p.fgr.comp.tree + p.fgr.comp.geophyt + 
                            p.fun.comp.sla + p.fun.comp.ldmc +
                            p.tax.abund +
                            p.fun.abund.lai + p.fun.abund.bm +
                            p.tax.div.0 + p.tax.div.1 + p.tax.div.2 +
                            p.fun.div.0 + p.fun.div.1 + p.fun.div.2 +
                            i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                            i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                            i.tax.abund.sort + i.tax.abund.det +
                            i.fun.abund +
                            i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                            i.fun.div.0 + i.fun.div.1 + i.fun.div.2,
                          data = d_f_path_plotlevel_z)


# perform model selection
mod_f_path_plot_sep_sel <- f.step.adj(mod_f_path_plot_sep, c("Inonat", "Idwcut", "Iharv",
                                                             "region"))


# manually check and adjust for collinear variables 
performance::check_collinearity(mod_f_path_plot_sep_sel)

update(mod_f_path_plot_sep_sel, . ~ . - p.tax.div.0) %>% AIC
update(mod_f_path_plot_sep_sel, . ~ . - p.tax.div.1) %>% AIC

mod_f_path_plot_sep_sel <- update(mod_f_path_plot_sep_sel, . ~ . - p.tax.div.1)
performance::check_collinearity(mod_f_path_plot_sep_sel)

update(mod_f_path_plot_sep_sel, . ~ . - i.tax.abund.det) %>% AIC
update(mod_f_path_plot_sep_sel, . ~ . - i.tax.abund.sort) %>% AIC

mod_f_path_plot_sep_sel <- update(mod_f_path_plot_sep_sel, . ~ . - i.tax.abund.sort)
performance::check_collinearity(mod_f_path_plot_sep_sel)

update(mod_f_path_plot_sep_sel, . ~ . - p.fun.comp.sla) %>% AIC
update(mod_f_path_plot_sep_sel, . ~ . - p.fun.comp.ldmc) %>% AIC

mod_f_path_plot_sep_sel <- update(mod_f_path_plot_sep_sel, . ~ . - p.fun.comp.ldmc)
performance::check_collinearity(mod_f_path_plot_sep_sel)

# SEM

l_modlist_f_plot_sep_1 <- list(
  lm(p.fgr.comp.forb ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fgr.comp.grass ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fgr.comp.geophyt ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fun.comp.sla ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.tax.abund ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.tax.div.0 ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(p.fun.div.2 ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(i.fun.comp.2 ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(i.tax.abund.det ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  lm(i.fun.abund ~ Inonat + Iharv + Idwcut + region, 
     data = d_f_path_plotlevel_z),
  mod_f_path_plot_sep_sel
)

# perform stepwise SEM adaption
l_modlist_f_plot_sep_sel <- f.SEM.step(l_modlist_f_plot_sep_1, d_f_path_plotlevel_z, 
                                       corr.errors = c("Inonat ~~ Iharv", 
                                                       "Inonat ~~ Idwcut", 
                                                       "Idwcut ~~ Iharv"))



# summarise pathways

d_semcoefs_f_plot_sep <- sem.coefs(l_modlist_f_plot_sep_sel, d_f_path_plotlevel_z,
                                   standardize = "none",
                                   corr.errors = c("Inonat ~~ Iharv", 
                                                   "Inonat ~~ Idwcut", 
                                                   "Idwcut ~~ Iharv"))

names(d_semcoefs_f_plot_sep)[6] <- "sign.code"

d_semcoefs_f_plot_sep <- d_semcoefs_f_plot_sep %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate_at(vars(response, predictor), ~gsub("~~ ", "", .))

# add significance
d_semcoefs_f_plot_sep <- d_semcoefs_f_plot_sep %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_f_plot_sep <- data.frame()

vars <- unique(d_semcoefs_f_plot_sep$response)
vars <- vars[!vars %in% c("Inonat", "Idwcut", "Iharv")]


for (var_i in c("Inonat", "Idwcut", "Iharv", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_f_plot_sep %>% 
    filter(predictor == var_i) %>% 
    filter(!response %in% c("Inonat", "Idwcut", "Iharv")) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop.comm_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop.comm_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_f_plot_sep,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop.comm_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop.comm_log")
  }
  
  pathways_f_plot_sep <- pathways_f_plot_sep %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}

pathways_f_plot_sep_sum <- pathways_f_plot_sep %>% 
  filter(predictor %in% c("Inonat", "Iharv", "Idwcut")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# grasslands community level ---------------------------------------------------

# prepare scaled data set
d_g_path_plotlevel_z <- d_comm %>% 
  select(plot, region, RW, HW, 
         a.herb2.prop.comm_log, LUI, G_std, M_std, F_std, 
         p.tax.comp.1, p.tax.comp.2, p.tax.comp.3, 
         p.fgr.comp.grass, p.tax.comp.legume,
         p.fun.comp.sla, p.fun.comp.ldmc,
         p.fun.comp.N, p.fun.comp.P, p.fun.comp.Ca, p.fun.comp.Mg,
         p.fun.comp.lignin, p.fun.comp.prim.fiber,
         p.tax.abund,
         p.fun.abund,
         p.tax.div.0, p.tax.div.1, p.tax.div.2,
         p.fun.div.0, p.fun.div.1, p.fun.div.2,
         i.tax.comp.1, i.tax.comp.2, i.tax.comp.3,
         i.fun.comp.1, i.fun.comp.2, i.fun.comp.3,
         i.tax.abund.sort, i.tax.abund.det,
         i.fun.abund,
         i.tax.div.0, i.tax.div.1, i.tax.div.2,
         i.fun.div.0, i.fun.div.1, i.fun.div.2) %>% 
  filter_all(~!is.na(.)) %>%  # remove all entries with NAs in any of the variables
  mutate_at(vars(-c(plot, region, RW, HW, a.herb2.prop.comm_log)), 
            ~ f_std(.)) 

# LUI combined -----------------------------------------------------------------.

# define global model
mod_g_path_plot <- lm(a.herb2.prop.comm_log ~ LUI + region +
                        p.tax.comp.1 + p.tax.comp.2 + p.tax.comp.3 + 
                        p.fgr.comp.grass + p.tax.comp.legume +
                        p.fun.comp.sla + p.fun.comp.ldmc +
                        p.fun.comp.N + p.fun.comp.P + p.fun.comp.Ca + p.fun.comp.Mg +
                        p.fun.comp.lignin + p.fun.comp.prim.fiber +
                        p.tax.abund +
                        p.fun.abund +
                        p.tax.div.0 + p.tax.div.1 + p.tax.div.2 +
                        p.fun.div.0 + p.fun.div.1 + p.fun.div.2 +
                        i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                        i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                        i.tax.abund.sort + i.tax.abund.det +
                        i.fun.abund +
                        i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                        i.fun.div.0 + i.fun.div.1 + i.fun.div.2,
                      data = d_g_path_plotlevel_z)

# perform model selection
mod_g_path_plot_sel <- f.step.adj(mod_g_path_plot, c("LUI", "region"))

# manually check and adjust for collinear variables 
performance::check_collinearity(mod_g_path_plot_sel)

update(mod_g_path_plot_sel, . ~ . - p.fun.div.0 - p.fun.div.1) %>% AIC
update(mod_g_path_plot_sel, . ~ . - p.fun.div.0 - p.fun.div.2) %>% AIC
update(mod_g_path_plot_sel, . ~ . - p.fun.div.1 - p.fun.div.2) %>% AIC

mod_g_path_plot_sel <- update(mod_g_path_plot_sel, . ~ . - p.fun.div.1 - p.fun.div.2)
performance::check_collinearity(mod_g_path_plot_sel)

update(mod_g_path_plot_sel, . ~ . - i.tax.abund.det) %>% AIC
update(mod_g_path_plot_sel, . ~ . - i.tax.abund.sort) %>% AIC

mod_g_path_plot_sel <- update(mod_g_path_plot_sel, . ~ . - i.tax.abund.sort)
performance::check_collinearity(mod_g_path_plot_sel)

update(mod_g_path_plot_sel, . ~ . - p.tax.div.2) %>% AIC
update(mod_g_path_plot_sel, . ~ . - p.tax.div.0) %>% AIC

mod_g_path_plot_sel <- update(mod_g_path_plot_sel, . ~ . - p.tax.div.2)
performance::check_collinearity(mod_g_path_plot_sel)

# SEM

l_modlist_g_plot_1 <- list(
  lm(p.tax.comp.2 ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fgr.comp.grass ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fgr.comp.legume ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.comp.lignin ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.comp.prim.fiber ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.tax.abund ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.tax.div.0 ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.div.0 ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(i.tax.comp.2 ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  lm(i.tax.abund.det ~ LUI + region, 
     data = d_g_path_plotlevel_z),
  mod_g_path_plot_sel
)

# perform stepwise SEM adaption
l_modlist_g_plot_sel <- f.SEM.step(l_modlist_g_plot_1, d_g_path_plotlevel_z)



# summarise pathways

d_semcoefs_g_plot <- sem.coefs(l_modlist_g_plot_sel, d_g_path_plotlevel_z,
                               standardize = "none")

names(d_semcoefs_g_plot)[6] <- "sign.code"

# add significance
d_semcoefs_g_plot <- d_semcoefs_g_plot %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_g_plot <- data.frame()

vars <- unique(d_semcoefs_g_plot$response)
vars <- vars[!vars %in% c("LUI")]

for (var_i in c("LUI", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_g_plot %>% 
    filter(predictor == var_i) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop.comm_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop.comm_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_g_plot,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop.comm_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop.comm_log")
  }
  
  pathways_g_plot <- pathways_g_plot %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}

pathways_g_plot_sum <- pathways_g_plot %>% 
  filter(predictor == "LUI") %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# LUI split --------------------------------------------------------------------.

# define global model
mod_g_path_plot_sep <- lm(a.herb2.prop.comm_log ~ M_std + F_std + G_std + region +
                            p.tax.comp.1 + p.tax.comp.2 + p.tax.comp.3 + 
                            p.fgr.comp.grass + p.tax.comp.legume +
                            p.fun.comp.sla + p.fun.comp.ldmc +
                            p.fun.comp.N + p.fun.comp.P + p.fun.comp.Ca + p.fun.comp.Mg +
                            p.fun.comp.lignin + p.fun.comp.prim.fiber +
                            p.tax.abund +
                            p.fun.abund +
                            p.tax.div.0 + p.tax.div.1 + p.tax.div.2 +
                            p.fun.div.0 + p.fun.div.1 + p.fun.div.2 +
                            i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                            i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                            i.tax.abund.sort + i.tax.abund.det +
                            i.fun.abund +
                            i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                            i.fun.div.0 + i.fun.div.1 + i.fun.div.2,
                          data = d_g_path_plotlevel_z)


# perform model selection
mod_g_path_plot_sep_sel <- f.step.adj(mod_g_path_plot_sep, c("M_std", "F_std", "G_std",
                                                             "region"))


# manually check and adjust for collinear variables
performance::check_collinearity(mod_g_path_plot_sep_sel)

update(mod_g_path_plot_sep_sel, . ~ . - i.fun.div.1 -  i.fun.div.0) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - i.fun.div.2 -  i.fun.div.0) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - i.fun.div.2 -  i.fun.div.1) %>% AIC

mod_g_path_plot_sep_sel <- update(mod_g_path_plot_sep_sel, . ~ . - i.fun.div.1 - i.fun.div.0)
performance::check_collinearity(mod_g_path_plot_sep_sel)

update(mod_g_path_plot_sep_sel, . ~ . - p.fun.div.1 -  p.fun.div.0) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - p.fun.div.2 -  p.fun.div.0) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - p.fun.div.2 -  p.fun.div.1) %>% AIC

mod_g_path_plot_sep_sel <- update(mod_g_path_plot_sep_sel, . ~ . - p.fun.div.2 - p.fun.div.1)
performance::check_collinearity(mod_g_path_plot_sep_sel)

update(mod_g_path_plot_sep_sel, . ~ . - i.tax.abund.det) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - i.tax.abund.sort) %>% AIC

mod_g_path_plot_sep_sel <- update(mod_g_path_plot_sep_sel, . ~ . - i.tax.abund.sort)
performance::check_collinearity(mod_g_path_plot_sep_sel)

update(mod_g_path_plot_sep_sel, . ~ . - p.tax.div.2) %>% AIC
update(mod_g_path_plot_sep_sel, . ~ . - p.tax.div.0) %>% AIC

mod_g_path_plot_sep_sel <- update(mod_g_path_plot_sep_sel, . ~ . - p.tax.div.0)
performance::check_collinearity(mod_g_path_plot_sep_sel)

# SEM

l_modlist_g_plot_sep_1 <- list(
  lm(p.tax.comp.2 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fgr.comp.grass ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fgr.comp.legume ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.comp.lignin ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.comp.prim.fiber ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.tax.abund ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.tax.div.2 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(p.fun.div.0 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.tax.comp.2 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.fun.comp.3 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.tax.abund.det ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.fun.abund ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.tax.div.1 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  lm(i.fun.div.2 ~ M_std + F_std + G_std + region, 
     data = d_g_path_plotlevel_z),
  mod_g_path_plot_sep_sel
)

# perform stepwise SEM adaption
l_modlist_g_plot_sep_sel <- f.SEM.step(l_modlist_g_plot_sep_1, d_g_path_plotlevel_z, 
                                       corr.errors = c("M_std ~~ F_std", 
                                                       "M_std ~~ G_std",
                                                       "F_std ~~ G_std"))


# summarise pathways

d_semcoefs_g_plot_sep <- sem.coefs(l_modlist_g_plot_sep_sel, d_g_path_plotlevel_z,
                                   standardize = "none",
                                   corr.errors = c("M_std ~~ F_std", 
                                                   "M_std ~~ G_std",
                                                   "F_std ~~ G_std"))

names(d_semcoefs_g_plot_sep)[6] <- "sign.code"

d_semcoefs_g_plot_sep <- d_semcoefs_g_plot_sep %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate_at(vars(response, predictor), ~gsub("~~ ", "", .))


# add significance
d_semcoefs_g_plot_sep <- d_semcoefs_g_plot_sep %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

# var_i <- "F_std"
pathways_g_plot_sep <- data.frame()

vars <- unique(d_semcoefs_g_plot_sep$response)
vars <- vars[!vars %in% c("M_std", "F_std", "G_std")]

for (var_i in c("M_std", "F_std", "G_std", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_g_plot_sep %>% 
    filter(predictor == var_i) %>% 
    filter(!response %in% c("M_std", "F_std", "G_std")) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop.comm_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop.comm_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_g_plot_sep,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop.comm_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop.comm_log")
  }
  
  pathways_g_plot_sep <- pathways_g_plot_sep %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}

pathways_g_plot_sep_sum <- pathways_g_plot_sep %>% 
  filter(predictor %in% c("M_std", "F_std", "G_std")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# forests species level --------------------------------------------------------

d_f_path_plotlevel_z <- d_spec %>% 
  select(plot, region, RW, HW, Host_Name_std,
         a.herb2.prop_log, ForMI, Inonat, Iharv, Idwcut, 
         p.fgr.comp.grass, p.fgr.comp.forb,
         p.fgr.comp.tree, p.fgr.comp.geophyt,
         p.fun.comp.ldmc, p.fun.comp.sla,
         p.tax.abund, 
         p.fun.abund,
         i.tax.comp.iszero,
         i.tax.comp.1, i.tax.comp.2, i.tax.comp.3,
         i.fun.comp.1, i.fun.comp.2, i.fun.comp.3,
         i.tax.abund,
         i.fun.abund,
         i.tax.div.0, i.tax.div.1, i.tax.div.2,
         i.fun.div.0, i.fun.div.1, i.fun.div.2) %>% 
  filter_all(~!is.na(.)) %>%  # remove all entries with NAs in any of the variables
  mutate_at(vars(-c(a.herb2.prop_log, ForMI, Inonat, Iharv, Idwcut)),
            ~ f_std(.))

# ForMI combined ---------------------------------------------------------------.

# define global model
mod_f_path_sample <- glmmTMB(a.herb2.prop_log ~ ForMI + region +
                               p.fgr.comp.grass + p.fgr.comp.forb +
                               p.fgr.comp.tree + p.fgr.comp.geophyt +
                               p.tax.abund + 
                               p.fun.comp.ldmc + p.fun.comp.sla +
                               p.fun.abund +
                               i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                               i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                               i.tax.abund +
                               i.fun.div.0 + i.fun.div.1 + i.fun.div.2 +
                               i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                               i.fun.abund + i.tax.comp.iszero +
                               (1 | plot) + (1 | Host_Name_std),
                             REML = F, 
                             data = d_f_path_samplelevel_z)

# perform model selection
mod_f_path_sample_sel <- f.step.adj(mod_f_path_sample, c("ForMI", "region"),
                                    fixcomb = "i.tax.comp.iszero")

# manually check and adjust for collinear variables
performance::check_collinearity(mod_f_path_sample_sel)

# SEM

l_modlist_f_sample_1 <- list(
  glmmTMB(p.fgr.comp.forb ~ ForMI + region +
            (1 | plot), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(p.fgr.comp.tree ~ ForMI + region +
            (1 | plot), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(p.fgr.comp.geophyt ~ ForMI + region +
            (1 | plot), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(p.tax.abund ~ ForMI + region +
            (1 | plot) + (1 | Host_Name_std),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.iszero ~ ForMI + region +
            (1 | plot) + (1 | Host_Name_std), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.1 ~ ForMI + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.2 ~ ForMI + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.abund ~ ForMI + region +
            (1 | plot) + (1 | Host_Name_std),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.fun.abund ~ ForMI + region +
            (1 | plot) + (1 | Host_Name_std),
          data = d_f_path_samplelevel_z),
  mod_f_path_sample_sel
)

# perform stepwise SEM adaption
l_modlist_f_sample_sel <- f.SEM.step(l_modlist_f_sample_1, d_f_path_samplelevel_z,
                                     fixcomb = "i.tax.comp.iszero")


# summarise pathways

d_semcoefs_f_sample <- sem.coefs(l_modlist_f_sample_sel, d_f_path_samplelevel_z,
                                 standardize = "none")

names(d_semcoefs_f_sample)[6] <- "sign.code"

# adjust model coefficients for logistic regressions
d_semcoefs_f_sample <- f_adjust_semcoefs(d_semcoefs_f_sample, l_modlist_f_sample_sel)

d_semcoefs_f_sample <- d_semcoefs_f_sample %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response))

# add significance
d_semcoefs_f_sample <- d_semcoefs_f_sample %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_f_sample <- data.frame()

vars <- unique(d_semcoefs_f_sample$response)
vars <- vars[!vars %in% c("ForMI")]

for (var_i in c("ForMI", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_f_sample %>% 
    filter(predictor == var_i) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_f_sample,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop_log")
  }
  
  pathways_f_sample <- pathways_f_sample %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}

pathways_f_sample_sum <- pathways_f_sample %>% 
  filter(predictor %in% c("ForMI")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()


# ForMI split ------------------------------------------------------------------.

# define global model
mod_f_path_sample_sep <- glmmTMB(a.herb2.prop_log ~ Inonat + Iharv + Idwcut + region +
                                   p.fgr.comp.grass + p.fgr.comp.forb +
                                   p.fgr.comp.tree + p.fgr.comp.geophyt +
                                   p.tax.abund + 
                                   p.fun.comp.ldmc + p.fun.comp.sla +
                                   p.fun.abund +
                                   i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                                   i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                                   i.tax.abund +
                                   i.fun.div.0 + i.fun.div.1 + i.fun.div.2 +
                                   i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                                   i.fun.abund + i.tax.comp.iszero +
                                   (1 | plot) + (1 | Host_Name_std),
                                 REML = F, 
                                 data = d_f_path_samplelevel_z)

# perform model selection
mod_f_path_sample_sep_sel <- f.step.adj(mod_f_path_sample_sep, c("Inonat", "Iharv", 
                                                                 "Idwcut", "region"),
                                        fixcomb = "i.tax.comp.iszero")

# manually check and adjust for collinear variables
performance::check_collinearity(mod_f_path_sample_sep_sel)

# SEM

l_modlist_f_sample_sep_1 <- list(
  glmmTMB(p.fgr.comp.forb ~ Inonat + Iharv + Idwcut + region +
            (1 | plot), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(p.fgr.comp.tree ~ Inonat + Iharv + Idwcut + region +
            (1 | plot), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(p.tax.abund ~ Inonat + Iharv + Idwcut + region +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.iszero ~ Inonat + Iharv + Idwcut + region+
            (1 | plot) + (1 | Host_Name_std), 
          family = binomial(),
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.2 ~ Inonat + Iharv + Idwcut + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.comp.3 ~ Inonat + Iharv + Idwcut + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_f_path_samplelevel_z),
  glmmTMB(i.tax.abund ~ Inonat + Iharv + Idwcut + region +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_f_path_samplelevel_z),
  glmmTMB(i.fun.abund ~ Inonat + Iharv + Idwcut + region +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_f_path_samplelevel_z),
  mod_f_path_sample_sep_sel
)

# perform stepwise SEM adaption
l_modlist_f_sample_sep_sel <- f.SEM.step(l_modlist_f_sample_sep_1, data = d_f_path_samplelevel_z,
                                         corr.errors = c("Inonat ~~ Idwcut",
                                                         "Inonat ~~ Iharv",
                                                         "Iharv ~~ Idwcut"),
                                         fixcomb = "i.tax.comp.iszero")


# summarise pathways

d_semcoefs_f_sample_sep <- sem.coefs(l_modlist_f_sample_sep_sel, d_f_path_samplelevel_z,
                                     standardize = "none",
                                     corr.errors = c("Inonat ~~ Iharv", 
                                                     "Inonat ~~ Idwcut", 
                                                     "Idwcut ~~ Iharv"))

names(d_semcoefs_f_sample_sep)[6] <- "sign.code"

# adjust model coefficients for logistic regressions
d_semcoefs_f_sample_sep <- f_adjust_semcoefs(d_semcoefs_f_sample_sep, l_modlist_f_sample_sep_sel)

d_semcoefs_f_sample_sep <- d_semcoefs_f_sample_sep %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate_at(vars(response, predictor), ~gsub("~~ ", "", .)) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response))

# add significance
d_semcoefs_f_sample_sep <- d_semcoefs_f_sample_sep %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

# var_i <- "F_std"
pathways_f_sample_sep <- data.frame()

vars <- unique(d_semcoefs_f_sample_sep$response)
vars <- vars[!vars %in% c("Inonat", "Iharv", "Idwcut")]

for (var_i in c("Inonat", "Iharv", "Idwcut", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_f_sample_sep %>% 
    filter(predictor == var_i) %>% 
    filter(!response %in% c("Inonat", "Iharv", "Idwcut")) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_f_sample_sep,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop_log")
  }
  
  pathways_f_sample_sep <- pathways_f_sample_sep %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}

pathways_f_sample_sep_sum <- pathways_f_sample_sep %>% 
  filter(predictor %in% c("Inonat", "Iharv", "Idwcut")) %>%
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# grasslands species level -----------------------------------------------------

d_f_path_plotlevel_z <- d_spec %>% 
  select(plot, region, RW, HW, Host_Name_std,
         a.herb2.prop_log, LUI, M_std, F_std, G_std, 
         p.fgr.comp.grass, p.tax.comp.legume,
         p.fun.comp.ldmc, p.fun.comp.sla,
         p.tax.abund, 
         p.fun.abund,
         i.tax.comp.iszero,
         i.tax.comp.1, i.tax.comp.2, i.tax.comp.3,
         i.fun.comp.1, i.fun.comp.2, i.fun.comp.3,
         i.tax.abund,
         i.fun.abund,
         i.fun.div.0, i.fun.div.1, i.fun.div.2,
         i.tax.div.0, i.tax.div.1, i.tax.div.2,) %>% 
  filter_all(~!is.na(.)) %>%  # remove all entries with NAs in any of the variables
  mutate_at(vars(-c(a.herb2.prop_log, LUI, M_std, F_std, G_std)),
            ~ f_std(.))

# LUI combined -----------------------------------------------------------------.

# define global model
mod_g_path_sample <- glmmTMB(a.herb2.prop_log ~ LUI + region +
                               p.fgr.comp.grass + p.tax.comp.legume +
                               p.tax.abund + 
                               p.fun.comp.ldmc + p.fun.comp.sla +
                               p.fun.abund +
                               i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                               i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                               i.tax.abund +
                               i.fun.div.0 + i.fun.div.1 + i.fun.div.2 +
                               i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                               i.fun.abund + i.tax.comp.iszero +
                               (1 | plot) + (1 | Host_Name_std),
                             REML = F,
                             data = d_g_path_samplelevel_z)

# perform model selection
mod_g_path_sample_sel <- f.step.adj(mod_g_path_sample, c("LUI", "region"), 
                                    fixcomb = "i.tax.comp.iszero")

# manually check and adjust for collinear variables
performance::check_collinearity(mod_g_path_sample_sel)

update(mod_g_path_sample_sel, ". ~ . - i.tax.div.0") %>% AIC
update(mod_g_path_sample_sel, ". ~ . - i.tax.div.1") %>% AIC

mod_g_path_sample_sel <- update(mod_g_path_sample_sel, ". ~ . -  i.tax.div.0")
performance::check_collinearity(mod_g_path_sample_sel)

update(mod_g_path_sample_sel, ". ~ . - p.fun.comp.ldmc") %>% AIC
update(mod_g_path_sample_sel, ". ~ . - p.fun.comp.sla") %>% AIC

mod_g_path_sample_sel <- update(mod_g_path_sample_sel, ". ~ . - p.fun.comp.ldmc")
performance::check_collinearity(mod_g_path_sample_sel)


# SEM

l_modlist_g_sample_1 <- list(
  glmmTMB(p.fgr.comp.grass ~ LUI + region +
            (1 | plot), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(p.fgr.comp.legume ~ LUI + region +
            (1 | plot), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(p.fun.comp.sla ~ LUI + region + 
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  glmmTMB(i.tax.comp.iszero ~ LUI + region +
            (1 | plot) + (1 | Host_Name_std), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(i.fun.comp.2 ~ LUI + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  glmmTMB(i.fun.comp.3 ~ LUI + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  glmmTMB(i.tax.abund ~ LUI + region + 
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  glmmTMB(i.tax.div.1 ~ LUI + region + 
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  mod_g_path_sample_sel
)

# perform stepwise SEM adaption
l_modlist_g_sample_sel <- f.SEM.step(l_modlist_g_sample_1, d_g_path_samplelevel_z)


# summarise pathways

d_semcoefs_g_sample <- sem.coefs(l_modlist_g_sample_sel, d_g_path_samplelevel_z,
                                 standardize = "none")

names(d_semcoefs_g_sample)[6] <- "sign.code"

# adjust model coefficients for logistic regressions
d_semcoefs_g_sample <- f_adjust_semcoefs(d_semcoefs_g_sample, l_modlist_g_sample_sel)

d_semcoefs_g_sample <- d_semcoefs_g_sample %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response))

# add significance
d_semcoefs_g_sample <- d_semcoefs_g_sample %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_g_sample <- data.frame()

vars <- unique(d_semcoefs_g_sample$response)
vars <- vars[!vars %in% c("LUI")]

for (var_i in c("LUI", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_g_sample %>% 
    filter(predictor == var_i) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_g_sample,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop_log")
  }
  
  pathways_g_sample <- pathways_g_sample %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
  
}
pathways_g_sample_sum <- pathways_g_sample %>% 
  filter(predictor %in% c("LUI")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# LUI split --------------------------------------------------------------------.

# define global model
mod_g_path_sample_sep <- glmmTMB(a.herb2.prop_log ~ G_std + M_std + F_std + region +
                                   p.fgr.comp.grass + p.tax.comp.legume +
                                   p.tax.abund + 
                                   p.fun.comp.ldmc + p.fun.comp.sla +
                                   p.fun.abund +
                                   i.tax.div.0 + i.tax.div.1 + i.tax.div.2 +
                                   i.tax.comp.1 + i.tax.comp.2 + i.tax.comp.3 +
                                   i.tax.abund +
                                   i.fun.div.0 + i.fun.div.1 + i.fun.div.2 +
                                   i.fun.comp.1 + i.fun.comp.2 + i.fun.comp.3 +
                                   i.fun.abund + i.tax.comp.iszero +
                                   (1 | plot) + (1 | Host_Name_std),
                                 REML = F, 
                                 data = d_g_path_samplelevel_z)

# perform model selection
mod_g_path_sample_sep_sel <- f.step.adj(mod_g_path_sample_sep, c("G_std", "M_std", 
                                                                 "F_std", "region"),
                                        fixcomb = "i.tax.comp.iszero")

# manually check and adjust for collinear variables
performance::check_collinearity(mod_g_path_sample_sep_sel)

update(mod_g_path_sample_sep_sel, . ~ . - p.fun.comp.ldmc) %>% AIC
update(mod_g_path_sample_sep_sel, . ~ . - p.fun.comp.sla) %>% AIC

mod_g_path_sample_sep_sel <- update(mod_g_path_sample_sep_sel, . ~ . - p.fun.comp.ldmc)
performance::check_collinearity(mod_g_path_sample_sep_sel)

# SEM

l_modlist_g_sample_sep_1 <- list(
  glmmTMB(p.fgr.comp.grass ~ M_std + F_std + G_std + region +
            (1 | plot), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(p.fgr.comp.legume ~ M_std + F_std + G_std + region +
            (1 | plot), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(p.fun.comp.sla ~ M_std + F_std + G_std + region + 
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  glmmTMB(i.tax.comp.iszero ~ M_std + F_std + G_std + region +
            (1 | plot) + (1 | Host_Name_std), 
          family = binomial(),
          data = d_g_path_samplelevel_z),
  glmmTMB(i.fun.comp.2 ~ M_std + F_std + G_std + region + i.tax.comp.iszero +
            (1 | plot) + (1 | Host_Name_std), 
          data = d_g_path_samplelevel_z),
  mod_g_path_sample_sep_sel
)

# perform stepwise SEM adaption
l_modlist_g_sample_sep_sel <- f.SEM.step(l_modlist_g_sample_sep_1, data = d_g_path_samplelevel_z,
                                         corr.errors = c("M_std ~~ F_std",
                                                         "M_std ~~ G_std",
                                                         "F_std ~~ G_std"),
                                         fixcomb = "i.tax.comp.iszero")


# summarise pathways

d_semcoefs_g_sample_sep <- sem.coefs(l_modlist_g_sample_sep_sel, d_g_path_samplelevel_z,
                                     standardize = "none",
                                     corr.errors = c("M_std ~~ F_std", 
                                                     "M_std ~~ G_std",
                                                     "F_std ~~ G_std"))

names(d_semcoefs_g_sample_sep)[6] <- "sign.code"

# adjust model coefficients for logistic regressions
d_semcoefs_g_sample_sep <- f_adjust_semcoefs(d_semcoefs_g_sample_sep, l_modlist_g_sample_sep_sel)

d_semcoefs_g_sample_sep <- d_semcoefs_g_sample_sep %>% 
  mutate_if(is.factor, ~ as.character(.)) %>% 
  mutate_at(vars(response, predictor), ~gsub("~~ ", "", .)) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response))

# add significance
d_semcoefs_g_sample_sep <- d_semcoefs_g_sample_sep %>% 
  mutate(sign = ifelse(p.value <= 0.05, T, F))

pathways_g_sample_sep <- data.frame()

vars <- unique(d_semcoefs_g_sample_sep$response)
vars <- vars[!vars %in% c("M_std", "F_std", "G_std")]

for (var_i in c("M_std", "F_std", "G_std", "regionHAI", "regionSCH", vars)){
  pathways_var_i <- data.frame()
  
  
  effects_i <- d_semcoefs_g_sample_sep %>% 
    filter(predictor == var_i) %>% 
    filter(!response %in% c("M_std", "F_std", "G_std")) %>% 
    mutate(path = paste(predictor, response, sep = " > "))
  
  
  pathways_var_i <- effects_i %>% 
    filter(response == "a.herb2.prop_log") %>% 
    bind_rows(pathways_var_i, .)
  
  effects_i <- effects_i %>% 
    filter(response != "a.herb2.prop_log")
  
  while (nrow(effects_i) > 0){
    effects_i <- effects_i %>% 
      left_join(d_semcoefs_g_sample_sep,
                by = c(response = "predictor"),
                suffix = c(".a", ".b")) %>% 
      rowwise() %>% 
      mutate(estimate = estimate.a * estimate.b,
             sign = all(c(sign.a, sign.b))) %>% 
      filter(! grepl(response.b, path)) %>% # exclude reciprocal circles
      ungroup() %>% 
      select(-c(response, estimate.a, estimate.b, sign.a, sign.b)) %>% 
      rename(response = response.b) %>% 
      mutate(path = paste(path, response, sep = " > "))
    
    pathways_var_i <- effects_i %>% 
      filter(response == "a.herb2.prop_log") %>% 
      bind_rows(pathways_var_i, .)
    
    effects_i <- effects_i %>% 
      filter(response != "a.herb2.prop_log")
  }
  
  pathways_g_sample_sep <- pathways_g_sample_sep %>% 
    bind_rows(pathways_var_i %>% 
                select(predictor, response, estimate, sign, path))
  
}


pathways_g_sample_sep_sum <- pathways_g_sample_sep %>% 
  filter(predictor %in% c("M_std", "F_std", "G_std")) %>% 
  mutate(groups = ifelse(grepl("> p.", path) &
                           grepl("> i.", path),
                         "both", 
                         ifelse(grepl("> p.", path), "plants", 
                                ifelse(grepl("> i.", path), "insects", 
                                       "direct")))) %>% 
  group_by(predictor, groups) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup()

# plotting (Fig. 4) ------------------------------------------------------------

pathways_sum <- pathways_g_plot_sum %>% 
  mutate(system = "Grassland",
         level = "Community level") %>% 
  bind_rows(pathways_g_plot_sep_sum %>% 
              mutate(system = "Grassland",
                     level = "Community level")) %>% 
  bind_rows(pathways_g_sample_sum %>% 
              mutate(system = "Grassland",
                     level = "Species level")) %>% 
  bind_rows(pathways_g_sample_sep_sum %>% 
              mutate(system = "Grassland",
                     level = "Species level")) %>% 
  bind_rows(pathways_f_plot_sum %>% 
              mutate(system = "Forest",
                     level = "Community level")) %>% 
  bind_rows(pathways_f_plot_sep_sum %>% 
              mutate(system = "Forest",
                     level = "Community level")) %>% 
  bind_rows(pathways_f_sample_sum %>% 
              mutate(system = "Forest",
                     level = "Species level")) %>% 
  bind_rows(pathways_f_sample_sep_sum %>% 
              mutate(system = "Forest",
                     level = "Species level")) %>% 
  mutate(system = factor(system, levels = c("Forest", "Grassland")),
         level = factor(level, levels = c("Community level", "Species level")),
         predictor = factor(predictor, levels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                                      "LUI", "G_std", "M_std", "F_std")),
                            labels = rev(c("ForMI", "Inonat", "Iharv", "Idwcut",
                                           "LUI", "Grazing", "Mowing", "Fertilisation"))),
         groups = factor(groups, levels = rev(c("direct", "plants", "both", "insects")),
                         labels = rev(c("Direct", "Plants", "Both", "Herbivores"))))


pathways_sum %>% 
  ggplot(aes(x = predictor, y = estimate, fill = groups)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 3.5, lty = 1) +
  geom_vline(xintercept = c(1.5, 2.5), lty = 3) +
  coord_flip() +
  facet_wrap(level ~ system, scales = "free_y") +
  scale_fill_manual(values = c(Direct = "#FFD966",
                               Plants = "#A9D18E",
                               Both = "#96755A",
                               Herbivores = "#F4B183"),
                    guide = guide_legend(reverse = TRUE),
                    name = "Pathway")

# plotting (Fig. 5 and supplementary) ------------------------------------------

# function to extract indirect pathways
f_ind_pathways <- function(data, luivars){
  covars <- unique(data$predictor)
  covars <- covars[!covars %in% c("regionHAI", "regionSCH", luivars)]
  
  out <- data.frame()
  
  for (var_i in luivars){
    for (covar_i in covars){
      out <- data %>% 
        filter(predictor == var_i & 
                 grepl(covar_i, path)) %>% 
        summarise(estimate = sum(estimate)) %>% 
        mutate(var = var_i,
               covar = covar_i) %>% 
        bind_rows(out, .)
      
    }
  }
  out
}

# combined ---------------------------------------------------------------------.

pathways_indirect <- f_ind_pathways(pathways_g_plot, c("LUI")) %>% 
  mutate(system = "Grassland",
         level = "Community level") %>% 
  bind_rows(f_ind_pathways(pathways_g_sample, c("LUI")) %>% 
              mutate(system = "Grassland",
                     level = "Species level")) %>% 
  bind_rows(f_ind_pathways(pathways_f_plot, c("ForMI")) %>% 
              mutate(system = "Forest",
                     level = "Community level")) %>% 
  bind_rows(f_ind_pathways(pathways_f_sample, c("ForMI")) %>% 
              mutate(system = "Forest",
                     level = "Species level")) %>% 
  mutate(trophy = ifelse(substr(covar, 1, 2) == "p.", "Plants", "Herbivores"),
         covar = factor(covar, levels = names(varlabs)))

p_indirect <- list()
for (level_i in c("Community level", "Species level")){
  for (system_i in c("Forest", "Grassland")){
    textcol_i <- pathways_indirect %>% 
      filter(system == system_i,
             level == level_i) %>% 
      select(covar) %>% 
      distinct() %>% 
      deframe()
    
    textcol_i <- ifelse(substr(textcol_i, 1, 1) == "p", "#70AD47", "#ED7D31")
    
    p_indirect[[paste(system_i, level_i)]] <-
      pathways_indirect %>% 
      filter(system == system_i,
             level == level_i) %>% 
      ggplot(aes(x = covar, y = estimate, fill = trophy)) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_col() +
      facet_grid(~  var, scales = "free_y") +
      coord_flip() +
      scale_fill_manual(values = c(Plants = "#A9D18E",
                                   Herbivores = "#F4B183")) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(color = textcol_i),
            strip.text = element_text(margin = margin(2,0,2,0, "mm"))) +
      scale_x_discrete(labels = varlabs)
  }
}

p_indirect[c("Forest Community level", "Grassland Community level")] <-
  lapply(p_indirect[c("Forest Community level", "Grassland Community level")],
         function(x) x + 
           theme(plot.margin = margin(40, 2, 22, 2, "pt")))

p_indirect[c("Forest Species level", "Grassland Species level")] <-
  lapply(p_indirect[c("Forest Species level", "Grassland Species level")],
         function(x) x + 
           theme(plot.margin = margin(40, 2, 2, 2, "pt")))


plot_grid(plot_grid(plotlist = p_indirect, labels = letters[1:4],
                    rel_heights = c(4.2, 4)),
          ggplot_labeller("Summed pathways estimates", 0, 14, NA),
          ncol = 1, rel_heights = c(20, 1))

# split ------------------------------------------------------------------------.

pathways_indirect_sep <- f_ind_pathways(pathways_g_plot_sep, c("M_std", "F_std", "G_std")) %>% 
  mutate(system = "Grassland",
         level = "Community level") %>% 
  bind_rows(f_ind_pathways(pathways_g_sample_sep, c("M_std", "F_std", "G_std")) %>% 
              mutate(system = "Grassland",
                     level = "Species level")) %>% 
  bind_rows(f_ind_pathways(pathways_f_plot_sep, c("Inonat", "Iharv", "Idwcut")) %>% 
              mutate(system = "Forest",
                     level = "Community level")) %>% 
  bind_rows(f_ind_pathways(pathways_f_sample_sep, c("Inonat", "Iharv", "Idwcut")) %>% 
              mutate(system = "Forest",
                     level = "Species level")) %>% 
  mutate(trophy = ifelse(substr(covar, 1, 2) == "p.", "Plants", "Herbivores"),
         covar = factor(covar, levels = names(varlabs)),
         var = factor(var, levels = c("Inonat", "Iharv", "Idwcut", 
                                      "G_std", "M_std", "F_std"),
                      labels = c("Inonat", "Iharv", "Idwcut",
                                 "Grazing", "Mowing", "Fertilisation")))

p_indirect <- list()
for (level_i in c("Community level", "Species level")){
  for (system_i in c("Forest", "Grassland")){
    textcol_i <- pathways_indirect_sep %>% 
      filter(system == system_i,
             level == level_i) %>% 
      select(covar) %>% 
      distinct() %>% 
      deframe()
    
    textcol_i <- ifelse(substr(textcol_i, 1, 1) == "p", "#70AD47", "#ED7D31")
    
    p_indirect[[paste(system_i, level_i)]] <-
      pathways_indirect_sep %>% 
      filter(system == system_i,
             level == level_i) %>% 
      ggplot(aes(x = covar, y = estimate, fill = trophy)) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_col() +
      facet_grid(~  var, scales = "free_y") +
      coord_flip() +
      scale_fill_manual(values = c(Plants = "#A9D18E",
                                   Herbivores = "#F4B183")) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(color = textcol_i),
            strip.text = element_text(margin = margin(2,0,2,0, "mm"))) +
      scale_x_discrete(labels = varlabs)
  }
}

p_indirect[c("Forest Community level", "Grassland Community level")] <-
  lapply(p_indirect[c("Forest Community level", "Grassland Community level")],
         function(x) x + 
           theme(plot.margin = margin(40, 2, 22, 2, "pt")))

p_indirect[c("Forest Species level", "Grassland Species level")] <-
  lapply(p_indirect[c("Forest Species level", "Grassland Species level")],
         function(x) x + 
           theme(plot.margin = margin(40, 2, 2, 2, "pt")))


plot_grid(plot_grid(plotlist = p_indirect, labels = letters[1:4],
                    rel_heights = c(4.2, 4)),
          ggplot_labeller("Summed pathways estimates", 0, 14, NA),
          ncol = 1, rel_heights = c(20, 1))

# detailled SEM plots ----------------------------------------------------------

# Forest community level combined ----------------------------------------------

# create frame for plotting, add rows of potentially possible connections:
d_semcoefs_f_plot2 <- expand.grid(predictor = unique(d_semcoefs_f_plot$predictor),
                                  response = unique(d_semcoefs_f_plot$response)) %>% 
  mutate_all(~as.character(.)) %>% 
  left_join(d_semcoefs_f_plot, by = c("predictor", "response")) %>% 
  mutate(exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude)) %>% 
  filter(!exclude)


var_order <- sapply(l_modlist_f_plot_sel,
                    function(x) as.character(formula(x)[2]))

d_R2_f_plot <- lapply(l_modlist_f_plot_sel, 
                      function(x) data.frame(response = as.character(formula(x)[2]),
                                             R2 = rsquared(x)$R.squared,
                                             stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")))


d_semcoefs_f_plot2 <- d_semcoefs_f_plot2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_f_plot) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "ForMI",
                                                  var_order, "R2")),
         response = factor(response, levels = c("ForMI", var_order)))


in_direct_paths_f_plot <- pathways_f_plot %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_f_plot2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_f_plot2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[2] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_f_plot2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_f_plot %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = .5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - .5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - .5,yend = length(predictorlevels) - .5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Forest community level split -------------------------------------------------

# create frame for plotting, add rows of potentially possible connections:
d_semcoefs_f_plot_sep2 <- expand.grid(predictor = unique(d_semcoefs_f_plot_sep$predictor),
                                      response = unique(d_semcoefs_f_plot_sep$response)) %>% 
  mutate_all(~as.character(.)) %>% 
  left_join(d_semcoefs_f_plot_sep, by = c("predictor", "response")) %>% 
  mutate(exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude)) %>% 
  filter(!exclude)


var_order <- sapply(l_modlist_f_plot_sep_sel,
                    function(x) as.character(formula(x)[2]))

d_R2_f_plot_sep <- lapply(l_modlist_f_plot_sep_sel, 
                          function(x) data.frame(response = as.character(formula(x)[2]),
                                                 R2 = rsquared(x)$R.squared,
                                                 stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")))


d_semcoefs_f_plot_sep2 <- d_semcoefs_f_plot_sep2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_f_plot_sep) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "Idwcut", "Iharv", "Inonat",
                                                  var_order, "R2")),
         response = factor(response, levels = c("Idwcut", "Iharv", "Inonat", var_order)))


in_direct_paths_f_plot_sep <- pathways_f_plot_sep %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_f_plot_sep2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_f_plot_sep2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[3] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_f_plot_sep2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_f_plot_sep %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = 2.5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - .5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - .5,yend = length(predictorlevels) - .5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))


# Grassland community level combined -------------------------------------------

# create frame for plotting, add rows of potentially possible connections:
d_semcoefs_g_plot2 <- expand.grid(predictor = unique(d_semcoefs_g_plot$predictor),
                                  response = unique(d_semcoefs_g_plot$response)) %>% 
  mutate_all(~as.character(.)) %>% 
  left_join(d_semcoefs_g_plot, by = c("predictor", "response")) %>% 
  mutate(exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude)) %>% 
  filter(!exclude)


var_order <- sapply(l_modlist_g_plot_sel,
                    function(x) as.character(formula(x)[2]))

d_R2_g_plot <- lapply(l_modlist_g_plot_sel, 
                      function(x) data.frame(response = as.character(formula(x)[2]),
                                             R2 = rsquared(x)$R.squared,
                                             stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")))


d_semcoefs_g_plot2 <- d_semcoefs_g_plot2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_g_plot) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "LUI",
                                                  var_order, "R2")),
         response = factor(response, levels = c("LUI", var_order)))


in_direct_paths_g_plot <- pathways_g_plot %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_g_plot2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_g_plot2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[2] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5


textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_g_plot2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_g_plot %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = .5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - .5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - .5,yend = length(predictorlevels) - .5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Grassland community level split ----------------------------------------------

# create frame for plotting, add rows of potentially possible connections:
d_semcoefs_g_plot_sep2 <- expand.grid(predictor = unique(d_semcoefs_g_plot_sep$predictor),
                                      response = unique(d_semcoefs_g_plot_sep$response)) %>% 
  mutate_all(~as.character(.)) %>% 
  left_join(d_semcoefs_g_plot_sep, by = c("predictor", "response")) %>% 
  mutate(exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude)) %>% 
  filter(!exclude)


var_order <- sapply(l_modlist_g_plot_sep_sel,
                    function(x) as.character(formula(x)[2]))

d_R2_g_plot_sep <- lapply(l_modlist_g_plot_sep_sel, 
                          function(x) data.frame(response = as.character(formula(x)[2]),
                                                 R2 = rsquared(x)$R.squared,
                                                 stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")))


d_semcoefs_g_plot_sep2 <- d_semcoefs_g_plot_sep2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_g_plot_sep) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "F_std", "M_std", "G_std",
                                                  var_order, "R2")),
         response = factor(response, levels = c("F_std", "M_std", "G_std", var_order)))


in_direct_paths_g_plot_sep <- pathways_g_plot_sep %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_g_plot_sep2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_g_plot_sep2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[3] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))

d_semcoefs_g_plot_sep2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_g_plot_sep %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = 2.5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - .5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - .5,yend = length(predictorlevels) - .5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Forest species level combined ------------------------------------------------

d_semcoefs_f_sample2 <- expand.grid(predictor = unique(d_semcoefs_f_sample$predictor),
                                    response = unique(d_semcoefs_f_sample$response),
                                    stringsAsFactors = F) %>% 
  left_join(d_semcoefs_f_sample, by = c("predictor", "response")) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response),
         exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude),
         exclude = ifelse(grepl("i.fun.comp|i.tax.comp", predictor) &
                            response == "i.tax.comp.iszeroyes", TRUE, exclude)) %>% 
  filter(!exclude)

var_order <- sapply(l_modlist_f_sample_sel,
                    function(x) as.character(formula(x)[2]))
var_order[var_order == "i.tax.comp.iszero"] <- "i.tax.comp.iszeroyes"
var_order[var_order == "p.fgr.comp.grass"] <- "p.fgr.comp.grass1"
var_order[var_order == "p.fgr.comp.forb"] <- "p.fgr.comp.forb1"
var_order[var_order == "p.fgr.comp.tree"] <- "p.fgr.comp.tree1"
var_order[var_order == "p.fgr.comp.geophyt"] <- "p.fgr.comp.geophyt1"



d_R2_f_sample <- lapply(l_modlist_f_sample_sel, 
                        function(x) data.frame(response = as.character(formula(x)[2]),
                                               R2m = my_rsq(x)$Marginal,
                                               R2c = my_rsq(x)$Conditional,
                                               stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")),
         response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response))

d_semcoefs_f_sample2 <- d_semcoefs_f_sample2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_f_sample) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "ForMI",
                                                  var_order, "R2c", "R2m")),
         response = factor(response, levels = c( "ForMI", var_order)))

in_direct_paths_f_sample <- pathways_f_sample %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_f_sample2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_f_sample2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()

lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[2] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_f_sample2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_f_sample %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = .5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - 1.5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - 1.5,yend = length(predictorlevels) - 1.5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Forest species level split ---------------------------------------------------

d_semcoefs_f_sample_sep2 <- expand.grid(predictor = unique(d_semcoefs_f_sample_sep$predictor),
                                        response = unique(d_semcoefs_f_sample_sep$response),
                                        stringsAsFactors = F) %>% 
  left_join(d_semcoefs_f_sample_sep, by = c("predictor", "response")) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response),
         exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude),
         exclude = ifelse(grepl("i.fun.comp|i.tax.comp", predictor) &
                            response == "i.tax.comp.iszeroyes", TRUE, exclude)) %>% 
  filter(!exclude)

var_order <- sapply(l_modlist_f_sample_sep_sel,
                    function(x) as.character(formula(x)[2]))
var_order[var_order == "i.tax.comp.iszero"] <- "i.tax.comp.iszeroyes"
var_order[var_order == "p.fgr.comp.grass"] <- "p.fgr.comp.grass1"
var_order[var_order == "p.fgr.comp.forb"] <- "p.fgr.comp.forb1"
var_order[var_order == "p.fgr.comp.tree"] <- "p.fgr.comp.tree1"
var_order[var_order == "p.fgr.comp.geophyt"] <- "p.fgr.comp.geophyt1"



d_R2_f_sample_sep <- lapply(l_modlist_f_sample_sep_sel, 
                            function(x) data.frame(response = as.character(formula(x)[2]),
                                                   R2m = my_rsq(x)$Marginal,
                                                   R2c = my_rsq(x)$Conditional,
                                                   stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")),
         response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.fgr.comp.forb", "p.fgr.comp.forb1", response),
         response = ifelse(response == "p.fgr.comp.tree", "p.fgr.comp.tree1", response),
         response = ifelse(response == "p.fgr.comp.geophyt", "p.fgr.comp.geophyt1", response))

d_semcoefs_f_sample_sep2 <- d_semcoefs_f_sample_sep2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_f_sample_sep) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "Idwcut", "Iharv", "Inonat",
                                                  var_order, "R2c", "R2m")),
         response = factor(response, levels = c("Idwcut", "Iharv", "Inonat", var_order)))

in_direct_paths_f_sample_sep <- pathways_f_sample_sep %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_f_sample_sep2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_f_sample_sep2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[3] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_f_sample_sep2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_f_sample_sep %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = 2.5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - 1.5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - 1.5,yend = length(predictorlevels) - 1.5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Grassland species level combined ---------------------------------------------

d_semcoefs_g_sample2 <- expand.grid(predictor = unique(d_semcoefs_g_sample$predictor),
                                    response = unique(d_semcoefs_g_sample$response),
                                    stringsAsFactors = F) %>% 
  left_join(d_semcoefs_g_sample, by = c("predictor", "response")) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response),
         exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude),
         exclude = ifelse(grepl("i.fun.comp|i.tax.comp", predictor) &
                            response == "i.tax.comp.iszeroyes", TRUE, exclude)) %>% 
  filter(!exclude)

var_order <- sapply(l_modlist_g_sample_sel,
                    function(x) as.character(formula(x)[2]))
var_order[var_order == "i.tax.comp.iszero"] <- "i.tax.comp.iszeroyes"
var_order[var_order == "p.fgr.comp.grass"] <- "p.fgr.comp.grass1"
var_order[var_order == "p.tax.comp.legume"] <- "p.tax.comp.legume1"


d_R2_g_sample <- lapply(l_modlist_g_sample_sel, 
                        function(x) data.frame(response = as.character(formula(x)[2]),
                                               R2m = my_rsq(x)$Marginal,
                                               R2c = my_rsq(x)$Conditional,
                                               stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")),
         response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response))

d_semcoefs_g_sample2 <- d_semcoefs_g_sample2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_g_sample) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "LUI",
                                                  var_order, "R2c", "R2m")),
         response = factor(response, levels = c( "LUI", var_order)))

in_direct_paths_g_sample <- pathways_g_sample %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_g_sample2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_g_sample2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[2] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_g_sample2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_g_sample %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = .5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - 1.5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - 1.5,yend = length(predictorlevels) - 1.5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))

# Grassland species level split ------------------------------------------------

d_semcoefs_g_sample_sep2 <- expand.grid(predictor = unique(d_semcoefs_g_sample_sep$predictor),
                                        response = unique(d_semcoefs_g_sample_sep$response),
                                        stringsAsFactors = F) %>% 
  left_join(d_semcoefs_g_sample_sep, by = c("predictor", "response")) %>% 
  mutate(response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response),
         exclude = ifelse(is.na(estimate), TRUE, FALSE),
         exclude = ifelse(substr(predictor, 1, 1) == "p" &
                            substr(response, 1, 1) == "i", FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            (grepl("abund", response) | 
                               grepl("div", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("div", response), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("comp", predictor) &
                            grepl("comp", response) &
                            (grepl(".fun", response) |
                               (grepl(".fgr", response) & !grepl(".fun", predictor)) |
                               (grepl(".tax", response) & grepl(".tax", predictor))), 
                          FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("abund", predictor) &
                            grepl("abund", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(substr(predictor, 1, 1) == substr(response, 1, 1) &
                            grepl("div", predictor) &
                            grepl("div", response) &
                            (grepl(".tax.", predictor) |
                               grepl(".fun.", response)), FALSE, exclude),
         exclude = ifelse(predictor == response, TRUE, exclude),
         exclude = ifelse(grepl("i.fun.comp|i.tax.comp", predictor) &
                            response == "i.tax.comp.iszeroyes", TRUE, exclude)) %>% 
  filter(!exclude)

var_order <- sapply(l_modlist_g_sample_sep_sel,
                    function(x) as.character(formula(x)[2]))
var_order[var_order == "i.tax.comp.iszero"] <- "i.tax.comp.iszeroyes"
var_order[var_order == "p.fgr.comp.grass"] <- "p.fgr.comp.grass1"
var_order[var_order == "p.tax.comp.legume"] <- "p.tax.comp.legume1"



d_R2_g_sample_sep <- lapply(l_modlist_g_sample_sep_sel, 
                            function(x) data.frame(response = as.character(formula(x)[2]),
                                                   R2m = my_rsq(x)$Marginal,
                                                   R2c = my_rsq(x)$Conditional,
                                                   stringsAsFactors = F)) %>% 
  do.call(rbind, .) %>% 
  gather(predictor, sign.code, -response) %>% 
  mutate(sign.code = as.character(formatC(sign.code, digits = 2, format = "f")),
         response = ifelse(response == "i.tax.comp.iszero", "i.tax.comp.iszeroyes", response),
         response = ifelse(response == "p.fgr.comp.grass", "p.fgr.comp.grass1", response),
         response = ifelse(response == "p.tax.comp.legume", "p.tax.comp.legume1", response))

d_semcoefs_g_sample_sep2 <- d_semcoefs_g_sample_sep2 %>% 
  mutate(sign.code = as.character(sign.code),
         response = as.character(response)) %>% 
  bind_rows(d_R2_g_sample_sep) %>% 
  mutate(predictor = factor(predictor, levels = c("regionHAI", "regionSCH",
                                                  "F_std", "M_std", "G_std",
                                                  var_order, "R2c", "R2m")),
         response = factor(response, levels = c( "F_std", "M_std", "G_std", var_order)))

in_direct_paths_g_sample_sep <- pathways_g_sample_sep %>% 
  mutate(response = ifelse(substr(path, regexpr(" > ", path) + 3, regexpr(" > ", path) + 14) == "a.herb2.prop",
                           "direct", "indirect")) %>% 
  group_by(predictor, response) %>% 
  summarise(estimate = sum(estimate)) %>% 
  ungroup() %>% 
  spread(response, estimate) %>% 
  rowwise() %>% 
  mutate(net_effect = sum(c(direct, indirect), na.rm = T)) %>% 
  gather(response, estimate, -predictor)

responselevels <- c(levels(droplevels(d_semcoefs_g_sample_sep2$response)),
                    "direct", "indirect", "net_effect")
predictorlevels <- levels(droplevels(d_semcoefs_g_sample_sep2$predictor))


lines_h <- data.frame(predictorlevels) %>% 
  mutate(predictorlevels = as.character(predictorlevels),
         nextlevel = c(predictorlevels[-1], NA)) %>% 
  rowwise() %>% 
  mutate(line = ifelse(grepl("region", predictorlevels) &
                         !grepl("region", nextlevel), 1, NA),
         line = ifelse(substr(predictorlevels, 1, 2) != "p." &
                         substr(nextlevel, 1, 2) == "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "p." &
                         substr(nextlevel, 1, 2) != "p.",
                       1, line),
         line = ifelse(substr(predictorlevels, 1, 2) == "i." &
                         substr(nextlevel, 1, 2) != "i.",
                       1, line),
         line = ifelse(grepl("comp", predictorlevels) &
                         !grepl("comp", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("abund", predictorlevels) &
                         !grepl("abund", nextlevel),
                       min(line, 2, na.rm = T), line),
         line = ifelse(grepl("tax", predictorlevels) &
                         !grepl("tax", nextlevel),
                       min(line, 3, na.rm = T), line)) %>% 
  ungroup()


lines_v <- lines_h[-c(1:3), ]

h_lines_1 <- which(lines_h$line == 1) + .5
h_lines_2 <- which(lines_h$line == 2) + .5
h_lines_3 <- which(lines_h$line == 3) + .5
v_lines_1 <-  which(lines_v$line == 1) + .5
v_lines_1 <- c(v_lines_1, v_lines_1[3] + 1)
v_lines_2 <-  which(lines_v$line == 2) + .5
v_lines_3 <-  which(lines_v$line == 3) + .5

textcol_x <- ifelse(substr(responselevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(responselevels, 1, 2) == "i.", "#ED7D31", 1))
textcol_y <- ifelse(substr(predictorlevels, 1, 2) == "p.", "#70AD47", 
                    ifelse(substr(predictorlevels, 1, 2) == "i.", "#ED7D31", "black"))


d_semcoefs_g_sample_sep2 %>% 
  mutate(response = as.character(response),
         predictor = as.character(predictor)) %>% 
  bind_rows(in_direct_paths_g_sample_sep %>% 
              mutate(sign.code = formatC(estimate, digits = 2, format = "f"),
                     predictor = as.character(predictor))) %>% 
  mutate(sign.code = ifelse(is.na(sign.code), "", sign.code),
         response = factor(response, levels = responselevels),
         predictor = factor(predictor, levels = predictorlevels)) %>% 
  ggplot(aes(x = response, y = predictor, fill = estimate, label = sign.code)) +
  geom_tile() +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "grey25") +
  geom_tile() +
  annotate(geom = "rect", xmin = 2.5, xmax = length(responselevels) - 2.5, 
           ymin = length(predictorlevels) - 1.5, ymax = Inf, 
           fill = "white") +
  geom_text() +
  geom_vline(xintercept = v_lines_1, size = 2) +
  geom_hline(yintercept = h_lines_1, size = 2) +
  geom_vline(xintercept = v_lines_2, lty = 2) +
  geom_hline(yintercept = h_lines_2, lty = 2) +
  geom_vline(xintercept = v_lines_3, lty = 3) +
  geom_hline(yintercept = h_lines_3, lty = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_fill_gradient2(na.value = "grey80", name = "Estimate") +
  scale_x_discrete(labels = varlabs) +
  scale_y_discrete(labels = varlabs) +
  xlab("Response") +
  ylab("Predictor") +
  annotate("segment", x = 1, xend = -10, 
           y = length(predictorlevels) - 1.5,yend = length(predictorlevels) - 1.5, 
           size = 2) +
  annotate("segment", y = 1, yend = -10, 
           x = length(responselevels) - 2.5, xend = length(responselevels) - 2.5, 
           size = 2) +
  coord_cartesian(xlim = c(1, length(responselevels)),
                  ylim = c(1, length(predictorlevels)),
                  clip="off") +
  theme(legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_text(color = textcol_x),
        axis.text.y = element_text(color = textcol_y))



# Spatial autocorrelation ######################################################
################################################################################.

# Check residuals (variograms) -------------------------------------------------

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.comm.formi, data = d_comm) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.comm.formi.sep, data = d_comm) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.comm.lui, data = d_comm) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.comm.lui.sep, data = d_comm) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.sample.formi, data = d_spec) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.sample.formi.sep, data = d_spec) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.sample.lui, data = d_spec) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

lapply(c("Alb", "Hai", "Sch"), f_vario, 
       mod = mod.sample.lui.sep, data = d_spec) %>% 
  do.call(rbind, .) %>% 
  mutate(dist = dist / 1000) %>% 
  ggplot(aes(x = dist, y = gamma, col = region)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ region, scales = "free_x") +
  theme(panel.spacing = unit(20, "mm"),
        legend.position = "none") +
  xlab("Distance [km]")

# Check random effect (bubble plots) -------------------------------------------

f_ranefplots(mod.sample.formi)
f_ranefplots(mod.sample.formi.sep)
f_ranefplots(mod.sample.lui)
f_ranefplots(mod.sample.lui.sep)

# additional figures ###########################################################
################################################################################.

col.plgroups <- c(tree = "#1b9e77",
                  forb = "#d95f02",
                  geophyt = "#7570b3",
                  grass = "#e7298a",
                  fern = "#66a61e",
                  legume = "red")

labs.plgroups <- c(tree = "Trees",
                   forb = "Forbs",
                   geophyt = "Geophytes",
                   grass = "Grasses",
                   fern = "Ferns",
                   legume = "Legumes")

# shares of plant groups in forests (along Inonat gradient) --------------------

p1 <- d_spec %>% 
  filter(!is.na(Inonat)) %>% 
  group_by(plot, Inonat, pl.group) %>% 
  summarise(coversum = sum(cover)) %>% 
  group_by(plot, Inonat) %>% 
  mutate(coverprop = coversum / sum(coversum)) %>% 
  ungroup() %>% 
  select(-coversum) %>% 
  spread(pl.group, coverprop, fill = 0) %>% 
  gather(pl.group, coverprop, -c(plot, Inonat)) %>% 
  mutate(pl.group = factor(pl.group, levels = names(labs.plgroups))) %>% 
  ggplot(aes(x = Inonat, y = coverprop, col = pl.group, fill = pl.group)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_color_manual(values = col.plgroups, labels = labs.plgroups, name = "Plant group") +
  scale_fill_manual(values = col.plgroups, labels = labs.plgroups, name = "Plant group") +
  scale_y_continuous(labels = paste0(seq(0, 100, 25), "%"), name = "Percentage of sampled community") +
  theme(legend.position = "none")

p2 <-
  d_spec %>% 
  filter(!is.na(Inonat)) %>% 
  group_by(plot, Inonat, pl.group) %>% 
  summarise(a.herb2.prop = sum(a.herb2.prop * prop.cover / sum(prop.cover))) %>% 
  mutate(pl.group = factor(pl.group, levels = names(labs.plgroups))) %>% 
  ggplot(aes(x = Inonat, y = a.herb2.prop, col = pl.group, fill = pl.group)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_color_manual(values = col.plgroups, labels = labs.plgroups, name = "Plant group") +
  scale_fill_manual(values = col.plgroups, labels = labs.plgroups, name = "Plant group") +
  scale_y_log10(name = "Mean herbivory rate", breaks = c(0.001, 0.010, 0.1), labels = paste0(c(0.1, 1, 10), "%"))

plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.2), labels = c("a", "b"))

# shares of damage types and plant groups overall ------------------------------

meanshares <-
  d_spec %>% 
  mutate(pl.group = ifelse(!pl.group %in% c("grass", "legume") & system == "grassland", "forb", pl.group)) %>% 
  group_by(system, plot, pl.group) %>% 
  summarise(coversum = sum(cover)) %>% 
  group_by(system, plot) %>% 
  mutate(coverprop = coversum / sum(coversum)) %>% 
  select(-coversum) %>%
  spread(pl.group, coverprop, fill = 0) %>%
  gather(pl.group, coverprop, -c(system, plot)) %>%
  group_by(system, pl.group) %>% 
  summarise(meanshare = mean(coverprop)) %>% 
  filter(meanshare != 0) %>% 
  mutate(pl.group = factor(pl.group, levels = names(labs.plgroups))) %>% 
  arrange(system, pl.group) %>% 
  group_by(system) %>% 
  mutate(meanshare_cumsum = cumsum(meanshare),
         meanshare_cumsum = meanshare_cumsum - meanshare / 2)

d_colplot <- d_spec %>% 
  select(system, pl.group, a.chewing.prop, a.scraping.prop, a.sucking.prop, 
         a.mines.prop, a.gallstot.prop) %>% 
  mutate(pl.group = ifelse(!pl.group %in% c("grass", "legume") & system == "grassland", "forb", pl.group)) %>% 
  group_by(system, pl.group) %>% 
  summarise_all(~mean(.)) %>% 
  ungroup() %>% 
  gather(damage.type, a.prop, -c(system, pl.group)) %>% 
  group_by(system, pl.group) %>% 
  mutate(a.prop = a.prop / sum(a.prop)) %>% 
  ungroup() %>% 
  left_join(meanshares, by = c("system", "pl.group")) %>%
  mutate(damage.type = gsub("a\\.", "", damage.type),
         damage.type = gsub("\\.prop", "", damage.type),
         damage.type = factor(damage.type, 
                              levels = c("chewing", "scraping", "sucking", "mines", "gallstot"),
                              labels = c("Chewing", "Scraping", "Sucking", "Mines", "Galls")),
         pl.group = factor(pl.group, levels = names(labs.plgroups),
                           labels = labs.plgroups),
         system = factor(system, labels = c("Forests", "Grasslands"))) 

d_forest <- d_colplot %>% 
  filter(system == "Forests") %>% 
  ggplot(aes(x = meanshare_cumsum, y = a.prop, fill = damage.type)) +
  geom_col(aes(width = meanshare), col = 1) +
  scale_x_continuous(breaks = meanshares$meanshare_cumsum[meanshares$system == "forest"],
                     labels = labs.plgroups[labs.plgroups != "Legumes"]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = paste0(seq(0, 100, 25), "%"), name = "Mean percentage of total damage") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"), name = "Damage type")

d_grassland <- d_colplot %>% 
  filter(system == "Grasslands") %>% 
  ggplot(aes(x = meanshare_cumsum, y = a.prop, fill = damage.type)) +
  geom_col(aes(width = meanshare), col = 1) +
  scale_x_continuous(breaks = meanshares$meanshare_cumsum[meanshares$system == "grassland"],
                     labels = c("Forbs", "Grasses", "Legumes")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(labels = paste0(seq(0, 100, 25), "%"), name = "Mean percentage of total damage") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"), name = "Damage type")

plot_grid(plot_grid(ggplot() + theme_nothing(), ggplot() + theme_nothing(), labels = c("a", "b")),
          plot_grid(d_forest, d_grassland, align = "h"),
          rel_heights = c(.2, 1), nrow = 2)

