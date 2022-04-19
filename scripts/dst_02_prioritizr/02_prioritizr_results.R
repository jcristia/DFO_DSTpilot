
# build official results scenarios

library(prioritizr)
library(tidyverse)
library(scales)
library(glue)


########################################
# paths and data
########################################

root <- r'(C:\Users\jcristia\Documents\GIS\DFO\DST_pilot)'
gdb <- file.path(root, 'spatial/03_working/dst_grid.gdb')
fc <- 'dst_grid1km_PUVSP'
out_path <- file.path(root, 'scripts/dst_02_prioritizr/outputs')
targets_csv <- read_csv(file.path(root, 'scripts/dst_01_preprocess/DSTpilot_spatialData - naming_scheme.csv'))
gdb_mpa <- file.path(root, 'spatial/03_working/dst_mpas.gdb')
mpa_pu_fc <- 'mpas_rcas_marxan'
zoning <- read_csv(file.path(root, 'scripts/dst_02_prioritizr/Zoning frameworks - zones_features.csv'))


########################################
# pus, features, targets, mpas locked-in
########################################

# read in feature class containing Planning Units and area of each feature in each PU
pus <- st_read(gdb, layer = fc)

# features list
features <- names(pus)[grepl('eco', names(pus)) | grepl('hu', names(pus)) ]

# pus with MPAs to lock in
mpas <- st_read(gdb_mpa, layer=mpa_pu_fc)
mpas <- st_drop_geometry(mpas)
threshold <- 150000 # set a threshold of overlap so that a sliver doesn't get locked in
mpas <- mpas[mpas$Shape_Area > threshold, ]
mpas <- select(mpas, uID)
pus <- left_join(
  pus, 
  mpas, 
  by='uID', 
  keep=TRUE)
pus <- mutate(pus, lockedin = case_when(uID.y > 0 ~ TRUE, is.na(uID.y) ~ FALSE))
pus <- pus %>% select(-uID.y) %>% rename(uID=uID.x)

# create dataframe of locked constraints for zones (different format)
# must be numeric value in 'status' field
# MPAs are all in zone 1
locked_df <- pus %>% 
  select(c('uID', 'lockedin')) %>% 
  st_drop_geometry() %>%
  rename(pu=uID) %>%
  filter(lockedin==TRUE)
locked_df$status <- 1
locked_df <- select(locked_df, -lockedin)
locked_df$zone <- 'zone1'

# targets
# need to be in same order as features
targets_csv <- rename(
  targets_csv,
  target_low = `Coarse-filter feature; Low Target Range (10%)`,
  target_med = `Coarse-filter feature; Medium Target Range (20%)`,
  target_hig = `Coarse-filter feature; High Target Range (30%)`)
targets_csv <- select(targets_csv, c(processed_file, target_low, target_med, target_hig))
targets <- left_join(
  as.data.frame(features), 
  targets_csv, 
  by=c('features'='processed_file'), 
  keep=TRUE)
# percent to proportion
targets <- mutate_at(targets, c('target_low', 'target_med', 'target_hig'), function(x)(x/100))





########################################
# function to evaluate each solution
########################################

evaluate_solution <- function(p1, s1) {
  
  # get solution columns for the first solution
  # this accounts for multiple zones
  sols <- grep('solution_1', colnames(s1), value=TRUE)
  
  # evaluate the performance of the first solution
  eval_numsel <- eval_n_summary(p1, s1[,sols])
  eval_cost <- eval_cost_summary(p1, s1[,sols])
  eval_targets <- eval_target_coverage_summary(p1, s1[,sols])
  eval_boundary <- eval_boundary_summary(p1, s1[,sols])
  eval_feat <- eval_feature_representation_summary(p1, s1[,sols])
  
  # set up table
  s2 <- tibble::rowid_to_column(s1, 'id_index')
  
  # Format for output
  # If a solution has zones then it the zones column is in list type since a
  # target can correspond to multiple zones
  # Place zone id in new column
  if("zone" %in% colnames(eval_targets))
  {
    eval_targets$zoning <- vapply(eval_targets$zone, FUN.VALUE = character(1), paste, sep = " & ")
    eval_targets <- select(eval_targets, -zone)
  }
  
  # output solutions, target evaluation, and summary evaluation to csv
  s2 <- as_tibble(s2)
  s2 <- select(s2, -'Shape')
  write.csv(s2, 'solution.csv', row.names=FALSE)
  write.csv(eval_targets, 'eval_targets.csv', row.names=FALSE)
  write.csv(eval_boundary, 'eval_boundary.csv', row.names=FALSE)
  write.csv(eval_feat, 'eval_featurerep.csv', row.names=FALSE)
  df_summ <- tibble(zone=eval_numsel$summary, total_selected=eval_numsel$cost, cost=eval_cost$cost)
  write.csv(df_summ, 'eval_summary.csv', row.names=FALSE)
  
}



##############################################################
# solution 101:
# 
# no boundary penalty
# min_shortfall_objective to identify the ones that fall short
# However, it will spend budget to exceed targets, so it will
# look different than min_set_objective and we can't compare
# it to Marxan scenarios.
# See Jeff’s Jan 25 answer here: 
# https://github.com/prioritizr/prioritizr/issues/224
# The thing on the shopping list and his point 2.
##############################################################

out_folder <- 's101_minshortfall'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# set the cost to zero for pus already covered by MPAs
pus$COST[pus$lockedin == TRUE] <- 0

# zone targets matrix
zone_names = c('zone1', 'zone2', 'zone3', 'zone4')
tz <- matrix(NA, nrow=length(features), ncol=4, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_4zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# Create zones object
z <- zones(
  features, features, features, features,
  zone_names=zone_names,
  feature_names=features
)

# total allowable cost for min_shortfall_objective
total_cost <- sum(pus$COST) + 1000

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1)

sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)



##############################################################
# solution 102:
# 
# no boundary penalty
# min_set_objective
# adjust targets of ones that were not met in previous solution
##############################################################

out_folder <- 's102_minset_adjtarg'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# set the cost to zero for pus already covered by MPAs
pus$COST[pus$lockedin == TRUE] <- 0

# zone targets matrix
zone_names = c('zone1', 'zone2', 'zone3', 'zone4')
tz <- matrix(NA, nrow=length(features), ncol=4, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_4zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# adjust targets manually
tz2 <- tz
tz2['hu_ot_dredgingsites',2] <- 0.0
tz2['hu_tr_vesseltraffic_d',3] <- 0.28
tz2['hu_co_fishing_dive_d',4] <- 0.49
tz2['hu_ot_underwaterinfrastructure',3]<-0.62
tz2['hu_rf_fishing_groundfish',4]<-0.62
tz2['hu_tr_portsandterminals',2]<-0.64
tz2['hu_rf_fishing_crab',4]<-0.64
tz2['hu_ot_citypopulation_d',2]<-0.65
tz2['hu_rf_fishing_prawnandshrimp',4]<-0.68
tz2['hu_ot_floatingstructures',2]<-0.69

# Create zones object
z <- zones(
  features, features, features, features,
  zone_names=zone_names,
  feature_names=features
)

# read in previous solution as start solution
pr <- read_rds(file.path(root, 'scripts/dst_02_prioritizr/outputs/s12_minsetobj_nobound/s12_minsetobj_nobound.rds'))
sols <- grep('solution_1', colnames(pr), value=TRUE)

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz2) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, start_solution = pr[,sols], time_limit=14400)

sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)



##############################################################
# solution 103:
# 
# same as above but with boundary penalty
##############################################################

out_folder <- 's103_minset_adjtarg_bp'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# set the cost to zero for pus already covered by MPAs
pus$COST[pus$lockedin == TRUE] <- 0

# zone targets matrix
zone_names = c('zone1', 'zone2', 'zone3', 'zone4')
tz <- matrix(NA, nrow=length(features), ncol=4, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_4zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# adjust targets manually
tz2 <- tz
tz2['hu_ot_dredgingsites',2] <- 0.0
tz2['hu_tr_vesseltraffic_d',3] <- 0.28
tz2['hu_co_fishing_dive_d',4] <- 0.49
tz2['hu_ot_underwaterinfrastructure',3]<-0.62
tz2['hu_rf_fishing_groundfish',4]<-0.62
tz2['hu_tr_portsandterminals',2]<-0.64
tz2['hu_rf_fishing_crab',4]<-0.64
tz2['hu_ot_citypopulation_d',2]<-0.65
tz2['hu_rf_fishing_prawnandshrimp',4]<-0.68
tz2['hu_ot_floatingstructures',2]<-0.69

# Create zones object
z <- zones(
  features, features, features, features,
  zone_names=zone_names,
  feature_names=features
)

# Get boundary length data for all pus
pus_bd <- boundary_matrix(pus)

# Re-scale boundary length data to match order of magnitude of costs to avoid
# numerical issues and reduce run times.
pus_bd@x <- rescale(pus_bd@x, to=c(10,30))

# create zone matrix which favors clumping planning units that are
# allocated to the same zone together - this is the default
zb <- diag(4)

# read in previous solution as start solution
pr <- read_rds(file.path(root, 'scripts/dst_02_prioritizr/outputs/s102_minset_adjtarg/s102_minset_adjtarg.rds'))
sols <- grep('solution_1', colnames(pr), value=TRUE)

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz2) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=0.000001, data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, start_solution = pr[,sols], time_limit=28800)

sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)




##############################################################
# solution 104:
# 
# same as 101, but without MPAs locked in
# We want to see how it changes meetings targets
# no boundary penalty
# min_shortfall_objective to identify the ones that fall short
##############################################################

out_folder <- 's104_minshortfall_nolockedin'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# set the cost to zero for pus already covered by MPAs
#pus$COST[pus$lockedin == TRUE] <- 0

# zone targets matrix
zone_names = c('zone1', 'zone2', 'zone3', 'zone4')
tz <- matrix(NA, nrow=length(features), ncol=4, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_4zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# Create zones object
z <- zones(
  features, features, features, features,
  zone_names=zone_names,
  feature_names=features
)

# total allowable cost for min_shortfall_objective
total_cost <- sum(pus$COST) + 1000

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=28800)

sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)



