
# explore and build up scenarios

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
#pus <- mutate(pus, lockedin = ifelse(uID.y > 0, TRUE, FALSE))
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
# solution 1:
# 
# Just ECO features
# no zones
# no BLM 
# low targets
# no locked in
##############################################################

out_folder <- 's01_base'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# remove hus
features_eco <- features[startsWith(features, 'eco_')]
targets_eco <- targets[startsWith(targets$features, 'eco'), ]

p <- problem(pus, features_eco, cost_column='COST') %>%
  add_min_set_objective() %>%
  add_relative_targets(targets_eco$target_low) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1)
#add_gap_portfolio(number_solutions = 1, pool_gap = 0)
#add_cbc_solver(gap=0) %>% # says this is fastest open source solver
#add_lpsymphony_solver(gap=0) %>% # however, this performs well and looks cleaner
#add_cuts_portfolio(number_solutions = 2) # when not using gurobi
s <- solve(p)
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)



##############################################################
# solution 2:
# 
# lockedin mpas
#
# For lockedin scenarios, I am not changing the cost of MPA
# pus to zero. I want to see how much more efficient a solution
# is when we don't have to lock in mpas.

# However, I do do this for when I add in a boundary penalty.
##############################################################

out_folder <- 's02_lockedin'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# remove hus
features_eco <- features[startsWith(features, 'eco_')]
targets_eco <- targets[startsWith(targets$features, 'eco'), ]

p <- problem(pus, features_eco, cost_column='COST') %>%
  add_min_set_objective() %>%
  add_relative_targets(targets_eco$target_low) %>%
  add_locked_in_constraints('lockedin') %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1)
s <- solve(p)
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)



##############################################################
# solution 3:
# 
# 2 zone problem
##############################################################

out_folder <- 's03_zones2'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# zone targets matrix
zone_names = c('zone1', 'zone2')
tz <- matrix(NA, nrow=length(features), ncol=2, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_2zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# Create a zones object
# You need to duplicate the feature data for the number of zones that you have.
# This is a bit repetitive for us. It's set up this way in case your planning
# units cover different extents per zone.
z <- zones(
  features, features,
  zone_names=zone_names,
  feature_names=features
)

# The cost column needs to be duplicated. This is for if a planning unit has
# a different cost of protecting based on which zone it gets allocated to.
p <- problem(pus, z, cost_column=c('COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=1800)
s <- solve(p)
evaluate_solution(p, s)

# irreplaceability
# gap_replace <- 95 # set super high for testing
# sols <- grep('solution_1', colnames(s), value=TRUE) # get all zones of 1st solution
# rc1 <- p %>%
#   add_gurobi_solver(gap=gap_replace, verbose=TRUE, threads = parallel::detectCores(TRUE)-1) %>%
#   eval_replacement_importance(s[,sols])
# saveRDS(rc1, glue('{out_folder}_replaceImp.rds'))
# rc1 <- tibble::rowid_to_column(rc1, 'id_index') # add id from index for join
# rc1 <- as.data.frame(rc1)
# rc1 <- select(rc1, -geometry) # drop shape field
# s2 <- tibble::rowid_to_column(s, 'id_index')
# s2 <- left_join(s2, rc1, by='id_index')
# s2 <- as_tibble(s2)
# s2 <- select(s2, -'Shape')
# write.csv(s2, 'irreplaceability.csv', row.names=FALSE)




##############################################################
# solution 4:
# 
# 2 zone problem
# manual locked in constraints for 1 zone
##############################################################

out_folder <- 's04_zones2_lockedin'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# zone targets matrix
zone_names = c('zone1', 'zone2')
tz <- matrix(NA, nrow=length(features), ncol=2, dimnames=list(features,zone_names))
for(feat in features){
  t <- targets[targets$features==feat,'target_low']
  z <- zoning[zoning$feature==feat, 'zone_20220315_2zones'][[1]]
  tz[feat, z] <- t
}
tz <- replace_na(tz, 0.0) # features are only targeted in one zone

# Create zones object
z <- zones(
  features, features,
  zone_names=zone_names,
  feature_names=features
)

p <- problem(pus, z, cost_column=c('COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=1800)
#s <- solve(p)
#evaluate_solution(p, s)

# can't meet targets with hus set at 70% and MPAs locked in to zone 1.
# You get an error during solve, and I'm not sure why it doesn't at least run 
# and tell me that targets were not met in the targets spreadsheet.
# Well actually, this is a Gurobi thing. It makes sense.

# From the documentation, it looks like for add_min_set_objective() that:
# "targets for features will always be met", which is a weird wording. I'm
# understanding that to mean that the solver will fail it if can't meet them,
# instead of outputting the best solution.

# I investigated the hu data overlap with MPAs and dredging is the only one that
# completely overlaps. Try setting that at 0.
#tz2 <- tz
#tz2['hu_ot_dredgingsites',2] <- 0.0
#p2 <- p %>% add_relative_targets(tz2)
#s2 <- solve(p2, force=TRUE)
# still no

# Keep decreasing zone 2 targets until I can get it to run.
#tz2[,2] <- ifelse(tz2[,2]>0, 0.52, tz2[,2])
#p3 <- p %>% add_relative_targets(tz2)
#s3 <- solve(p3, force=TRUE)
# ok, so it will run with dredging set to 0 and all other hu targets at a max of 52%

# Now try the original problem but with add_min_shortfall_objective:
total_cost <- sum(pus$COST)
p4 <- p %>% add_min_shortfall_objective(total_cost)
s4 <- solve(p4, force=TRUE)
evaluate_solution(p4, s4)
# Oh wow, this is super useful!!!
# First, there's likely a few numerical errors, since a few eco features show
# up as not being met, but I think we can ignore these. Otherwise, it looks like
# it is only 3 hu features that are an issue: dredging, ports and terminals, and
# dive fishing.
# Ports and terminals gets super close, so we shouldn't worry about that.
# Dive fishing also gets fairly close.





##############################################################
# solution 5:
# 
# 4 zone problem
# no constraints
##############################################################

out_folder <- 's05_zones4'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

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

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=1800)
#s <- solve(p)
# nope doesn't work
# try with min shortfall objective

total_cost <- sum(pus$COST)
p2 <- p %>% add_min_shortfall_objective(total_cost)
s2 <- solve(p2, force=TRUE)
evaluate_solution(p2, s2)
# its just vessel traffic and recreation groundfish fishing that can't be met






##############################################################
# solution 6:
# 
# 4 zone problem
# locked in mpa constraints
##############################################################

out_folder <- 's06_zones4_constraints'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

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

# from the previous problem, I already know that I can't meet it with the normal
# objective, so switch to min_shortfall
total_cost <- sum(pus$COST)
p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=1800)
s <- solve(p, force=TRUE)
evaluate_solution(p, s)

# many more not meeting targets, but most are not significant








##############################################################
# BOUNDARY STUFF
##############################################################


##############################################################
# solution 7 and 8
# 
# Add boundary penalty and find "sweet spot" starting point
# 
# Boundary penalty, general:
# High penalty value favors more clumping.
# Edge_factor of 0.5 means that edges receive half the penalty,
# this is so that you do not over penalize pus that do not have
# neighbors.
# Behind the scenes there is a boundary matrix that gets 
# created that calculates the amount of shared boundary of the
# solution.
# Penalty: essentially, the length of the boundary gets added 
# to the cost.Therefore, the units of the boundary should be 
# somewhat in relation to the cost. If there is a penalty of 
# 1 then the boundaries as is get added in.
#
# Since boundary penalty is added the cost, a penalty too small results in
# having no effect and a boundary penalty too big results in just one large
# reserve.
#
# Adding a boundary penalty really slows down processing time.
# So it's not ideal to start blind and just test a ton of values
# See this vignette and notes below for calibrating:
# https://prioritizr.net/articles/calibrating_trade-offs_tutorial.html#preliminary-processing
# 
##############################################################

out_folder <- 's07_boundary_cohon1'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

# set the cost to zero for pus already covered by MPAs
# These are already protected and there shouldn't be a cost to select them.
# ALSO, this is important for properly configuring the boundary penalty.
# "This is so that the total cost estimates of the prioritization reflects the 
# total cost of establishing new protected areas – not just total land value. In
# other words, we want the total cost estimate for a prioritization to reflect 
# the cost of implementing conservation actions."
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

# Get boundary length data for all pus
pus_bd <- boundary_matrix(pus)
# This is a matrix. Off-diagonal is the length of boundary shared with another cell.
# In my case, this will be None or 1000.
# The diagonal will be the length of boundary that doesn't have neighbors.
# This will be None,1000,2000, or 3000.

# Re-scale boundary length data to match order of magnitude of costs to avoid
# numerical issues and reduce run times.
# My costs are all based on an area of 100. This was scaled down from 1,000,000
# with edges of 1000. Therefore, consider my edges to be 10.
pus_bd@x <- rescale(pus_bd@x, to=c(10,30))

# create zone matrix which favors clumping planning units that are
# allocated to the same zone together - note that this is the default
zb <- diag(4)


########## Determining boundary penalty sweet spot

# FROM p87 in the Marxan Best Practices manual:
# See document or Evernote for plot.
# "A third way to pick values for BLM is to use a weighting method":
# "First, set BLM to 0 and optimise cost to find the lowest cost solution 
# possible. Calculate the cost and boundary length for that solution and plot it
# as point X, the minimum cost solution (Figure 8.4). Then set all costs to 
# zero and BLM to 1 to find point Y (the minimum possible boundary solution). 
# Calculate the slope of line “a” connecting those two solutions in objective 
# space: (Cost(X) – Cost(Y))/(Boundary(X)‐Boundary(Y)). This is the estimated
# trade‐off curve with two points. Use the absolute value of the slope of line
# “a” as the BLM and reset all costs back to their original values. In this 
# example, the BLM would be 0.00088525. This value represents a “sweet spot” on 
# the trade‐off curve between minimizing cost and minimizing boundary length.
# Small changes in BLM around this “sweet spot” value are likely to make the 
# largest changes in spatial patterns of selected reserve networks. BLM values 
# of 10, 1, and even 0.1 are all so much higher than this “sweet spot” on the 
# curve that they all yield similar reserve networks—clustering reserves to the 
# maximum possible extent." 
# "Run Marxan again to locate point Z. With three solutions, the trade‐off curve
# is estimated as dashed lines “b” and “c.” If the solutions represented by 
# point Z are more clustered than desired, the process can be repeated with line
# “c” in order to find a new value for BLM."

# Then I follow the code from the Prioritizr vignette which does the above and
# refers to Cohon et al. 1979 as establishing this method.


# generate ideal prioritization based on cost criteria
total_cost <- sum(pus$COST) + 1000
p_cost <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
s_cost <- solve(p_cost, force=TRUE)
saveRDS(s_cost, glue('{out_folder}_s_cost.rds'))

# Generate ideal prioritization based on spatial fragmentation criteria
# Set costs to zero so that they aren't considered and boundary to 1 so that it
# is the only thing considered.

# create cost column of zeros
pus$COSTZERO <- 0
total_cost <- sum(pus$COST) + 1000
p_costzero <- problem(pus, z, cost_column=c('COSTZERO', 'COSTZERO', 'COSTZERO', 'COSTZERO')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=1, data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
s_costzero <- solve(p_costzero, force=TRUE)
saveRDS(s_costzero, glue('{out_folder}_s_costzero.rds'))

# generate problem formulation with costs and boundary penalties for
# calculating performance metrics
total_cost <- sum(pus$COST) + 1000
p_base <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_boundary_penalties(penalty=1, data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions()

# calculate performance metrics for ideal cost prioritization
s_cost_metrics <- tibble(
  total_cost1 = eval_cost_summary(p_base, s_cost[, c("solution_1_zone1", "solution_1_zone2", "solution_1_zone3", "solution_1_zone4")])$cost,
  total_boundary_length =
    eval_boundary_summary(p_base, s_cost[, c("solution_1_zone1", "solution_1_zone2", "solution_1_zone3", "solution_1_zone4")])$boundary)
# calculate performance metrics for ideal boundary length prioritization
s_costzero_metrics <- tibble(
  total_cost1 = eval_cost_summary(p_base, s_costzero[, c("solution_1_zone1", "solution_1_zone2", "solution_1_zone3", "solution_1_zone4")])$cost,
  total_boundary_length =
    eval_boundary_summary(p_base, s_costzero[, c("solution_1_zone1", "solution_1_zone2", "solution_1_zone3", "solution_1_zone4")])$boundary)

# calculate penalty value based on Cohon et al. 1979
cohon_penalty <- abs(
  (s_cost_metrics$total_cost1 - s_costzero_metrics$total_cost1) /
    (s_cost_metrics$total_boundary_length - s_costzero_metrics$total_boundary_length))
# round to 5 decimal places to avoid numerical issues during optimization
cohon_penalty <- round(cohon_penalty, 5)
saveRDS(cohon_penalty, 'cohon_penalty.rds')
cohon_penalty <- read_rds('cohon_penalty.rds')


# generate prioritization using penalty value calculated using Cohon et al. 1979
# I get 5 sweet spot values: 1 for overall and 4 for zones.
# Try at least the first 2 to get a sense of how things look.

# Overall penalty value:
total_cost <- sum(pus$COST) + 1000
p_cohon1 <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=cohon_penalty[1], data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
s_cohon1 <- solve(p_cohon1, force=TRUE)
saveRDS(s_cohon1, glue('{out_folder}_sweetspot.rds'))
evaluate_solution(p_cohon1, s_cohon1)

##################
# so now lower the boundary until I see some changes
# it doesn't have to be perfect. Just judge the sensitivity.
out_folder <- 's08_boundary_cohon2'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)
total_cost <- sum(pus$COST) + 1000
p_cohon2 <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=(cohon_penalty[1]/100), data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
s_cohon2 <- solve(p_cohon2, force=TRUE)
saveRDS(s_cohon2, glue('{out_folder}.rds'))
evaluate_solution(p_cohon2, s_cohon2)
# This is a good one to hang out to to show progression.

out_folder <- 's09_boundary_cohon3'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)
total_cost <- sum(pus$COST) + 1000
p_cohon3 <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=(cohon_penalty[1]/1000), data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
s_cohon3 <- solve(p_cohon3, force=TRUE)
saveRDS(s_cohon3, glue('{out_folder}.rds'))
evaluate_solution(p_cohon3, s_cohon3)

##################
# try smaller boundary penalty
out_folder <- 's10_boundary'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)
total_cost <- sum(pus$COST) + 1000
p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_shortfall_objective(total_cost) %>%
  add_relative_targets(tz) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=0.000001, data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=21600)
sink('log.txt')
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)

# This one looks good for now.



##############################################################
# solution 11:
# 
# Test reducing the targets and using min_set_objective and
# see how different it is to min_shortfall
# I'm adjusting based on the shortfalls in solution #10
##############################################################
out_folder <- 's11_minsetobj'
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

# adjust targets
tz2 <- tz
tz2['hu_ot_dredgingsites',2] <- 0.0
tz2['hu_tr_vesseltraffic_d',3] <- 0.28
tz2['hu_co_fishing_dive_d',4] <- 0.49
tz2['hu_ot_underwaterinfrastructure',3]<-0.62
tz2['hu_rf_fishing_groundfish',4]<-0.63
tz2['hu_tr_portsandterminals',2]<-0.64
tz2['hu_rf_fishing_crab',4]<-0.64
tz2['hu_ot_citypopulation_d',2]<-0.65
tz2['hu_rf_fishing_prawnandshrimp',4]<-0.68
tz2['hu_ot_floatingstructures',2]<-0.69
tz2['hu_co_fishing_trap_d',4]<-0.69
tz2['hu_co_fishing_gill_salmon_d',4]<-0.69
tz2['hu_ot_loghandlingstorage',2]<-0.69
tz2['hu_rf_fishing_anadromous',4]<-0.69
tz2['eco_fish_salmon',1]<-0.39
tz2['hu_co_fishing_seine_salmon_d',4]<-0.69
tz2['eco_fish_sixgillshark',1]<-0.39


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

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz2) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_boundary_penalties(penalty=0.000001, data=pus_bd, zone=zb, edge_factor=c(0.5, 0.5, 0.5, 0.5)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.01, threads = parallel::detectCores(TRUE)-1, time_limit=3600)
sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)

#!!!!
# So this is way different. min_shortfall was doing what I thought it was doing.
# Basically, when I give it all that budget, it seems to use most of it up to 
# exceed targets, so the overall solution is going to look quite different.
# See Jeff’s Jan 25 answer here: https://github.com/prioritizr/prioritizr/issues/224
# The thing on the shopping list and his point 2.
# Basically, min_shortfall can allow me to see that targets that aren’t met, 
# but it is going to spend all that extra budget  to exceed targets.

#############


#############
# same thing but without boundary penalty. I will compare this to Marxan to see
# if we get about the same thing
out_folder <- 's12_minsetobj_nobound'
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
# tz2['hu_co_fishing_trap_d',4]<-0.69
# tz2['hu_co_fishing_gill_salmon_d',4]<-0.69
# tz2['hu_ot_loghandlingstorage',2]<-0.69
# tz2['hu_rf_fishing_anadromous',4]<-0.69
# tz2['eco_fish_salmon',1]<-0.39
# tz2['hu_co_fishing_seine_salmon_d',4]<-0.69
# tz2['eco_fish_sixgillshark',1]<-0.39


# Create zones object
z <- zones(
  features, features, features, features,
  zone_names=zone_names,
  feature_names=features
)

p <- problem(pus, z, cost_column=c('COST', 'COST', 'COST', 'COST')) %>%
  add_min_set_objective() %>%
  add_relative_targets(tz2) %>%
  add_manual_locked_constraints(locked_df) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1, time_limit=3600)
sink('log.txt', split=TRUE)
s <- solve(p, force=TRUE)
sink()
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s)







##############################################################
# Miscelanous stuff
# 
##############################################################

#### runtime
#rt <- attr(s, 'runtime')
#### if it reached optimality or if it reached the time limit I set:
#st <- attr(s, 'status')
#### There's no way to save the gap acheived, so you need to note this after a run.


#### Calculate irreplaceability
# This is so slow so only run it on
# the final solution, if at all.

# gap_replace <- 95 # set super high for testing
# out_folder <- 's03_zones2'
# out_dir <- file.path(out_path, out_folder)
# setwd(out_dir)
# sols <- grep('solution_1', colnames(s), value=TRUE) # get all zones of 1st solution
# rc1 <- p %>%
#   add_gurobi_solver(gap=gap_replace, verbose=TRUE, threads = parallel::detectCores(TRUE)-1) %>%
#   eval_replacement_importance(s[,sols])
# saveRDS(rc1, glue('{out_folder}_replaceImp.rds'))
# rc1 <- tibble::rowid_to_column(rc1, 'id_index') # add id from index for join
# rc1 <- as.data.frame(rc1)
# rc1 <- select(rc1, -geometry) # drop shape field
# s2 <- tibble::rowid_to_column(s, 'id_index')
# s2 <- left_join(s2, rc1, by='id_index')
# s2 <- as_tibble(s2)
# s2 <- select(s2, -'Shape')
# write.csv(s2, 'irreplaceability.csv', row.names=FALSE)










