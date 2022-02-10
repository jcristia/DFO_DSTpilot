
# just get some data in and refamiliarize yourself with Prioritizr

library(prioritizr)
library(dplyr)
library(glue)


########################################
# paths and data
########################################

root <- r'(C:\Users\jcristia\Documents\GIS\DFO\DST_pilot)'
gdb <- file.path(root, 'spatial/03_working/dst_grid.gdb')
fc <- 'dst_grid1km_PUVSP'
out_path <- file.path(root, 'scripts/dst_02_prioritizr/outputs')

# read in feature class
# contains Planning Units and area of each feature in each PU
pus <- st_read(gdb, layer = fc)

# make locked in/out fields logical TRUE/FALSE (fgdb doesn't have a logical data type)
#pus$locked_in <- as.logical(pus$locked_in)
#pus$locked_out <- as.logical(pus$locked_out)



########################################
# function to evaluate each solution and
# calculate irreplaceability
########################################

evaluate_solution <- function(p1, s1, eval_replace) {
  
  # evaluate the performance of the first solution
  eval_numsel <- eval_n_summary(p1, s1[,'solution_1'])
  eval_cost <- eval_cost_summary(p1, s1[,'solution_1'])
  eval_targets <- eval_target_coverage_summary(p1, s1[,'solution_1'])
  eval_boundary <- eval_boundary_summary(p1, s1[,'solution_1'])
  
  # set up table
  s2 <- tibble::rowid_to_column(s1, 'id_index')
  
  # calculate importance (irreplaceability) of the first solution
  # This is SLOW. Maybe don't run this every time.
  if (eval_replace == TRUE){
    rc1 <- p1 %>%
      add_gurobi_solver(gap=0, verbose=FALSE) %>%
      eval_replacement_importance(s1[,'solution_1'])
    rc1 <- tibble::rowid_to_column(rc1, 'id_index') # add id from index for join
    rc1 <- as.data.frame(rc1)
    rc1 <- rc1[c('id_index', 'rc')] # drop shape field
    s2 <- left_join(s2, rc1, by='id_index')
  }
  
  # output solutions, target evaluation, and summary evaluation to csv
  s2 <- as_tibble(s2)
  s2 <- select(s2, -'Shape')
  write.csv(s2, 'solution.csv', row.names=FALSE)
  write.csv(eval_targets, 'eval_targets.csv', row.names=FALSE)
  write.csv(eval_boundary, 'eval_boundary.csv', row.names=FALSE)
  df_summ <- tibble(total_selected=eval_numsel$cost, cost=eval_cost$cost)
  write.csv(df_summ, 'eval_summary.csv', row.names=FALSE)
  
}


##############################################################
# solution 1:
# 
# only eco features
# equal targets
##############################################################

out_folder <- 's01_base'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

features <- names(pus)[grepl('mpatt_eco', names(pus))]
targets <- 0.1

p <- problem(pus, features, cost_column='COST') %>%
  add_min_set_objective() %>%
  add_relative_targets(targets) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0, threads = parallel::detectCores(TRUE)-1) %>%
  add_gap_portfolio(number_solutions = 1, pool_gap = 0.0)
  #add_cbc_solver(gap=0) %>% # says this is fastest open source solver
  #add_lpsymphony_solver(gap=0) %>% # however, this performs well and looks cleaner
  #add_cuts_portfolio(number_solutions = 2) # when not using gurobi
s <- solve(p)
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s, eval_replace=FALSE)



##############################################################
# solution 2:
# 
# Add boundary penalty
# 
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
##############################################################

out_folder <- 's02_boundary10e7'
out_dir <- file.path(out_path, out_folder)
dir.create(out_dir)
setwd(out_dir)

features <- names(pus)[grepl('mpatt_eco', names(pus))]
targets <- 0.1

p <- problem(pus, features, cost_column='COST') %>%
  add_min_set_objective() %>%
  add_relative_targets(targets) %>%
  add_binary_decisions() %>%
  # test values starting with 10e-7 to get a sense of the sensitivity
  add_boundary_penalties(penalty=0.0000001, edge_factor=0.5) %>%
  add_gurobi_solver(gap = 0.01, threads = parallel::detectCores(TRUE)-1) %>%
  add_gap_portfolio(number_solutions = 1, pool_gap = 0.01)
s <- solve(p)
saveRDS(s, glue('{out_folder}.rds'))
evaluate_solution(p, s, eval_replace=FALSE)




##############################################################
# solutions 3-n:
# 
# Trying to speed up boundary penalty stuff and figure out
# how to properly do sensitivity analysis.
#
# Things to do:
# Create, scale and save boundary matrix.
# Do "sweet spot" test to find starting point.
# read in previous solution as a starting point.
##############################################################

pus_bd <- boundary_matrix(pus) # Then scale this. So that it is in the same
# order of magnitude as the costs (so 1-100).
# Perhaps for the "previous" solution it should first be the sweet spot one
# I use as a starting point.
pres_solu <- readRDS(glue('{out_folder}.rds'))

# perfect, there is code for picking BLM values as Marxan suggests:
# https://prioritizr.net/articles/calibrating_trade-offs_tutorial.html#cohon-et-al--1979-method

# Then, generate some solutions around this BLM values.
# see code in prioritizr calibrating trade-offs vignette for how to do
# multiple runs with lapply
























