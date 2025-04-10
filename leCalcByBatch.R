####################
### TESTING CODE ###
####################

#####################################################################################
### It will be an OFFICIALLY code for running LIFE TABLE by CHAINS, DRAWS in MCMC ###
#####################################################################################
library(dplyr)
library(data.table)
library(ncdf4)
library(purrr)

setwd("/Users/jing-zhongwang/Library/CloudStorage/GoogleDrive-clockwcc@gmail.com/我的雲端硬碟/NTUEPM/Master Thesis/Paper Writing/Data")

DeathCount <- fread("DeathCount_v5.csv") %>% 
  filter(! Urbanicity_AST2023 == 9) %>% 
  rename(COUNT_RAW = COUNT)
DeathCount_male <- DeathCount %>% filter(SEX == 1)
PPmale <- nc_open("smoothed/PPmale.nc")
PPd_male <- ncvar_get(PPmale, "posterior_predictive/d")
PPd_male


library(future.apply) # Parallel computation
library(parallel) # Detecting Cores

# Set parallel environment
no_cores <- detectCores() - 1
if (no_cores < 1) no_cores <- 1
plan(multisession, workers = no_cores) # background R session
# plan(multicore) # if in Linux or server

if(!is.data.table(DeathCount_male)){
  DeathCount_male_dt <- as.data.table(DeathCount_male)
}else{
  DeathCount_male_dt <- copy(DeathCount_male)
}


year <- 2000:2021
sex <- 1
moi <- levels(as.factor(DeathCount_male_dt$MOI))
Combinations <- data.table(expand.grid(YEAR = year, SEX = sex, MOI = moi))


n_rows_ppd <- dim(PPd_male)[1]
n_draws <- dim(PPd_male)[2]
n_chains <- dim(PPd_male)[3]
if(nrow(DeathCount_male_dt) != n_rows_ppd){
  stop("Number of rows in DeathCount_male does not match first dimension of PPd_male.")
}


# MAIN LOOP
LEsmoothed_list <- list()
##################
#For Testing
n_chains = 4
n_draws = 2
##################
for(chain in 1: n_chains){ 
  for(draw in 1: n_draws){ 
    cat("CHAIN =", chain, "DRAW =", draw, "\n")
    startTime <- Sys.time()
    
    current_counts <- PPd_male[1: 149644, draw, chain]
    test_dt <- copy(DeathCount_male_dt)
    test_dt[, COUNT := current_counts]
    
    # Using `future_lapply` to compute parallel
    # `future_lapply` passes each line of Combinations to the function
    # We need a helper function to receive these parameters and call leCompute_dt
    compute_for_row <- function(idx, combo_dt, data_to_use){
      row <- combo_dt[idx, ]
      result <- leCompute_dt(data_to_use, row$YEAR, row$SEX, row$MOI) %>% 
        mutate(CHAIN = chain, DRAW = draw)
      return(result)
    }
    
    # DO leCompute_dt
    results_list <- future_lapply(1:nrow(Combinations),
                                  compute_for_row,
                                  combo_dt = Combinations,
                                  data_to_use = test_dt)
    
    # Use rbindlist
    valid_results <- Filter(Negate(is.null), results_list)
    if(length(valid_results) > 0){
      LEsmoothed <- rbindlist(valid_results)
      
      LEsmoothed[, CHAIN := chain]
      LEsmoothed[, DRAW := draw]
      
      # Save the results
      LEsmoothed_list[[paste0("chain_", chain, "_draw_", draw)]] <- LEsmoothed
    } else {
      cat("  WARN: No valid results produced for chain", chain, "draw", draw, "\n")
    }
    results_list <- NULL
    
    costTime <- Sys.time() - startTime
    cat("  Time taken:", costTime, "\n")
  }
}

# Back to non-parallel environment
plan(sequential)


##############################################
##############################################
# LEsmoothed_list <- list()
# year <- 2000:2021
# sex <- 1
# moi <- levels(as.factor(DeathCount$MOI))
# Combinations <- expand.grid(YEAR = year, SEX = sex, MOI = moi)
# for(chain in 1: 2){
#   for(draw in 1: 10){
#     cat("CHAIN =", chain, "DRAW =", draw, "\n")
#     startTime <- Sys.time()
#     test <- DeathCount_male %>% mutate(COUNT = PPd_male[1: 149644, draw, chain])
#     LEsmoothed <- do.call(rbind, apply(Combinations, 1, function(row) leCompute(test, row["YEAR"], row["SEX"], row["MOI"])))
#     LEsmoothed <- LEsmoothed %>% mutate(CHAIN = chain, DRAW = draw)
#     LEsmoothed_list[[paste0("chain_", chain, "_draw_", draw)]] <- LEsmoothed
#     costTime <- Sys.time() - startTime
#     cat(costTime, "\n")
#   }
# }
# for(chain in 1: 1){
#   for(draw in 1: 100){
#     for(i in 2000: 2021){
#       for(j in 1: 2){
#         for(k in moi){
#           cat("YEAR =", i, "SEX =", j, "MOI =", k, "\n")
#           leCompute_dt(DeathCount_male %>% mutate(COUNT = PPd_male[1: 149644, draw, chain]), i, j, k)
#         }
#       }
#     }
#   }
# }
# 
# data_dt <- leTable_dt(DeathCount_male %>% mutate(COUNT = PPd_male[1: 149644, 9, 1]), 2000, 1, 64000360)
# lx70 <- lx70_vals

