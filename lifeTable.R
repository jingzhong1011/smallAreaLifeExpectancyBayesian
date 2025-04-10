##################
### Life Table ###
##################

# For now it proves as a robust and faster version for life table,
# reducing for about 3/4 time.
# It use data.table for boosting process, suggest using parallel computation.

library(data.table)

leTable_dt <- function(data_dt, year_val, sex_val, moi_val){
  if(!is.data.table(data_dt)){
    data_dt <- as.data.table(data_dt)
  }
  # Filter
  filtered_dt <- data_dt[YEAR == year_val & SEX == sex_val & MOI == moi_val,
                         .(SEX, YEAR, MOI, AGP, COUNT, POP)]
  return(filtered_dt)
}

### Kannisto-Thatcher life table extension in high age ###
ktExtension <- function(lx70){
  # make sure input of lx70 is a vector
  # For age groups >= 70, calculate hazard rate using the approximation
  # mu(x+1/2) ~ -log(1 - q) = -log(p) where p is the probability of
  # survival to the next age group and equals l(x+5) / lx

  mux <- (log(lx70) - log(c(lx70[-1], NA))) / 5
  mux <- case_when(is.na(mux) ~ 0, T ~ mux)
  
  valid_mux_indices <- !is.na(mux)
  mux_valid <- mux[valid_mux_indices]
  lx70_valid_for_mux <- lx70[valid_mux_indices]
  
  # # Check if the data is enable to do following steps
  # if(length(lx70) < 4 || length(mux_valid) < 3){
  #   stop("Insufficient data points (lx70) for ktExtension.")
  # }
  
  # Calculate lx for 1-year age groups from 70 to 84 (or up to available data)
  # The original code assumes 4 groups (70-74, 75-79, 80-84, 85+)
  # Let's replicate the structure assuming 4 groups input
  # Rep lengths depend on having 3 intervals for mux (70-74 -> 75-79, etc.)
  
  # Generate AGP sequence
  agp_seq <- 70:84 # Original code goes up to 85 for calculation setup
  
  # Create the data frame/data table structure
  ktData_list <- list(
    AGP = agp_seq,
    lx = lx70[rep(1: 3, each = 5)], # Replicate lx for each 5-year group up to 84
    int = rep(0:4, 3),             # Interval within the 5-year group
    mux = rep(mux_valid[1:3], each = 5) # Replicate mux
  )
  
  # Add the starting point for age 85+
  ktData_list$AGP <- c(ktData_list$AGP, 85)
  ktData_list$lx <- c(ktData_list$lx, lx70[4]) # lx at start of 85+ group
  ktData_list$int <- c(ktData_list$int, 0)     # Interval is 0 for the start of 85+
  ktData_list$mux <- c(ktData_list$mux, NA)    # Mux is not applicable here
  
  ktData <- as.data.table(ktData_list) # Convert to data.table
  
  # Calculate single-year lx up to 84, and qx / logitqx
  ktData[, lx := ifelse(AGP == 85, lx, lx * exp(mux * -int))] # Calculate lx within 5-year groups
  ktData[, dx := lx - shift(lx, type = "lead")] # dx = lx - l(x+1)
  ktData[AGP == 85, dx := NA] # dx for the last age (85) needs the next lx (from model)
  ktData[, qx := dx / lx]
  ktData[, logitqx := log(qx / (1 - qx))]
  ktData <- ktData[AGP < 85] # Keep only data up to age 84 for regression
  
  # Run regression on logit of probability of dying by logit
  y <- ktData$logitqx
  x <- ktData$AGP + 0.5 # Age at midpoint
  
  # Check for NA/Inf values before regression
  valid_reg <- !is.na(y) & is.finite(y) & !is.na(x) & is.finite(x)
  if(sum(valid_reg) < 2){
    stop("Not enough valid data points for logit qx regression in ktExtension.")
  }
  
  # Used weights = x
  mod <- lm(y[valid_reg] ~ x[valid_reg], weights = x[valid_reg])
  logA <- coef(mod)[1]
  B <- coef(mod)[2]
  
  # Calculate qx for age x >= 85
  # Expand to age 129
  ages_extrapolate <- 85: 129
  logitqx85 <- logA + B * (ages_extrapolate + 0.5)
  qx85 <- exp(logitqx85) / (1 + exp(logitqx85))
  
  # Prepare data table for ages 85+
  ktData85 <- data.table(
    AGP = ages_extrapolate,
    lx = NA_real_, # Initialize lx column
    qx = qx85
  )
  ktData85[1, lx := lx70[4]] # Starting lx for age 85
  
  # Calculate lx and dx for age x >= 85 using data.table efficiently
  # Calculate px first
  ktData85[, px := 1 - qx]
  
  # Calculate lx by iteration
  for(i in 2: nrow(ktData85)){
    ktData85$lx[i] <- ktData85$lx[i - 1] * ktData85$px[i - 1]
  }
  # Calculate dx
  ktData85[, dx := lx * qx]
  
  # Combine data < 85 and >= 85
  # We only need the calculated dx for the final ax calculation
  # The structure needs dx for ages 70-129
  dx_full <- c(ktData$dx, ktData85$dx) # Combine dx values
  agp_full <- 70: 129
  
  # Calculate ax for the original 5-year age groups (70-74, 75-79, 80-84, 85+)
  # Define the original age groups for aggregation
  agp_groups <- c(70, 75, 80, 85)
  # Assign each single year age (70-129) to its corresponding group start age
  group_map <- agp_groups[findInterval(agp_full, agp_groups)]
  
  # Create a temporary data table for aggregation
  temp_ax_calc <- data.table(agp_full, dx_full, group_map)
  
  # Midpoint calculation adjustment: yl = age - group_start + 0.5
  temp_ax_calc[, yl := agp_full - group_map + 0.5]
  
  # Calculate sum(dx * yl) and sum(dx) per group
  ax_calc <- temp_ax_calc[, .(sum_dx_yl = sum(dx_full * yl, na.rm = TRUE),
                              sum_dx = sum(dx_full, na.rm = TRUE)), by = group_map]
  
  # Calculate ax = sum(dx * yl) / sum(dx)
  ax_calc[, ax := ifelse(sum_dx == 0, 2.5, sum_dx_yl / sum_dx)] # Use 2.5 if sum(dx) is 0 (or handle as needed)
  
  # Return the ax values for the 70-74, 75-79, 80-84, 85+ groups
  # Ensure the order matches the expected output (corresponding to lx70 input)
  # The groups correspond to indices 16, 17, 18, 19 in the main leCalc data
  # Need to return ax for 70-74, 75-79, 80-84, 85+ (4 values)
  final_ax <- ax_calc[order(group_map)]$ax
  
  # The result should be 4 values
  if(length(final_ax) != 4){
    warning("ktExtension did not produce exactly 4 ax values. Check input lx70 length and logic.")
    return(rep(2.5, 4))
  }
  
  return(final_ax)
}


### Young age ax calculation ###
# This function is simple, no significant optimization needed
youngExtension <- function(mx5, sex) {
  # p.48 Preston et al, Demography: measuring and modeling population processes,
  # "The lower the mortality, the more heavily will infant deaths be
  # concentrated at the earliest stages of infancy". The contingency table
  # is based on Coale and Demeny (1983), who fitted a line to international
  # and intertemporal data
  ax5 <- numeric(2) # Pre-allocate vector
  m0 <- mx5[1] # Infant mortality rate (age 0)
  
  # Use ifelse for vectorized condition check (though only one value is checked)
  if(sex == 1){ # Male
    ax5[1] <- ifelse(m0 >= 0.107, 0.330, 0.045 + 2.684 * m0)
    ax5[2] <- ifelse(m0 >= 0.107, 1.352, 1.651 - 2.816 * m0)
  }else if(sex == 2) { # Female
    ax5[1] <- ifelse(m0 >= 0.107, 0.350, 0.053 + 2.800 * m0)
    ax5[2] <- ifelse(m0 >= 0.107, 1.361, 1.522 - 1.518 * m0)
  }else{
    stop("Invalid sex indicator.")
  }
  return(ax5)
}



### LE calculation using data.table ###
leCalc_dt <- function(data_dt){
  # Ensure input is a data.table
  if(!is.data.table(data_dt)){
    data_dt <- as.data.table(data_dt)
  }
  
  sex_val <- data_dt$SEX[1]
  n_rows <- nrow(data_dt) # Should be 19
  
  # Define nx values directly
  nx_vals <- c(1, 4, rep(5, n_rows - 3), NA)
  
  data_dt[, id := .I] # Row number
  data_dt[, ax := c(0.1, 0.4, rep(2.5, n_rows - 2))] # Initial ax guess
  data_dt[, nx := nx_vals]
  
  # Calculate mx, handle POP < 1 or mx == 0
  data_dt[, mx := ifelse(POP < 1, 1, COUNT / POP)]
  data_dt[mx == 0, mx := 1e-5] # Handle zero mortality case
  
  # Calculate initial qx, px, lx, dx
  data_dt[, qx := ifelse(AGP == 85, 1, (nx * mx) / (1 + (nx - ax) * mx))]
  data_dt[qx > 1, qx := 0.99999]
  data_dt[is.na(qx), qx := 1] # Ensure last qx is 1 if nx is NA
  data_dt[, px := 1 - qx]
  data_dt[, lx := 1e5] # Initialize lx
  # Calculate lx using cumulative product of lagged px
  # data_dt[id > 1, lx := 1e5 * cumprod(shift(px, fill=1)[id > 1])] # Careful with indexing
  px_shifted <- shift(data_dt$px, n=1, fill=1) # Get lagged px, fill first with 1
  lx_calc <- 1e5 * cumprod(px_shifted)
  data_dt[, lx := lx_calc]
  data_dt[, dx := lx * qx]
  
  # Apply extensions for ax
  lx70_vals <- data_dt[AGP >= 70, lx]
  mx5_vals <- data_dt[AGP < 5, mx]
  
  # Ensure enough values for ktExtension
  if(length(lx70_vals) >= 4){
    ax_kt <- ktExtension(lx70_vals)
    data_dt[AGP %in% c(70, 75, 80, 85), ax := ax_kt] # Assign the 4 ax values
  }else{
    warning("Not enough age groups >= 70 for ktExtension. Using default ax=2.5.")
    data_dt[16:19, ax := 2.5]
  }
  
  # Ensure enough values for youngExtension
  if(length(mx5_vals) >= 2){
    ax_young <- youngExtension(mx5_vals, sex_val)
    data_dt[AGP %in% c(0, 1), ax := ax_young] # Assign the 2 ax values
  }else{
    warning("Not enough age groups < 5 for youngExtension. Using default ax.")
  }
  
  
  ### Iteration using data.table ###
  iter <- 4 # Number of iterations

  
  for(i in seq(iter)){
    # Store ax from previous iteration to calculate ax.new
    nx <- data_dt$nx
    ax_old <- data_dt$ax
    dx_old <- data_dt$dx # Need dx from previous step for ax calculation
    
    # Calculate ax.new for ages 5 to 80 (indices 3 to 18)
    # Need dx[k-1], dx[k], dx[k+1]
    ax_new <- ax_old # Start with old values
    for(k in 3:(n_rows - 1)){ # Indices 3 to 18
      # Check bounds for dx indices
      dx_km1 <- if (k > 1) dx_old[k - 1] else 0 # Handle k = 3 case (k - 1 = 2)
      dx_k <- dx_old[k]
      dx_kp1 <- if (k < n_rows) dx_old[k + 1] else 0 # Handle k = 18 case (k + 1 = 19)
      
      # Avoid division by zero if dx_k is 0
      if(dx_k != 0){
        ax_new[k] <- (-5 / 24 * dx_km1 + 2.5 * dx_k + 5 / 24 * dx_kp1) / dx_k
        ax_new[k] <- pmax(pmin(ax_new[k], nx[k]), 0)
      }else{
        ax_new[k] <- 2.5 # Default if dx is zero
      }
    }
    
    # Keep the ax values from extensions for young and old ages unchanged
    ax_new[data_dt$AGP %in% c(0, 1)] <- data_dt$ax[data_dt$AGP %in% c(0, 1)]
    ax_new[data_dt$AGP %in% c(70, 75, 80, 85)] <- data_dt$ax[data_dt$AGP %in% c(70, 75, 80, 85)]
    
    # Update ax
    data_dt[, ax := ax_new]
    
    # Recalculate qx, px, lx, dx based on new ax
    # Note: nx is involved, use the nx column
    data_dt[, qx := (nx * mx) / (1 + (nx - ax) * mx)]
    # Handle specific cases for qx
    data_dt[AGP == 0, qx := (1 * mx) / (1 + (1 - ax) * mx)]
    data_dt[AGP == 1, qx := (4 * mx) / (1 + (4 - ax) * mx)]
    data_dt[qx > 1, qx := 0.999999] # Cap qx at slightly less than 1
    data_dt[AGP == 85, qx := 1]     # Last age group qx = 1
    data_dt[is.na(qx), qx:= 1] # Ensure NA qx (e.g. last row if nx was NA) is 1
    
    data_dt[, px := 1 - qx]
    
    # Recalculate lx
    px_shifted <- shift(data_dt$px, n=1, fill=1)
    lx_calc <- 1e5 * cumprod(px_shifted)
    data_dt[, lx := lx_calc]
    
    # Recalculate dx
    data_dt[, dx := lx * qx]
    
    # iterData[[i + 1]] <- data_dt[, .(AGP, ax, qx, lx, dx)]
  }
  
  # After iterations, calculate Lx, Tx, ex
  # Lx = n*l(x+n) + ax*dx  (general formula)
  # Lx[1] = 1*lx[2] + ax[1]*dx[1]
  # Lx[2] = 4*lx[3] + ax[2]*dx[2]
  # Lx[k] = 5*lx[k+1] + ax[k]*dx[k] for k=3,...,18
  # Lx[19] = lx[19] / mx[19]  (or lx[19] * ax[19] if ax represents avg lifetime in last interval)

  data_dt[, lx_lead := shift(lx, type = "lead")]
  
  data_dt[, Lx := 0]
  data_dt[id == 1, Lx := nx * lx_lead + ax * dx]
  data_dt[id == 2, Lx := nx * lx_lead + ax * dx]
  data_dt[id >= 3 & id < n_rows, Lx := nx * lx_lead + ax * dx]
  
  # Handle last row (id=19, AGP=85)
  # Ensure mx[19] is not zero
  last_mx <- data_dt[id == n_rows, mx]
  data_dt[id == n_rows, Lx := ifelse(last_mx == 0, 0, lx / last_mx)] # Avoid division by zero
  
  # Calculate Tx (cumulative sum of Lx from bottom up)
  data_dt[, Tx := rev(cumsum(rev(Lx)))]
  
  # Calculate ex
  data_dt[, ex := ifelse(lx == 0, 0, Tx / lx)] 
  
  # Calculate Standard Error and Confidence Intervals
  # spx = variance of qx
  # W_spx = contribution to variance of ex[0]
  # STx = variance of Tx
  # SeSE = standard error of ex
  
  data_dt[, spx := fifelse(id == n_rows,
                           4 / (COUNT * (mx ^ 2)), # Approximation for last interval variance
                           fifelse(COUNT == 0, 0, (qx ^ 2) * (1 - qx) / COUNT))]
  
  data_dt[, spx := fifelse(id == n_rows,
                           fifelse(COUNT == 0 | mx == 0, 0, 4 / (COUNT * (mx^2))), # Handle COUNT=0 or mx=0
                           fifelse(COUNT == 0, 0, (qx^2) * (1 - qx) / COUNT))]
  
  # Need ex_lead (ex[x+1]) for W_spx calculation
  data_dt[, ex_lead := shift(ex, type = "lead", fill = 0)]
  
  data_dt[, W_spx := fifelse(id < n_rows,
                             spx * (lx^2) * (((1 - ax) * nx + ex_lead)^2),
                             ((lx / 2)^2) * spx)]
  
  data_dt[, STx := rev(cumsum(rev(W_spx)))] # Variance of Tx (cumulative from bottom)
  data_dt[, SeSE := fifelse(lx == 0, 0, sqrt(STx / (lx^2)))] # Standard Error of ex
  
  # Confidence intervals
  z <- 1.96 # For 95% CI (Original code had z=0.95, which is likely wrong, z=1.96 for 95%)
  # Let's use 1.96 assuming 95% CI is desired.
  
  z <- qnorm(1 - (1 - 0.95) / 2) # Standard way to get z for 95% CI
  
  data_dt[, lowercl := ex - z * SeSE]
  data_dt[, uppercl := ex + z * SeSE]
  
  # Ensure lower CI is not negative
  data_dt[lowercl < 0, lowercl := 0]
  
  # Prepare output table
  # Add age group labels
  agp_labels <- c("0", "1-4", paste0(seq(5, 80, by = 5), "-", seq(9, 84, by = 5)), "85+")
  data_dt[, agp := agp_labels]
  # Select and order columns for the final output
  output_dt <- data_dt[, .(agp, age=AGP, deaths=COUNT, pops=POP, nx, ax, mx, qx, dx, lx, Lx, Tx, ex, lowercl, uppercl)]
  
  return(output_dt)
}


### LE compute by loop wrapper (using data.table functions) ###
leCompute_dt <- function(data_dt, year_val, sex_val, moi_val){
  # Ensure input is data.table
  if(!is.data.table(data_dt)){
    data_dt <- as.data.table(data_dt)
  }
  lt_dt <- leTable_dt(data_dt, year_val, sex_val, moi_val)
  
  # Check if filtering resulted in data
  if(nrow(lt_dt) == 0){
    warning(paste("No data found for YEAR =", year_val, "SEX =", sex_val, "MOI =", moi_val))
    return(NULL) # Return NULL or an empty data.table structure
  }
  # Check if lt_dt has the expected 19 rows
  if(nrow(lt_dt) != 19){
    warning(paste("Data for YEAR =", year_val, "SEX =", sex_val, "MOI =", moi_val,
                  "does not have the expected 19 age groups. Skipping calculation."))
    return(NULL)
  }
  
  result_dt <- leCalc_dt(lt_dt)
  result_dt[, YEAR := year_val]
  result_dt[, SEX := sex_val]
  result_dt[, MOI := moi_val]
  
  output_dt <- result_dt[, .(YEAR, SEX, MOI, agp, ax, dx, mx, lx, Lx, Tx, ex, lowercl, uppercl)]
  
  return(output_dt)
}

