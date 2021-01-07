#NA functions

#LF trick  function
nudge_zero <- function(x){
  if(identical(x,0)){
    x <- 0.1
  }
  return(x)
}

#function to do NA conversion
zero_NA <- function(x){
  x[x == 0] <- NA
  x
}

#functions from July 2020 to force monotonicity starting from cum series
make_vec1 <- function(x){
  x_out <- x -lag(x)
  
  x_out[1] <- x[1]
  
  return(x_out)
}

# A function factory for getting integer y-axis values.
# from: https://joshuacook.netlify.com/post/integer-values-ggplot-axis/
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


#function to find index that marks first sequence of length_use values.  Default length = 8 per Lloyd Provost 30 March 2020 
index_test <- function(x,
                       index,
                       length_use=8){
  x_check <- x[index:(index+length_use - 1)]
  if(all(x_check>0)){
    use_seq <- TRUE
    index_use <- index
  } else {
    use_seq <- FALSE
    index_use <- index + 1
  }
  return(list(use_seq,index_use))
}

#function to detect large data dumps
detect_outlier_dates <- function(
  data,
  ratio_limit = 6)
{
  model <- loess(New_var ~ as.numeric(datex), data = data, span = 0.25)
  

  outlier_index <- data$New_var > max(10, median(data$New_var, na.rm = TRUE)) &
                   (data$New_var / model$fitted) >= ratio_limit
  
  data$datex[outlier_index]
}

#function to force monotonicity in nominally monotone non-decreasing series
force_monotonicity <- function(x) {
  
  sum_original <- sum(x, na.rm = TRUE)
  
  # spread negative values backwards in series so that cumulative sum never decreases
  negative_index <- which(x < 0)
  
  while (any(negative_index)) {
    
    index <- tail(negative_index, 1)
    
    value <- abs(x[index])
    
    transferred <- 0
    
    x[index] <- 0
    
    while (transferred < value) {
      
      index <- index - 1
      
      transfer <- min(x[index], value - transferred)
      
      x[index] <- x[index] - transfer
      
      transferred <- transferred + transfer
    }
    
    negative_index <- which(x < 0)
  }
  
  stopifnot(sum(x, na.rm = TRUE) == sum_original)
  
  x
}

#function to identify changes in phase
model_phase_change <- function(
  data,
  y)
{
  data[[y]][!is.na(data[[y]]) & data[[y]] <= 0] <- NA
  
  data$log_10_var <- log10(data[[y]])
  
  data$log_10_var[is.infinite(data$log_10_var)] <- NA

  data$serial_day <- as.numeric(difftime(data$datex, min(data$datex, na.rm = TRUE), units = 'days')) + 1
  
  repeat {
    
    lm_raw <- try(lm(log_10_var ~ serial_day,
                     data = data))

    if ('try-error' %in% class(lm_raw)) break
    else {
      residual_diff <- diff(lm_raw$residuals)

      median_moving_range <- median(abs(residual_diff), na.rm=TRUE)

      data <- merge(data,
                           cbind(serial_day = lm_raw$model$serial_day,
                                 fitted_value = lm_raw$fitted.values),
                           by = 'serial_day',
                           all.x = TRUE)

      distance <- abs(data$log_10_var - data$fitted_value) / median_moving_range

      data$fitted_value <- NULL
      
      break
    }
  }
    
  if ('try-error' %in% class(lm_raw)) {
    
    output <- list(
      lm = NULL,
      exp_growth = FALSE,
      median_moving_range = NULL)
    
  } else {
    
    serial_day_coefficient <- lm_raw$coefficients['serial_day']

    growth_sign <- sign(serial_day_coefficient)
    
    serial_day_pvalue <- try(summary(lm_raw)$coefficients['serial_day', 'Pr(>|t|)'])
    
    exp_growth <- !is.na(serial_day_pvalue) &&
                  serial_day_pvalue < 0.05
    
    mean_var <- mean(data$log_10_var, na.rm = TRUE)
    
    # If exponential growth (up or down) is not detected, 
    # re-define median moving range based on mean_var as midline
    if (!exp_growth) {
      
      residuals_mean <- data$log_10_var - mean_var
      
      median_moving_range <- median(abs(residuals_mean), na.rm=TRUE)
      
    }
    
    output <- list(
      lm = lm_raw,
      mean_var = mean_var,
      growth_sign = growth_sign,
      exp_growth = exp_growth,
      median_moving_range = median_moving_range)
  }
  
  output
}

#function to determine phases with control limits and center lines
find_phase_dates <- function(
  data,
  adjust,
  ghost = TRUE)
  
{
  #message(sprintf(' -- %s-level: finding phase dates for %s', data$level[1], data$state[1]))
  message(sprintf(' -- %s-level: finding phase dates for %s', data$var_name[1], data$NAME[1]))
  result <- try({
    
    data <- data[order(data$datex), ]
    
    # Remove negative var by iteratively reducing the counts (to floor of 0) 
    # on preceding days.
    #data$New_var <- force_monotonicity(data$New_var)
    
    # Detect 'ghost' points - only run for raw series
    # Note ghost argument is TRUE by default, because 1) we want it TRUE for raw series, and
    # 2) for adjusted series, we want it TRUE for the first iteration of the function,
    # where the phases get estimated using the raw series. On the recursive call to find_phase_dates,
    # ghost will be set to FALSE since we don't want to apply ghosting to the adjusted series.
    if (ghost) {
      
      data$New_var_Dump <- NA
      
      ghost_dates <- detect_outlier_dates(data[c('datex', 'New_var')])
    
      ghost_index <- data$datex %in% ghost_dates
      
      if (any(ghost_index)) {
        
        data$New_var_Dump[ghost_index] <- data$New_var[ghost_index]
        data$New_var[ghost_index]      <- NA
        
      }
    }
    
    phase_data <- list()
    
    if (any(data$New_var > 0, na.rm = TRUE)) {
      
      date_phase_start <- min(data$datex[data$New_var > 0], na.rm = TRUE)
      
      data_var <- data[data$datex >= date_phase_start, ]
      
      data_var$weekday <- lubridate::wday(data_var$datex)
      
      date_max <- max(data_var$datex, na.rm = TRUE)
       
      cumulative_var <- cumsum(ifelse(is.na(data_var$New_var),
                                         0,
                                         data_var$New_var))
      
      epoch <- 1
      
      # phase is the number of phases detected so far (always increasing)
      # epoch_phase is the number of phases detected in the current epoch
      # (increments on every phase change, except new epochs, when it resets to 1)
      phase <- 1
      epoch_phase <- 1
            
      if (any(cumulative_var >= 8)) {
        
        while (date_phase_start <= (date_max - 4)) {
          
          message(' -- Epoch ', epoch, ', phase ', phase, ' overall, phase ', epoch_phase, ' within epoch: ', date_phase_start)

          phase_index <- data_var$datex >= date_phase_start
          
          # need 8+ cumulative values in phase to estimate midline and ucl
          cumulative_var_phase <- cumsum(data_var$New_var[phase_index])
          parameter_date <- data_var$datex[phase_index][which.max(cumulative_var_phase >= 8)]
          
          parameter_index <- data_var$datex >= date_phase_start & 
            data_var$datex <= max(date_phase_start + 20, parameter_date)
          
          test_index <- data_var$datex >= parameter_date
          
          if (epoch %in% c(1, 4)) {
            
            date_check_original <-
              date_check <- min(max(date_phase_start + 4, parameter_date), date_max)

            while (date_check <= date_max) {
              
              midline_index <- data_var$datex >= date_phase_start & data_var$datex <= min(date_check, date_phase_start + 20)
              
              midline <- mean(data_var$New_var[midline_index], na.rm = TRUE)
              
              residuals <- data_var$New_var - midline
              
              observed_values <- data_var$New_var
              
              check_index <- data_var$datex >= date_check_original & data_var$datex <= date_check
              
              residuals_check <- (observed_values - midline)[check_index]
              
              mmr <- median(abs(diff(residuals_check)))
              
              distance_check <- abs(residuals_check) / mmr
              
              ucl <- midline + 3 * sqrt(midline)
              lcl <- midline - 3 * sqrt(midline)
              
              if (!is.na(lcl) && lcl <= 0) lcl <- NA
              
              # find two consecutive dates that values cross above ucl, take the earliest
              above_ucl <- observed_values[check_index] > ucl & !is.na(observed_values[check_index])
              above_ucl_streak <- ave(above_ucl, cumsum(!above_ucl), FUN = cumsum)
              
              # In phase 1, we only want the date the series first exceeds ucl
              # otherwise, we want the first day of two consecutive days exceeding the ucl
              if (epoch_phase == 1) above_ucl_streak_index <- head(which(above_ucl_streak == 1, 1))
              else above_ucl_streak_index <- head(which(above_ucl_streak == 2), 1) - 1
              
              date_above_ucl <- data_var$datex[check_index][above_ucl_streak_index]

              # find the first date that values are above or below midline for 8 days in a row
              above_midline <- observed_values[check_index] > midline & !is.na(observed_values[check_index])
              above_midline_streak <- ave(above_midline, cumsum(!above_midline), FUN = cumsum)
              date_above_midline <- data_var$datex[check_index][head(which(above_midline_streak == 8), 1)]

              below_midline <- observed_values[check_index] < midline & !is.na(observed_values[check_index])
              below_midline_streak <- ave(below_midline, cumsum(!below_midline), FUN = cumsum)
              date_below_midline <- data_var$datex[check_index][head(which(below_midline_streak == 8), 1)]
              
              date_phase_end <- min(date_above_ucl, date_above_midline, date_below_midline) - 1
              
              # if no phase end condition was detected, date_phase_end will be -Inf
              # if date_phase was detected (finite), break the loop
              if (is.finite(date_phase_end)) {

                break
                
              } else date_check <- date_check + 1
            }
            
            if (!is.finite(date_phase_end)) {              
              date_phase_end <- date_max
            }
            
          } else if (epoch %in% c(2, 3)) {
            
            if (phase_change_result$exp_growth) {
              
              data_var$serial_day <- as.numeric(difftime(data_var$datex, date_phase_start, units = 'days')) + 1
  
              midline <- predict(phase_change_result$lm, data_var)
              
              data_var$serial_day <- NULL
            
            } else {
              
              midline <- phase_change_result$mean_var
              
            }
            
            date_check_original <-
              date_check <- min(max(date_phase_start + 21, parameter_date), date_max)
            
            mmr <- phase_change_result$median_moving_range
            
            observed_values <- log10(data_var$New_var)
            
            while (date_check <= date_max) {
              
              ucl <- midline + (3.14 * mmr)
              lcl <- midline - (3.14 * mmr)
              
              # In epoch 3 with no exponential growth, midline is a constant 
              # (and so is mmr), so ucl and lcl are also length 1.
              # This causes problems when indexing by check_index, so for epoch 3, we expand the constant
              # out to a length equal to the number of rows in data_var. Not the best solution, but easier
              # than adding conditional logic in the streak checks though.
              if (epoch == 3 && length(midline == 1)) {
                ucl <- rep(ucl, nrow(data_var))
                lcl <- rep(lcl, nrow(data_var))
              }
              
              check_index <- data_var$datex >= date_check_original & data_var$datex <= date_check
              
              above_ucl <- observed_values[check_index] > ucl[check_index] & !is.na(observed_values[check_index])
              above_ucl_streak <- ave(above_ucl, cumsum(!above_ucl), FUN = cumsum)
              date_above_ucl <- data_var$datex[check_index][head(which(above_ucl_streak == 2), 1) - 1]
              
              below_ucl <- observed_values[check_index] < lcl[check_index] & !is.na(observed_values[check_index])
              below_ucl_streak <- ave(below_ucl, cumsum(!below_ucl), FUN = cumsum)
              date_below_ucl <- data_var$datex[check_index][head(which(below_ucl_streak == 2), 1) - 1]
              
              residuals_check <- (observed_values - midline)[check_index]
              
              streak_sign <- 0
              streak_length <- 0
            
              for (residual_i in seq_along(residuals_check)) {
                
                residual <- residuals_check[residual_i]
                
                if (is.finite(residual)) {
                  if (sign(residual) == streak_sign) {
                    streak_length <- streak_length + 1
                  } else {
                    streak_sign <- sign(residual)
                    streak_length <- 0
                  }
                }
                
                streak_found <- streak_length == 8 
                  
                if (streak_found) {
                  date_streak <- data_var$datex[check_index][residual_i]
                  break
                } else if (residual_i == length(residuals_check)) {
                  date_streak <- numeric(0)
                  break
                }
              }
              
              date_phase_end <- min(date_below_ucl, date_above_ucl, date_streak) - 1
              
              if (is.finite(date_phase_end)) {
                
                break
                
              } else date_check <- date_check + 1
            }
            
            if (!is.finite(date_phase_end)) {              
              date_phase_end <- date_max
            }
              
            phase_days <- as.numeric(difftime(date_phase_end, date_phase_start, units = 'day')) + 1
            
            if (length(midline) > 1) {
              midline <- head(midline[phase_index], phase_days)
            }
            
            ucl     <- midline + (3.14 * mmr)
            lcl     <- midline - (3.14 * mmr)
            
            midline  <- 10 ^ midline
            ucl      <- 10 ^ ucl
            lcl      <- 10 ^ lcl
          }
          
          phase_index <- data_var$datex >= date_phase_start & data_var$datex <= date_phase_end
          
          if (epoch %in% c(1, 4)) {
            
            observed_values_phase <- observed_values[phase_index]
            
          } else if (epoch %in% c(2, 3)) {
            
            observed_values_phase <- 10 ^ observed_values[phase_index]
            mmr <- 10 ^ mmr
            
          }
          
          residuals_phase <- (observed_values_phase - midline)
          
          distance_phase <- abs(residuals_phase) / mmr
          
          phase_parameters <- list(
            phase = phase,
            epoch = epoch,
            midline = midline,
            lcl = lcl,
            ucl = ucl,
            start = date_phase_start,
            end   = date_phase_end)
          
          phase_data[[phase]] <- phase_parameters

          if (date_phase_end < (date_max - 4)) {
            
            model_index <- data_var$datex > date_phase_end & data_var$datex <= min(date_phase_end + 21, date_max, na.rm = TRUE)
            
            phase_change_result <- model_phase_change(
              data = data_var[model_index, ],
              y = 'New_var')
            
            phase <- phase + 1
            epoch_phase <- epoch_phase + 1
            
            if (phase_change_result$exp_growth && phase_change_result$growth_sign == 1) {
           
              if (epoch != 2) {
                epoch <- 2
                epoch_phase <- 1
              }
              
            } else {

              if (epoch == 2) {
                
                epoch <- 3
                epoch_phase <- 1
                
              } else if (epoch == 3) {
              
                # Epoch 4 requires:
                # 1. current phase LCL < 2, AND
                #   2a. phase change due to points below UCL, OR
                #   2b. phase change due to streak below midline (residual sign -1)
                # This excludes phase changes due to an uptick in cases.
                epoch_4 <- phase_parameters$lcl < 2 &&
                          (length(date_below_ucl) == 1 ||
                          (length(date_streak) == 1 && streak_sign == -1))
                
                if (epoch_4) {
                  epoch <- 4
                  epoch_phase <- 1
                }
              }
            }
          }
          
          date_phase_start <- date_phase_end + 1
        }
      } else {
        
        phase_data[[phase]] <- list(
          phase = phase,
          epoch = epoch,
          midline = NA,
          lcl = NA,
          ucl = NA,
          start = date_phase_start,
          end   = max(data_var$datex))
        
      }
    }
    
    for (phase_parameters in phase_data) {
      
      index <- data$datex >= phase_parameters$start & data$datex <= phase_parameters$end
      
      data$phase[index]   <- phase_parameters$phase
      
      data$epoch[index]   <- phase_parameters$epoch
      data$midline[index] <- phase_parameters$midline
      data$lcl[index]     <- phase_parameters$lcl
      data$ucl[index]     <- phase_parameters$ucl
      data[[paste0('DatePhase', phase_parameters$phase)]] <- phase_parameters$start

      # Only adjust series if 1) requested by user,
      # 2) epoch is 2 or 3,
      # 3) at least 21 days of records were observed in this phase.
      if (adjust && phase_parameters$epoch %in% c(2, 3) && sum(index) >= 21) {
        
        residuals <- log10(data$New_var) - log10(data$midline)

        data$weekday <- lubridate::wday(data$datex)
        
        data$residual_by_weekday[index] <- 
          ave(residuals[index],
              data$weekday[index],
              FUN = function(x) median(x, na.rm = TRUE))
        
        adjusted_var <- 10 ^ 
          (log10(data$New_var[index]) - data$residual_by_weekday[index])
        
        adjusted_var[is.na(data$New_var_Dump[index]) & (!is.finite(adjusted_var) | adjusted_var < 0)] <- 0
        
        # normalize to actual total var in phase
        adjusted_var <- adjusted_var * (sum(data$New_var[index], na.rm = TRUE) / sum(adjusted_var, na.rm = TRUE))
        
        data$New_var[index] <- round(adjusted_var)
        
        data$residual_by_weekday <- NULL
      }
    }
    
    if (adjust) {
      data$New_var[!is.finite(data$New_var) & is.na(data$New_var_Dump)] <- 0
      
      data <- find_phase_dates(data = data[c('GEO', 'NAME', 'datex', 'New_var', 'New_var_max', 'New_var_Dump', 'var_name')],
                               adjust = FALSE,
                               ghost = FALSE)
      
    }
    
    #the next conditions create the values for final Epoch, I don't understand why the country and state are different in the 
    #original code.   We don't need the values in the adaptation I am making here...
    #if (data$level[1] == 'state') {
      #data$EPOCH_last <- tail(data$epoch[!is.na(data$epoch)], 1)
      
      #data$EPOCH_last[is.infinite(data$EPOCH_last)] <- NA
   # } else if (data$level[1] == 'country') {
   #   data$EPOCH <- data$epoch
  #  }
    
    data <- data[data$datex >= phase_data[[1]]$start, ]
  })
  
  if ('try-error' %in% class(result)) {
    #browser()
    message(' -- Error encountered for ', data$var_name[1])
  }

  data
}
  

##############functions to read the JSON files from WI DPH website
JSON_fread_state <- function(state_url = state_url){
  list_JSON_state <- fromJSON(state_url, flatten=TRUE)
  
  df_JSON_state <- list_JSON_state$features
  
  names(df_JSON_state) <- gsub("properties.","",names(df_JSON_state))
  
  df_JSON_state$Date_time <- as.POSIXct(gsub("\\+00","",df_JSON_state$DATE))
  df_JSON_state$Date_reported <- as.Date(df_JSON_state$Date_time, tz = "US/Central")
  
  return(df_JSON_state)
}

JSON_fread <- function(state_url = state_url, county_url = county_url){
  list_JSON_county <- fromJSON(county_url, flatten=TRUE)
  
  df_JSON_county <- list_JSON_county$features
  
  list_JSON_state <- fromJSON(state_url, flatten=TRUE)
  
  df_JSON_state <- list_JSON_state$features[,c(1:67,109)]
  
  df_JSON <- rbind.data.frame(df_JSON_state,df_JSON_county)
  
  names(df_JSON) <- gsub("properties.","",names(df_JSON))
  
  df_JSON$Date_time <- as.POSIXct(gsub("\\+00","",df_JSON$DATE))
  df_JSON$Date_reported <- as.Date(df_JSON$Date_time, tz = "US/Central")
  
  return(df_JSON)
}


#function to create data frame for selected variables
make_df <- function(i, subset_list1 = subset_list, skip_start1 = skip_start, vars_common1 = vars_common){
  var_use <- subset_list1[i]
  
  df <- df1[,c(vars_common1,var_use)]
  
  names(df)[5] <- "var_cum"
  
  df_out <- df %>% 
    arrange(datex) %>% 
    filter(!is.na(var_cum)) %>% 
    mutate(New_var = force_monotonicity(c(var_cum[1],diff(var_cum)))) %>% 
    mutate(New_var_max = max(New_var, na.rm = TRUE), var_name = var_use) 
  
  if(var_use %in% skip_start1) {
    df_out <- df_out[-1,]
  } 
  
  return(df_out)
} 
