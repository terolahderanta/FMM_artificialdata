duration_sd_summary <- function(clust_list, names, dat){
  
  k <- clust_list[[1]]$cluster %>% unique() %>% length()
  
  duration_sds <- NULL
  
  for (i in 1:length(clust_list)) {
    duration_sds <- rbind(duration_sds,
                          cluster_duration_sd(clust_list[[i]], dat, k))
  }
  
  duration_summary <- apply(duration_sds, 1, FUN = summary) %>% t()
  
  names <- row.names(duration_sds)
  
  
  duration_summary <- 
    duration_summary %>% 
    as.data.frame() %>% 
    mutate(name = names, .before = 1) %>% 
    as_tibble()
  
  #names(duration_summary) <- names
  # 
  # duration_sds <- 
  #   duration_sds %>%
  #   as.data.frame() %>% 
  #   mutate(name = names, .before = 1)
  # 
  # names(duration_sds)[2:(k+1)] <- c(paste("k", 1:k, sep = ""))
  # 
  # return(duration_sds)
  
  return(duration_summary)
}

# Calculate sd of the durations in each cluster
# and return a vector of sds
cluster_duration_sd <- function(clust_obj, dat, k){
  
  ret_mat <- sapply(X = 1:k, FUN = function(X){
    dat %>% 
      filter(clust_obj$clusters == X) %>% 
      select(par1,par2,par3) %>% 
      summarise(across(.fns = sd)) |> 
      as.matrix()
  }
  )
  
  row.names(ret_mat) <- paste("par", 1:3,"_lambda",clust_obj$lambda_par,sep ="")
  
  return(ret_mat)
}

table_proximity_summary <- function(data, cl_obj, dig = 2, probs = c(0, 0.25, 0.5, 0.75, 1), frac = FALSE){
  
  coords <- data %>% select(x,y) %>% as.matrix()
  
  weights <- data %>% pull(w)
  
  wd <- weights_and_dist(coords = coords,
                         weights = weights,
                         cl_obj = cl_obj)
  
  sizes <- cl_sizes(weights = weights,
                    cl_obj = cl_obj,
                    frac = frac)
  
  res <- c(round(c(weighted.mean(wd$d, wd$w)*1, 
                   reldist::wtd.quantile(wd$d, q = probs, weight = wd$w)*1,
                   round(sd(sizes),digits = dig - 1),
                   min(sizes),
                   max(sizes)
  ), 
  digits = dig))
  
  n_res <- length(res)
  
  names(res)[1] <- "Mean"
  names(res)[n_res - 2] <- "S.D."
  names(res)[n_res - 1] <- "Min"
  names(res)[n_res] <- "Max"
  return(res)
}

# Create a proximity summary table for clusterings
proximity_summary <- function(clust_list, names, dat, probs = c(0, 0.25, 0.5, 0.75, 1)){
  summ_table <- table_proximity_summary(data = dat,
                                        cl_obj = clust_list[[1]],
                                        probs = probs)
  
  if(length(clust_list) <= 1) return(summ_table)
  
  for (i in 2:length(clust_list)) {
    summ_table <- rbind(summ_table,
                        table_proximity_summary(data = dat,
                                                cl_obj = clust_list[[i]],
                                                probs = probs))
  }
  
  summ_table <- 
    summ_table %>% 
    as.data.frame() %>% 
    mutate(name = names, .before = 1) %>% 
    as_tibble()
  
  return(summ_table)
}


# Returns a summary data table of session durations on each AP
get_durations_summary <- function(dat){
  
  # Minimum and maximum of durations
  #dur_min <- dat %>% pull(duration) %>% min()
  #dur_max <- dat %>% pull(duration) %>% max()
  
  # divide the durations to 10 partition groups
  #n_part <- 10
  #width_part <- (dur_max - dur_min)/n_part
  #bound_part <- (0:n_part) * width_part
  #partitions <- tibble(L = bound_part[-(length(bound_part))],
  #                     U = bound_part[-1])
  
  
  ret_dat <- dat %>% 
    # Group by AP
    group_by(ap_id) %>% 
    
    # Get summary statistics
    mutate(mean_duration = mean(duration), 
           sd_duration = sd(duration),
           q_10 = quantile(duration, probs = 0.1),
           q_90 = quantile(duration, probs = 0.9),
           n = n()) %>%
    
    # Return only summary information
    select(-duration, - start_time, -end_time) %>% 
    
    # Only one row per AP
    distinct() %>% 
    
    # Order data with lowest mean to highest
    arrange(mean_duration) %>% 
    
    # Ungroup the data
    ungroup() %>%
    
    # Focus on the center: xlim = c(120.9, 121.9), ylim = c(30.6,31.7)
    filter(between(x, 120.9, 121.9)) %>%
    filter(between(y, 30.6, 31.7))
  
  # Count the maximum number of simultaneous connections for each AP
  max_simult <- sapply(X = ret_dat %>% pull(ap_id), FUN = peak_in_ap, dat = dat)
  
  
  ret_dat <-  ret_dat %>% 
    
    # Add those to the data frame
    add_column(max_simult = max_simult)
  
  return(ret_dat)
}

#' Calculates distances from ap to its center and weight assigned to that
weights_and_dist <- function(coords, weights, cl_obj){
  k <- nrow(cl_obj$centers)
  
  points <- which(cl_obj$clusters != 99)
  #n <- nrow(coords)
  
  centers <- cl_obj$centers %>% as_tibble() %>% select(x,y) %>% as.matrix()
  
  
  # d <- matrix(0, ncol = k, nrow = n)
  # w <- matrix(0, ncol = k, nrow = n)
  
  d <- matrix(0, ncol = k, nrow = length(points))
  w <- matrix(0, ncol = k, nrow = length(points))
  
  for (i in 1:length(points)) {
    d[i,] <- apply(centers, MARGIN = 1, FUN =  function(x) {
      pracma::haversine(x, coords[points[i],])
    })  
    w[i,] <- cl_obj$assign_frac[points[i],1:k] * weights[points[i]]  
  }
  
  return(list(d = d, w = w))
}

# Sizes of the clusters
cl_sizes <- function(weights, cl_obj, frac = FALSE){
  k <- nrow(cl_obj$centers)
  
  if(!frac){
    res <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) {sum(weights[cl_obj$cluster == x])}
    )
  } else {
    res <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) { sum(weights*cl_obj$assign_frac[,x]) }
    )
  }
  return(res)
}