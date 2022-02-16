# Function to simulate N gamma mixtures with different cluster sizes and
# different non-spatial attributes
simulate_gamma_mixture_n <- function(N, n_total, k, ...){
  
  dat_N <- list()
  plot_N <- list()
  
  for (i in 1:N) {
    
      # Simulate vector of cluster sizes such that 
      # the total number of points is n_total
      n <- simulate_cluster_sizes(n_total, k)
      
      # Simulate one gamma mixture
      dat_N[[i]] <- simulate_gamma_mixture(n = n, k = k,...) 
      
      # Create non-spatial attributes
      dat_N[[i]] <- simulate_nonspat_attributes(dat_N[[i]])
  }
  
  return(dat_N)
}


# Simulate normal distribution with different parameters
# for each non-spatial group
simulate_group_normal <- function(x){
  
  # Parameters for each group 
  # (4 groups each with mu vector and Sigma matrix)
  params <- tribble(
    ~group, ~mu, ~Sigma,
    1, c(1, 5, 7), diag(c(1,2,2)),
    2, c(4, 5, 6), diag(c(1,0.5,1)),
    3, c(8, 4, 1), diag(c(2,2,2)),
    4, c(5, 7, 3), diag(c(1,1,2)),
  )
  
  # Generate random sample of 100 observation for each category
  values <- apply(as.matrix(x), 1, function(x) {
    switch(
      x,
      "1" =
        mvrnorm(1, params[1,]$mu[[1]], params[1,]$Sigma[[1]]),
      "2" =
        mvrnorm(1, params[2,]$mu[[1]], params[2,]$Sigma[[1]]),
      "3" =
        mvrnorm(1, params[3,]$mu[[1]], params[3,]$Sigma[[1]]),
      "4" =
        mvrnorm(1, params[4,]$mu[[1]], params[4,]$Sigma[[1]]),
      rep(9999, 3)
    )
  })
  
  
  #return(values)
  return(data.frame(par1 = values[1,],
                    par2 = values[2,],
                    par3 = values[3,]))
}

# Simulate multinomial distribution with different probability parameters
# for each non-spatial group
simulate_group_multinomial <- function(x){
  
  # Probabilities for each group (4 groups each with 3 parameters)
  probs <- tibble(
    group1 = c(0.05, 0.55, 0.4),
    group2 = c(0.33, 0.34, 0.33),
    group3 = c(0.15, 0.15, 0.7),
    group4 = c(0.40, 0.50, 0.10)
  )
  
  # Generate random sample of 100 observation for each category
  values <- apply(as.matrix(x), 1, function(x) {
    switch(
      x,
      "1" =
        rmultinom(1, 100, probs$group1) |> t(),
      "2" =
        rmultinom(1, 100, probs$group2) |> t(),
      "3" =
        rmultinom(1, 100, probs$group3) |> t(),
      "4" =
        rmultinom(1, 100, probs$group4) |> t(),
      rep(9999, 3)
    )
  })
  
  return(data.frame(par1 = values[1,],
                    par2 = values[2,],
                    par3 = values[3,]))
}

# Function to create non-spatial attribute groups for the data set and
# simulating n variables for each point
simulate_nonspat_attributes <- function(dat, nonspat_dist = "normal"){
  
  # Divide the data randomly to four sections
  div_x <- runif(1, 0.2, 0.8)
  div_y <- runif(1, 0.2, 0.8)
  
  dat$div <- c(x = div_x, y = div_y)
  
  dat$Y <- dat$Y |> 
    # Create a variable to describe the section
    dplyr::mutate(non_spat_group = if_else(x < div_x,
                                    true = if_else(y > div_y, "1", "3"),
                                    false = if_else(y > div_y, "2", "4")))
  
  
  if(nonspat_dist == "multinomial"){
    # Simulate variables from multinomial distribution
    dat$Y <- dat$Y |>
      mutate(simulate_group_multinomial(non_spat_group))
  } else if(nonspat_dist == "normal"){
    # Simulate variables from normal distribution
    dat$Y <- dat$Y |>
      mutate(simulate_group_normal(non_spat_group))
  } else {
    stop("Choose an available distribution for non-spatial attributes!")
  }
  
  
  
  
  return(dat)
}

# Function to simulate cluster sizes based on the total number of points
# and categorical distribution
simulate_cluster_sizes <- function(n_total, k, probs = rep(1/k, k)){
  
  tibble(x = rcat(n_total, probs)) |> 
    group_by(x) |> 
    dplyr::summarise(count = n()) |> 
    pull(count)
  
}

#' simulate_gamma_mixture
#'
#' Generate data points from multivariate gamma distributions.
#' @param n Number of data points. Can be a constant or a vector specifying the size of different clusters.
#' @param k Number of clusters.
#' @param w_dist_params Lower and upper limit of weights. Default is ]1, 100[.
#' @param w_dist Distribution of weights (not used). Default unifrom.
#' @param shape_params Shape parameters of Gamma distributions, vector of length 4. Default is runif(4, 0, 15) which are sorted in pairs.
#' @param rate_params rate parameters of Gamma distributions, vector of length 4. Default is runif(4, 0, 100) which are sorted in pairs.
#' @param n_out Number of outliers. If not specified, no outliers are generated.
#' @param out_scale How far outlier points should be from the most extreme values in SE sense. Default is 5.
#' @param outgroup_alpha How close to the grid edge outgroup members are sampled (a scale parameter between 0 and 1). Default is 0.5.
#' @param scale_between_range Scale the output data points between a user defined range.
#' @param place_on_grid Place clusters uniformly on a k x k grid. If "TRUE", cluster centers are moved to grid centers. Default is "FALSE".
#' @param overlap_scale Control the overlap of data points from different clusters. If greater than zero, points are moved towards the center of the corresponding cluster. Default is zero.
#' @return A list containing the following components:
#' \itemize{
#' \item Y - Table of data points (x, y), weights (w) and the original cluster label (orig_group).
#' \item mu_true - Original cluster centers.
#' }
#' 
#' @keywords cluster simulation rpack gamma distribution 
#' @export
#' @examples
#' seed = Sys.time()
#'
#' seed = as.integer(seed)
#'
#' seed = seed %% 100000
#'
#' seed = 10376
#'
#' set.seed(seed)
#'
#' k = 10
#'
#' n = 500
#'
#' n_out = 30
#'
#' test_data = simulate_gamma_mixture(n, k, n_out=n_out, scale_between_range = c(0, 1), outgroup_alpha = 0.9, place_on_grid = T, 
#'                                  overlap_scale = 0.4)
#'
#' true_mu = test_data$mu_true
#'
#' test_data = test_data$Y
#'
#' plot(test_data$x, test_data$y, col=c(rep(rainbow(k), each=(round(n/k))), rep("black", n_out)), 
#'     pch=rep(16, n + n_out), cex=test_data$w/90)
#'
#' points(true_mu[ ,1], true_mu[ ,2], pch=4, cex=1.2, lwd=3)
#'
#' @export
#'
simulate_gamma_mixture <- function(n, k, w_dist_params = c(1, 100), w_dist = "uniform", shape_params = NULL,
                                   rate_params = NULL, n_out=NULL, out_scale=5, outgroup_alpha=0.5, scale_between_range=NULL, 
                                   place_on_grid=F, overlap_scale=0){
  
  if(is.null(shape_params)){
    
    shape_params <- runif(4, 0, 15)
    
    shape_params <- sort(shape_params)
  
    shape_params <- shape_params[c(1, 4, 2, 3)]
      
  }
  
  if(is.null(rate_params)){
    
    rate_params <- runif(4, 0, 100)
    
    rate_params <- sort(rate_params)
    
    rate_params <- rate_params[c(1, 4, 2, 3)]
    
  }
  
  Y <- NULL
  
  mu <- matrix(0, k, 2)
  
  on_edge <- sample(c(T, F), k, replace = T)
  
  for(i in 1:k){
    
    K <<- i
    
	if(length(n) > 1){
      
      n_sub <- n[i]
      
      xy_dist <- matrix(0, n_sub, ncol=k)
      
    }
    
    if(length(n) == 1){
      
      n_sub <- round(n/k)
      
      xy_dist <- matrix(0, round(sum(n)/k), ncol=k)
      
    }
    
    s <- c(runif(1, shape_params[1], shape_params[2]), runif(1, shape_params[3], shape_params[4]))
    
    r <- 1/c(runif(1, rate_params[1], rate_params[2]), runif(1, rate_params[3], rate_params[4]))
    
    move_dist <- runif(2, -10, 10)
    
    mu[K, ] <- s/r + move_dist
    
    Sigma <- rpack::random_sigma()
    
    xy <- rmvgamma(n_sub, s, r, Sigma)
    
    xy <- sweep(xy, 2, move_dist, "+")
    
    orig_group <- rep(i, n_sub)
    
    w <- floor(runif(n_sub, min=w_dist_params[1], max=w_dist_params[2]) + 1)
    
    Y <- rbind(Y, data.frame(x = xy[ , 1], y = xy[ , 2], w=w, orig_group = as.factor(orig_group)))
    
    # The weight depends on the Euclidean distance between the point and cluster center (gamma distribution mean) 
    
    if(on_edge[i] == T){
      
      w <- sort(w)
      
      #print(w)
      Y$w[Y$orig_group == i] <- w
      
      e_dist <- function(x) sqrt(sum((x - mu[i, ])^2))
      
      xy_dist[ , i] <- apply(Y[Y$orig_group == i, c("x", "y")], 1, e_dist)
      
      xy_dist[ , i] <- xy_dist[ , i]/sum(xy_dist[ , i])
      
      sort_dist <- sort(xy_dist[ , i])
      
      d <- Y$w[Y$orig_group == i]
      
      Y$w[Y$orig_group == i] <- d[match(xy_dist[ , i], sort_dist)] 
      
    }
    
  }
  
  if(place_on_grid == T){
    
    a <- k
    
    b <- k
    
    x_grid <- seq(min(Y$x), max(Y$x), length.out = a + 1)
    y_grid <- seq(min(Y$y), max(Y$y), length.out = b + 1)
    
    Midx <- (x_grid[-1] + head(x_grid, -1)) /2
    Midy <- (y_grid[-1] + head(y_grid, -1)) /2
    
    g = expand.grid(Midx, Midy)
    
    colnames(g) = c("x", "y")
    
    grid_points = sample(1:nrow(g), k)
    
    for(i in 1:k){
      
      Y$x[Y$orig_group == i] <- Y$x[Y$orig_group == i] - mean(Y$x[Y$orig_group == i]) + g[grid_points[i], "x"]
      
      Y$y[Y$orig_group == i] <- Y$y[Y$orig_group == i] - mean(Y$y[Y$orig_group == i]) + g[grid_points[i], "y"]

    }
    
    mu <- g[grid_points, ]
    
  }
  
  if(overlap_scale > 0){
    
    for(i in 1:k){
      
      d <- sqrt((Y$x[Y$orig_group == i] - mu[i , 1])^2 + (Y$y[Y$orig_group == i] - mu[i , 2])^2)
      
      kos <- abs(Y$x[Y$orig_group == i] - mu[i , 1])/d
      
      zin <- abs(Y$y[Y$orig_group == i] - mu[i , 2])/d
      
      Y$x[Y$orig_group == i] <- Y$x[Y$orig_group == i] + (Y$x[Y$orig_group == i] < mu[i , 1])*kos*d*overlap_scale - 
        (Y$x[Y$orig_group == i] > mu[i , 1])*kos*d*overlap_scale
      
      Y$y[Y$orig_group == i] <- Y$y[Y$orig_group == i] + (Y$y[Y$orig_group == i] < mu[i , 2])*zin*d*overlap_scale - 
        (Y$y[Y$orig_group == i] > mu[i , 2])*zin*d*overlap_scale
      
    }
    
  }
  
  if(!is.null(n_out)){
    
    sim_out_data <- simulate_unif_grid(n_out = n_out, alpha = outgroup_alpha, x = Y$x, y = Y$y)
    
    w <- floor(runif(n_out, min=w_dist_params[1], max=w_dist_params[2]) + 1)
    
    X <- data.frame(x = sim_out_data$x, y = sim_out_data$y, w = w, orig_group = as.factor(k+1))
    
    Y <- rbind(Y, X)
    
  }

  if(!is.null(scale_between_range)){
    
    mu[ , 1] <- (scale_between_range[2] - scale_between_range[1])*((mu[ , 1] - min(Y$x))/
                                                                    (max(Y$x) - min(Y$x))) + scale_between_range[1]
     
    mu[ , 2] <- (scale_between_range[2] - scale_between_range[1])*((mu[ , 2] - min(Y$y))/
                                                                    (max(Y$y) - min(Y$y))) + scale_between_range[1]
     
    Y$x <- (scale_between_range[2] - scale_between_range[1])*((Y$x - min(Y$x))/(max(Y$x) - min(Y$x))) + scale_between_range[1]
     
    Y$y <- (scale_between_range[2] - scale_between_range[1])*((Y$y - min(Y$y))/(max(Y$y) - min(Y$y))) + scale_between_range[1]
     
  }  
  
  # Add column to identify outliers easily
  Y <- Y |> 
    mutate(is_outlier = if_else(orig_group == k + 1, TRUE, FALSE))
  
  return(list(Y = Y, mu_true = mu))
  
}