clust_with_params_list <- function(dat_list, lambda_par, lambda,...){
  n_list <- length(dat_list)
  ret_list <- list()
  
  for (i in 1:n_list) {
    ret_list <- append(
      ret_list,
      list(clust_with_params(
        coords = dat_list[[i]]$Y |> dplyr::select(x, y),
        params = dat_list[[i]]$Y %>% dplyr::select(par1,par2,par3),
        weights = dat_list[[i]]$Y %>% pull(w),
        lambda_par = lambda_par,
        lambda = lambda,
        ...
      ))
    )
  }
  
  return(ret_list)
}


# Cluster with params and different lambdas, return a list of clusterings
clust_with_params <- function(coords, 
                              params, 
                              weights,
                              lambda_par, 
                              lambda, ...) {
  clust_list <- list()
  
  for (lambda_pari in lambda_par) {
    clust_list$lambda_par <- lambda_pari
    
    clust_list <- append(clust_list,
                         list(alt_alg(coords = coords, 
                                      params = params,
                                      weights = weights, 
                                      lambda_params = lambda_pari,
                                      lambda = lambda, ...)))
  }
  
  
  
  return(clust_list)
}

#' Different colors for the plots
#' @export
c_col <- c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
           "chocolate","darkred","yellow3","darkgreen","bisque4","magenta",
           "royalblue","tomato4","steelblue1",
           "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
           "coral3","beige","brown4","indianred1","lightgreen","orchid")

#' Calculate the mode of vector v
#' @param v A vector.
#' @return The mode.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Calculate the difference between two points
#' @param x1 1. point. (coordinates first)
#' @param x2 2. point. (coordinates first)
#' @param lambda Relative weight of coordinates vs. params.  
#' @return Squared Euclidean distance.
#' @export
euc_dist2_par <- function(x1, x2){
  n_coords <- 2
  x1_coord <- x1[1:n_coords]
  x2_coord <- x2[1:n_coords]
  x1_par <- x1[-(1:n_coords)]
  x2_par <- x2[-(1:n_coords)]
  
  lambda <- 0.001
  
  # Coordinate part
  lambda * sum((x1_coord - x2_coord) ^ 2) +
    
    # Parameter part
    (1 - lambda) * sum((x1_par - x2_par) ^ 2)
}

#' Calculate the difference between two points
#' @param x1 1. point. (coordinates first)
#' @param x2 2. point. (coordinates first)
#' @param lambda Relative weight of coordinates vs. params.  
#' @return Squared Euclidean distance.
#' @export
euc_dist3_par <- function(x1, x2){
  n_coords <- 2
  x1_coord <- x1[1:n_coords]
  x2_coord <- x2[1:n_coords]
  x1_par <- x1[-(1:n_coords)]
  x2_par <- x2[-(1:n_coords)]
  
  lambda <- 0.001
  
  # Coordinate part
  lambda * sum(abs(x1_coord - x2_coord) ^ 3) +
    
    # Parameter part
    (1 - lambda) * sum((x1_par - x2_par) ^ 2)
}

#' Calculate the difference between two points
#' @param x1 1. point. (coordinates first)
#' @param x2 2. point. (coordinates first)
#' @return Squared Euclidean distance.
#' @export
scaled_euc2_dist <- function(x1, x2, sigma = 1){
  n_coords <- 2
  
  # spatial coordinates
  x1_coord <- x1[1:n_coords]
  x2_coord <- x2[1:n_coords]
  
  # data geometry
  x1_par <- x1[-(1:n_coords)]
  x2_par <- x2[-(1:n_coords)]
  
  # the more data values differ from each other, the further the spatial points are from each other
  sum((x1_coord - x2_coord) ^ 2 / gaussian_scaler(x1_par, x2_par, sigma)) 
}


#' Calculate the difference between two points
#' @param x1 1. point. (coordinates first)
#' @param x2 2. point. (coordinates first)
#' @return Squared Euclidean distance.
#' @export
scaled_euc3_dist <- function(x1, x2, sigma = 0.2){
  n_coords <- 2
  
  # spatial coordinates
  x1_coord <- x1[1:n_coords]
  x2_coord <- x2[1:n_coords]
  
  # data geometry
  x1_par <- x1[-(1:n_coords)]
  x2_par <- x2[-(1:n_coords)]
  
  # the more data values differ from each other, the further the spatial points are from each other
  sum(abs(x1_coord - x2_coord) ^ 3 / gaussian_scaler(x1_par, x2_par, sigma)) 
}


#' Calculate the Euclidean distance between two points
#' @param x1 1. point.
#' @param x2 2. point.
#' @return Euclidean distance.
#' @export
euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#' Calculate the squared Euclidean distance between two points
#' @param x1 1. point.
#' @param x2 2. point.
#' @return Squared Euclidean distance.
#' @export
euc_dist2 <- function(x1, x2) sum((x1 - x2) ^ 2)

#' Calculate the medoid of the data points
#' @param data A data.frame.
#' @param w Weights of the data points.
#' @param d A distance metric.
#' @export
#' @return The medoid.
medoid <- function(data,
                   w = rep(1, nrow(data)),
                   d = euc_dist2) {
  n <- nrow(data)
  if (n < 1) {
    stop("Tried to calculate medoid from zero number of points! (rpack)")
  }
  if (n == 1) {
    return(data[1, ])
  }
  
  w_dists <- sapply(
    1:n,
    FUN = function(x) {
      sum(w[-x] * apply(data[-x, ], 1, FUN = d, x2 = data[x, ]))
    }
  )
  
  return(data[which.min(w_dists), ])
}

#' Calculate the medoid from distance matrix
#' @param dist_mat Distance matrix for the data points.
#' @param ids Ids for the points in distance matrix. Uses all of the points by default.  
#' @param w Weights of the data points.
#' @export
#' @return The id for the medoid.
medoid_dist_mat <- function(dist_mat,
                            ids = 1:nrow(dist_mat),
                            w = rep(1, nrow(dist_mat))) {
  
  # Exceptions
  n <- nrow(dist_mat)
  if (n < 1 | length(ids) == 0) {
    stop("Tried to calculate medoid from zero number of points! (rpack)")
  }
  if (n == 1 | length(ids) == 1) {
    return(ids[1])
  }
  
  # Weighted distances from the given set of points
  wdists <- dist_mat[ids, ids] * w[ids]
  
  # Calculate column sums
  wdist_to_centers <- colSums(wdists)
  
  return(ids[which.min(wdist_to_centers)])
}

#' Calculate the medoid from distance matrix
#' @param dist_mat Distance matrix for the data points.
#' @param ids Ids for the points in distance matrix. Uses all of the points by default.  
#' @param w Weights of the data points.
#' @export
#' @return The id for the medoid.
medoid_dist_to_centers <- function(dist_to_centers,
                                   ids = 1:nrow(dist_to_centers),
                                   w = rep(1, nrow(dist_to_centers))) {
  
  # Exceptions
  n <- nrow(dist_to_centers)
  if (n < 1 | length(ids) == 0) {
    stop("Tried to calculate medoid from zero number of points! (rpack)")
  }
  
  #cat("length of w-vector:")
  #print(length(w[ids]))
  #cat(paste("dist_to_center dimension in medoid:"))
  #print(dim(dist_to_centers[ids,]))
  
  
  # Weighted distances from the given set of points
  wdists <- dist_to_centers[ids, ] * w[ids]
  
  # Calculate column sums
  if(length(ids) == 1){
    wdist_to_centers <- wdists
  } else {
    wdist_to_centers <- colSums(wdists)
  }
  return(which.min(wdist_to_centers))
}

#' Kmeans++
#'
#' Implementation of the K-means++ algorithm. Whereas normal kmeans selects all the initial center
#' cluster centers randomly, kmeans++ randomly selects only the first center. For each
#' consecutive center, the probability of selection is weighed by the distance to already selected
#' centers.
#'
#' Implementation adapted from one by Hans Werner (https://stat.ethz.ch/pipermail/r-help/2012-January/300051.html).
#'
#' See following article for more information on kmeans++:
#'
#'   Arthur, D., and Vassilvitskii, S. "k-means++: The advantages of careful seeding."
#'   Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms. Society
#'   for Industrial and Applied Mathematics, 2007.
#'
#' @param X A matrix or a data frame containing the objects, one per row.
#' @param k Number of clusters.
#' @export
kmpp <- function(X, k) {
  
  if (!is.matrix(X)) X <- as.matrix(X)  # X must be a matrix
  
  n <- nrow(X)
  ncoords <- ncol(X)
  C <- numeric(k)                # initialize centers to zero
  C[1] <- sample(1:n, size = 1)  # select first element randomly
  
  for (i in 2:k) {
    dm <- pracma::distmat(X, matrix(X[C, ], ncol = ncoords))
    pr <- apply(dm, 1, min)
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  if(length(unique(C)) == k){
    cl <- stats::kmeans(X, X[C, ])  
    cl$initial_centers <- X[C,]
  } else {
    cl <- kmpp(X, k)
  }
  #cl <- stats::kmeans(X, X[C, ])
  #cl$initial_centers <- X[C,]
  
  return(cl)
}



#' Full alternating algorithm
#'
#' @param coords Coordinates of the data points.
#' @param weights Weights of the points in a vector.
#' @param k Number of clusters.
#' @param N Number of starting values.
#' @param range Limits for the cluster size in a list.
#' @param capacity_weights Different weights for capacity limits.
#' @param d Distance function used in clustering.
#' @param center_init Options to initialize center locations. Default is "random" and other choice is "kmpp". 
#' @param lambda Outgroup parameter.
#' @param frac_memb Can points be partially allocated?
#' @param place_to_point Place the cluster head in a point?
#' @param fixed_centers Possible fixed center locations.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to.
#' @param dist_mat Distance matrix for all the points. 
#' @param print_output Different types of printing outputs, "progress" is default and "steps" stepwise-print.
#'
#' @return Clustering object with allocation, center locations and the value of the objective function
#' @export
alt_alg <- function(coords, 
                    weights, 
                    k,
                    params = NULL,
                    N = 10, 
                    range = c(min(weights)/2, sum(weights)),
                    capacity_weights = weights, 
                    d = euc_dist2,
                    center_init = "random", 
                    lambda = NULL,
                    frac_memb = FALSE, 
                    place_to_point = TRUE, 
                    fixed_centers = NULL, 
                    gurobi_params = NULL,
                    multip_centers = rep(1, nrow(coords)),
                    dist_mat = NULL,
                    print_output = "progress",
                    normalization = TRUE,
                    lambda_fixed = NULL,
                    lambda_params = NULL){
  
  # Check arguments
  assertthat::assert_that(is.matrix(coords) || is.data.frame(coords), msg = "coords must be a matrix or a data.frame!")
  
  assertthat::assert_that(nrow(coords) >= k, msg = "must have at least k coords points!")
  assertthat::assert_that(is.numeric(weights), msg = "weights must be an numeric vector!")
  assertthat::assert_that(is.numeric(capacity_weights), msg = "capacity weights must be an numeric vector!")
  assertthat::assert_that(length(weights) == nrow(coords), msg = "coords and weight must have the same number of rows!")
  assertthat::assert_that(length(capacity_weights) == nrow(coords), msg = "coords and capacity weights must have the same number of rows!")
  assertthat::assert_that(is.numeric(k), msg = "k must be a numeric scalar!")
  assertthat::assert_that(length(k) == 1, msg = "k must be a numeric scalar!")
  
  assertthat::assert_that(is.numeric(range))
  
  if(!purrr::is_null(lambda)) assertthat::is.number(lambda)
  if(!purrr::is_null(lambda_fixed)) assertthat::is.number(lambda_fixed)
  
  assertthat::assert_that(is.logical(normalization), msg = "normalization must be TRUE or FALSE!")
  assertthat::assert_that(is.logical(frac_memb), msg = "frac_memb must be TRUE or FALSE!")
  assertthat::assert_that(is.logical(place_to_point), msg = "place_to_point must be TRUE or FALSE!")
  
  # Calculate distance matrix
  if(is.null(dist_mat) & place_to_point){
    
    # Print information about the distance matrix
    n <- nrow(coords)
    cat(paste("Creating ", n, "x", n ," distance matrix... ", sep = ""))
    temp_mat_time <- Sys.time()
    
    if(is.null(params)) {
      # Calculate distances with distance metric d
      dist_mat <- apply(
        X = coords,
        MARGIN = 1,
        FUN = function(x) {
          apply(
            X = coords,
            MARGIN = 1,
            FUN = d,
            x2 = x
          )
        }
      )
    } else {
      # Calculate distances with distance metric d
      # dist_mat <- apply(
      #   X = cbind(coords, params),
      #   MARGIN = 1,
      #   FUN = function(x) {
      #     apply(
      #       X = cbind(coords, params),
      #       MARGIN = 1,
      #       FUN = d,
      #       x2 = x, 
      #       lambda = lambda_params
      #     )
      #   }
      # )
      
      # Calculate distance matrix based on spatial distance
      dist_mat_coords <- apply(
        X = coords,
        MARGIN = 1,
        FUN = function(x) {
          apply(
            X = coords,
            MARGIN = 1,
            FUN = hav.dist2,
            x2 = x
          )
        }
      ) |> scale()
      
      # Calculate distance matrix based on parameter distance
      dist_mat_params <- apply(
        X = params,
        MARGIN = 1,
        FUN = function(x) {
          apply(
            X = params,
            MARGIN = 1,
            FUN = rpack::euc_dist,
            x2 = x
          )
        }
      ) |> scale()
      
      # Combine the two distance matrix with lambda
      dist_mat <- lambda_params * dist_mat_coords +
        (1 - lambda_params) * dist_mat_params
      
    }
    
    cat(paste("Matrix created! (", format(round(Sys.time() - temp_mat_time)) ,")\n\n", sep = ""))
    
    # Normalizing distances
    if(normalization){
      dist_mat <- dist_mat / max(dist_mat)
    }
    
  } else if(place_to_point){
    
    # Normalizing distances
    if(normalization){
      dist_mat <- dist_mat / max(dist_mat)
    }
    
  } else {
    # If no distance matrix is used
    dist_mat <- NULL
  }
  
  if(normalization) {
    # Normalization for the capacity weights
    max_cap_w <- max(capacity_weights)
    capacity_weights <- capacity_weights / max_cap_w
    range <- range / max_cap_w
    
    # Normalization for the weights
    weights <- weights / max(weights)
  }
  
  # Print the information about run
  if(print_output == "progress"){
    cat(paste("Progress (N = ", N,"):\n", sep = ""))
    cat(paste("______________________________\n"))
    progress_bar <- 0
  } 
  
  # Total iteration time
  temp_total_time <- Sys.time()
  
  for (i in 1:N) {
    
    if(print_output == "steps"){
      cat(paste("\nIteration ", i, "/", N, "\n---------------------------\n", sep = ""))
      temp_iter_time <- Sys.time()
    }
    
    # One clustering
    temp_clust <- capacitated_LA(coords = coords,
                                 weights = weights,
                                 k = k,
                                 params = params,
                                 ranges = range,
                                 capacity_weights = capacity_weights,
                                 lambda = lambda,
                                 d = d,
                                 dist_mat = dist_mat,
                                 center_init = center_init,
                                 place_to_point = place_to_point,
                                 frac_memb = frac_memb,
                                 fixed_centers = fixed_centers,
                                 gurobi_params = gurobi_params,
                                 multip_centers = multip_centers,
                                 print_output = print_output,
                                 lambda_fixed = lambda_fixed)
    
    # Save the first iteration as the best one
    if(i == 1){
      min_obj <- temp_clust$obj
      best_clust <- temp_clust
    }
    
    # Print the number of completed laps
    if(print_output == "progress") {
      if((floor((i / N) * 30) > progress_bar)) {
        cat(paste0(rep("#", floor((
          i / N
        ) * 30) - progress_bar), collapse = ""))
        progress_bar <- floor((i / N) * 30)
      }
    } else if(print_output == "steps"){
      cat(paste("Iteration time: ", format(round(Sys.time() - temp_iter_time)), "\n", sep = ""))
    }
    
    # Save the iteration with the lowest value of objective function
    if(temp_clust$obj < min_obj){
      min_obj <- temp_clust$obj
      best_clust <-  temp_clust
    }
  }
  
  cat("\n\n")
  
  # Print total iteration time
  cat(paste("Total iteration time: ", format(round(Sys.time() - temp_total_time)),"\n", sep = ""))
  
  best_clust$lambda_par <- lambda_params
  best_clust$lambda <- lambda
  
  return(best_clust)
}

#' Capacitated location allocation  
#' 
#' Alternating algorithm for maximizing the optimization function in LA/clustering. Includes various constraints and extensions 
#' such as: Capacity limits, outliers, weights for the data points, different distance metrics, different
#' memberships, fixed points etc. 
#' 
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param ranges Lower and upper limits for the clustering
#' @param capacity_weights Different weights for capacity limits.
#' @param d The distance function.
#' @param center_init Options to initialize center locations. Default is "random" and other choice is "kmpp". 
#' @param fixed_centers Predetermined center locations.
#' @param lambda Outlier-parameter.
#' @param place_to_point if TRUE, cluster centers will be placed to a point.
#' @param frac_memb If TRUE memberships are fractional.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param dist_mat Distance matrix for all the points. 
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to. 
#' @param print_output Print options, default is NULL. One option is "steps" which prints information about steps.
#' @return A list containing cluster allocation, cluster center and the current value of the objective function.
#' @keywords internal
capacitated_LA <- function(coords,
                           weights,
                           k,
                           ranges,
                           params = NULL,  
                           capacity_weights = weights,
                           d = euc_dist2, 
                           center_init = NULL,
                           fixed_centers = NULL , 
                           lambda = NULL, 
                           place_to_point = TRUE, 
                           frac_memb = FALSE, 
                           gurobi_params = NULL, 
                           dist_mat = NULL,
                           multip_centers = rep(1, nrow(coords)),
                           print_output = NULL,
                           lambda_fixed = NULL){
  
  # Number of objects in data
  n <- nrow(coords)
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_centers), 0, nrow(fixed_centers))
  
  # Cluster centers with fixed centers first
  if(n_fixed > 0){
    
    if(place_to_point){
      
      fixed_center_ids <- prodlim::row.match(fixed_centers %>% as.data.frame(), coords %>% as.data.frame())
      #fixed_center_ids <- which(
      #  ((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))) & 
      #    ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2)))
      #)
      
    } 
  } else {
    fixed_center_ids <- NULL
  }
  
  # Initialize centers
  switch(center_init,
         # Pick initial centers randomly
         random = {
           
           if(place_to_point){
             # Don't sample the points in fixed_centers
             sample_points <- which(!(1:n %in% fixed_center_ids))
             center_ids <- c(fixed_center_ids, sample(sample_points, k - n_fixed))
             centers <- coords[center_ids,]
           } else {
             center_ids <- NULL
             centers <- coords[sample(1:n, k),]
           }
         },
         
         # Pick initial centers with k-means++
         kmpp = {
           
           # Replicate data points according to their weight
           weighted_coords <- apply(coords, 2, function(x) rep(x, round(weights)))
           
           # Run k-means++
           init_kmpp <- kmpp(weighted_coords, k)
           
           if(place_to_point){
             
             # Select closest points to the kmpp-result as the centers
             center_ids <- apply(X = init_kmpp$centers, MARGIN = 1, FUN = function(x){
               temp_dist_kmpp <- 
                 apply(X = coords, MARGIN = 1, FUN = function(y){d(x,y)})
               which.min(temp_dist_kmpp)
             })
             
             # Do not choose the same point twice
             ifelse(duplicated(center_ids),
                    sample((1:n)[-center_ids],size = 1),
                    center_ids)
             
             # Select the centers according to ids
             centers <- coords[center_ids,]
             
           } else {
             centers <- init_kmpp$centers
             center_ids <- NULL
           }
         },
         
         stop("No such choice for center initialization! (rpack)")
  )
  
  # Maximum number of laps
  max_laps <- 100
  
  for (iter in 1:max_laps) {
    # Old mu is saved to check for convergence
    old_centers <- centers
    
    # Print detailed steps
    if(print_output == "steps"){
      cat("A-step... ")
      temp_time <- Sys.time()
    }
    
    # Clusters in equally weighted data (Allocation-step)
    temp_alloc <- allocation_step(
      coords = coords,
      weights = weights,
      params = params,
      k = k,
      centers = centers,
      ranges = ranges,
      center_ids = center_ids,
      capacity_weights = capacity_weights,
      lambda = lambda,
      d = d,
      frac_memb = frac_memb,
      gurobi_params = gurobi_params,
      dist_mat = dist_mat,
      multip_centers = multip_centers
    )
    
    # Print detailed steps
    if(print_output == "steps"){
      cat(paste("Done! (", format(round(Sys.time() - temp_time, digits = 3), nsmall = 3), ")\n", sep = ""))
    }
    
    # Save the value of the objective function
    obj_min <- temp_alloc$obj_min
    
    # Print detailed steps
    if(print_output == "steps"){
      cat("L-step... ")
      temp_time <- Sys.time()
    }
    
    # Updating cluster centers (Parameter-step)
    temp_loc <- location_step(
      coords = coords,
      weights = weights,
      k = k,
      params = params,
      assign_frac = temp_alloc$assign_frac,
      fixed_centers = fixed_centers,
      d = d,
      place_to_point = place_to_point,
      dist_mat = dist_mat,
      lambda_fixed = lambda_fixed
    )
    
    # Print detailed steps
    if(print_output == "steps"){
      cat(paste("Done! (", format(round(Sys.time() - temp_time, digits = 3), nsmall = 3), ")\n", sep = ""))
    }
    
    centers <- temp_loc$centers
    
    center_ids <- temp_loc$center_ids
    
    # Print a message showing that max number of iterations was reached
    if(iter == max_laps & print_output == "steps"){
      cat(paste("WARNING! Reached maximum number of LA-iterations! Returning the clustering from last lap...\n",sep = ""))
    }
    
    # If nothing is changing, stop
    if(all(old_centers == centers)) break
  }
  
  # Save the assignments
  assign_frac <- temp_alloc$assign_frac
  
  # Hard clusters from assign_frac
  clusters <- apply(assign_frac, 1, which.max)
  
  # Cluster 99 is the outgroup
  clusters <- ifelse(clusters == (k+1), 99, clusters)
  
  # Return cluster allocation, cluster center and the current value of the objective function
  return(list(clusters = clusters, centers = centers, obj = obj_min, assign_frac = assign_frac))
}

#' Update cluster allocations by minimizing the objective function.
#'
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param centers The parameters (locations) that define the k distributions.
#' @param ranges Lower and upper limits for the clustering
#' @param center_ids Ids for the data points that are selected as centers.
#' @param capacity_weights Different weights for capacity limits.
#' @param lambda Outlier-parameter
#' @param d Distance function.
#' @param frac_memb If TRUE memberships are fractional.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param dist_mat Distance matrix for all the points.
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to. 
#' @return New cluster allocations for each object in data and the maximum of the objective function.
#' @keywords internal
allocation_step <- function(coords,
                            weights,
                            k,
                            centers,
                            ranges,
                            params = NULL,
                            center_ids = NULL,
                            capacity_weights = weights,
                            lambda = NULL,
                            d = euc_dist2,
                            frac_memb = FALSE,
                            gurobi_params = NULL,
                            dist_mat = NULL,
                            multip_centers = rep(1, nrow(coords))) {
  
  # Number of objects in data
  n <- nrow(coords)
  
  # Is there an outgroup cluster
  is_outgroup <- !is.null(lambda)
  
  # Number of range groups
  if(is.vector(ranges)){
    ranges <- matrix(data = ranges, nrow = 1, ncol = 2)
    g <- 1
  } else {
    g <- nrow(ranges) 
  }
  
  # Number of decision variables
  n_decision <- n * k + 
    # More than one range groups
    ifelse(g > 1, k * g, 0) +
    # Outliers
    ifelse(is_outgroup, n, 0)
  
  # Calculate the distances to centers (matrix C)
  if(is.null(dist_mat) | length(center_ids) == 0){
    
    
    C <- matrix(0, ncol = k, nrow = n)
    if (is.null(params)) {
      for (i in 1:k) {
        C[, i] <- apply(coords,
                        MARGIN = 1,
                        FUN = d,
                        x2 = centers[i, ])
      }
    } else {
      
      for (i in 1:k) {
        centers_i <- centers[i, ]
        params_i <- params[which((centers_i %>% pull(1)) == (coords %>% pull(1)) &
                                   (centers_i %>% pull(2)) == (coords %>% pull(2))),]
        
        C[, i] <- apply(cbind(coords, params),
                        MARGIN = 1,
                        FUN = d,
                        x2 = c(centers_i, params_i))
      }
    }
  } else {
    # Read distances from distance matrix
    C <- dist_mat[,center_ids]
  }
  
  # Use weighted distances
  C <- C * weights
  
  # Gurobi model
  model <- list()
  
  if(g == 1){
    
    # Constraints for the upper and lower limit
    const_LU <- rbind(
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = rep(1:k, each = n),
        j = 1:(n * k),
        x = rep(capacity_weights, k)
      ),
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = rep(1:k, each = n),
        j = 1:(n * k),
        x = rep(capacity_weights, k)
      )
    )
    
    # Add the constraints to the model
    model$A <- rbind(Matrix::spMatrix(
      nrow = n,
      ncol = n_decision,
      i = rep(1:n, times = ifelse(is_outgroup, k + 1, k)),
      j = rep(1:n_decision),
      x = rep(1, n_decision)
    ),
    const_LU)
    
    # Right hand side values
    model$rhs <- c(multip_centers,
                   rep(ranges[1, 2], k),
                   rep(ranges[1, 1], k))
    
    # Model sense
    model$sense      <- c(rep('=', n), 
                          rep('<', k), 
                          rep('>', k))
    
    
  } else {
    # Large number
    M <- 1000
    
    # In constraint matrices: first rangegroups, then clusters
    
    # Constraints for the first lower and upper limit
    const_LU <- rbind(
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = c(rep(1:k, each = n), 1:k),
        j = c(1:(n * k), (n * k) + 1:k),
        x = c(rep(capacity_weights, k), rep(M, k))
      ),
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = c(rep(1:k, each = n), 1:k),
        j = c(1:(n * k), (n * k) + 1:k),
        x = c(rep(capacity_weights, k), rep(-M, k))
      )
    )
    
    # Right hand side values for the first capacity constraints
    rhs_LU <- c(rep(ranges[1, 2] + M, k),
                rep(ranges[1, 1] - M, k))
    
    # Model sense for the first capacity constraints
    sense_LU     <- c(rep('<', k),
                      rep('>', k))
    
    # Constraints, rhs and sense for the rest of the lower and upper limits
    for(i in 2:g){
      const_LU <- rbind(
        const_LU,
        Matrix::spMatrix(
          nrow = k,
          ncol = n_decision,
          i = c(rep(1:k, each = n), 1:k),
          j = c(1:(n * k), (n * k + k * (i - 1)) + 1:k),
          x = c(rep(capacity_weights, k), rep(M, k))
        ),
        Matrix::spMatrix(
          nrow = k,
          ncol = n_decision,
          i = c(rep(1:k, each = n), 1:k),
          j = c(1:(n * k), (n * k + k * (i - 1)) + 1:k),
          x = c(rep(capacity_weights, k), rep(-M, k))
        )
      )
      
      rhs_LU <- c(rhs_LU,
                  rep(ranges[i, 2] + M, k),
                  rep(ranges[i, 1] - M, k))
      
      sense_LU <- c(sense_LU,
                    rep('<', k),
                    rep('>', k))
      
    }
    
    # Constraints for the cluster size group
    const_group <- Matrix::spMatrix(
      nrow = k,
      ncol = n_decision,
      i = rep(1:k, each = g),
      j = c((n * k) + 1:(k * g)),
      x = rep(1, k * g)
    )
    
    # Add all constraints to the model
    model$A <- rbind(
      Matrix::spMatrix(
        nrow = n,
        ncol = n_decision,
        i = rep(1:n, times = ifelse(is_outgroup, k + 1, k)),
        j = rep(1:(n * k)),
        x = rep(1, n * k)
      ),
      const_LU,
      const_group
    )
    
    # Right hand side values (multiple membership, upper and lower limits, cluster groups)
    model$rhs <- c(multip_centers,
                   rhs_LU,
                   rep(1, k))
    
    # Model sense
    model$sense <- c(rep('=', n),
                     sense_LU,
                     rep('=', k))
  }
  
  
  # Objective function
  obj_fn <- c(c(C),
              switch(g > 1, rep(0, k * g), NULL),
              switch(is_outgroup, lambda * weights, NULL))
  
  model$obj <- obj_fn
  
  # Minimization task
  model$modelsense <- 'min'
  
  # B = Binary, C = Continuous
  model$vtype <- ifelse(frac_memb, 'C', 'B')
  
  # Using timelimit-parameter to stop the optimization if time exceeds 10 minutes
  # and disabling the print output from gurobi.
  if(is.null(gurobi_params)){
    gurobi_params <- list()
    gurobi_params$TimeLimit <- 600
    gurobi_params$OutputFlag <- 0  
  }
  
  # Solving the linear program
  result <- gurobi::gurobi(model, params = gurobi_params)
  
  # Send error message if the model was infeasible
  if(result$status == "INFEASIBLE") {stop("Model was infeasible! (rpack)")}
  
  # Returns the assignments
  assign_frac <- Matrix::Matrix(matrix((result$x)[1:(ifelse(is_outgroup, n * k + n, n * k))],
                                       ncol = ifelse(is_outgroup, k + 1, k)), sparse = TRUE)
  
  # Returns the value of the objective function
  obj_min <- round(result$objval, digits = 5)
  
  # Clear space
  rm(model, result)
  
  return(list(assign_frac = assign_frac,
              obj_min = obj_min))
}


#' Updates the parameters (centers) for each cluster.
#'
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights The weights of the data points.
#' @param k The number of clusters.
#' @param assign_frac A vector of cluster assignments for each data point.
#' @param fixed_centers Predetermined center locations.
#' @param d The distance function.
#' @param place_to_point if TRUE, cluster centers will be placed to a point.
#' @param dist_mat Distance matrix for all the points.
#' @return New cluster centers.
#' @keywords internal
location_step <- function(coords,
                          weights,
                          k,
                          assign_frac,
                          params = NULL,
                          fixed_centers = NULL,
                          d = euc_dist2,
                          place_to_point = TRUE,
                          dist_mat = NULL,
                          lambda_fixed = NULL) {
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_centers), 0, nrow(fixed_centers))
  
  # Initialization of cluster center matrix
  centers <- matrix(0, nrow = k, ncol = ncol(coords))
  
  # Add fixed centers first
  if (n_fixed > 0) {
    centers[1:n_fixed, ] <- fixed_centers
  }
  
  if (place_to_point) {
    # Initialization of cluster id vector
    center_ids <- rep(0, k)
    
    # Add fixed centers first
    if (n_fixed > 0) {
      
      center_ids <- prodlim::row.match(fixed_centers %>% as.data.frame(), coords %>% as.data.frame())
      #center_ids <- c(which(((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))
      #) &
      #  ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2))
      #  )), rep(0, k - n_fixed))
    }
  }
  
  # Update center for each cluster
  if (place_to_point) {
    # If fixed servers are allowed to be released
    if (!is.null(lambda_fixed)) {
      # Potential new centers
      potential_center_ids <- rep(0, n_fixed)
      
      # Ids for the fixed servers
      fixed_center_ids <- which(((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))) &
                                  ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2))))
      
      # Most optimal center location for the fixed center clusters
      for (i in 1:n_fixed) {
        # Compute medoids only with points that are relevant in the cluster i
        relevant_cl <- assign_frac[, i] > 0.001
        
        relevant_ids <- which(relevant_cl)
        
        # Computing medoid ids for cluster i
        potential_center_ids[i] <-
          medoid_dist_mat(dist_mat = dist_mat,
                          ids = relevant_ids,
                          w = weights)
        
        # Combined distance to the potential center
        wdist_pot_center <- sum(dist_mat[relevant_ids, potential_center_ids[i]] *
                                  weights[relevant_ids])
        
        # Combined distance to the fixed center
        wdist_fixed_center <- sum(dist_mat[relevant_ids, fixed_center_ids[i]] *
                                    weights[relevant_ids])
        
        # TODO: add lambda_fixed to here
        if (wdist_fixed_center > wdist_pot_center + lambda_fixed) {
          center_ids[i] <- potential_center_ids[i]
        }
        
      }
      
      # Decide the rest of the center locations
      for (i in (n_fixed + 1):k) {
        # Compute medoids only with points that are relevant in the cluster i
        relevant_cl <- assign_frac[, i] > 0.001
        
        # Computing medoid ids for cluster i
        center_ids[i] <- medoid_dist_mat(dist_mat = dist_mat,
                                         ids = which(relevant_cl),
                                         w = weights)
      }
      
      # Decide centers from the ids
      centers <- coords[center_ids, ]
      
    } else {
      
      if(n_fixed < k){
        for (i in (n_fixed + 1):k) {
          # Compute medoids only with points that are relevant in the cluster i
          relevant_cl <- assign_frac[, i] > 0.001
          
          # Computing medoid ids for cluster i
          center_ids[i] <- medoid_dist_mat(dist_mat = dist_mat,
                                           ids = which(relevant_cl),
                                           w = weights)
        }
      }
      
      # Decide centers from the ids
      centers <- coords[center_ids, ]
    }
  } else {
    if(is.null(params)){
      
      for (i in (ifelse(n_fixed > 0, n_fixed + 1, 1)):k) {
        # Check whether euc_dist or euc_dist2 is used
        if (d(0, 2) == 2) {
          w_assign <- weights * assign_frac[, i]
          
          # If only one point belongs to the cluster, the Weiszfeld algorithm won't work
          if (sum(w_assign > 0) == 1) {
            centers[i,] <- coords %>%
              slice(which(w_assign > 0)) %>%
              unlist(., use.names = FALSE)
            
          } else {
            # Weighted median
            weiszfeld <-
              Gmedian::Weiszfeld(coords, weights = w_assign)$median
            
          }
          
        } else if (d(0, 2) == 4) {
          # Weighted mean
          centers[i,] <-
            colSums(coords * weights * assign_frac[, i]) / sum(assign_frac[, i] * weights)
        }
      }
    } else {
      # TODO: Tähän lokaatio parametreilla
      
      
      
    }
    center_ids <- NULL
  }
  
  
  return(list(centers = centers, center_ids = center_ids))
}