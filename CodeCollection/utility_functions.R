get_cols <- function(){
  return(c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
           "chocolate","darkred","yellow3","darkgreen","bisque4","magenta",
           "royalblue","tomato4","steelblue1",
           "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
           "coral3","beige","brown4","indianred1","lightgreen","orchid"))
}

plot_cl = function (coords, weights, clusters, centers, title = "", subtitle = NULL, outgroup_label=99, outgroup_legend = "outgroup") 
{
  
  #coords = test_dat[,1:2]
  centers = as.matrix(centers)
  coords = as.matrix(coords)
  x = coords[, 1]
  y = coords[, 2]
  k = nrow(centers)
  w = weights
  centers = as.matrix(centers)
  help_clusters = clusters
  clusters = as.factor(clusters)
  cl_sizes = apply(X = t(1:k), MARGIN = 2, FUN = function(x) {
    sum(w[clusters == x])
  })
  cl_sizes = cl_sizes[help_clusters]
  cl_sizes[is.na(cl_sizes)] = sum(w[help_clusters == outgroup_label])
  cl_sizes[help_clusters == outgroup_label] = 
    paste(outgroup_legend, " (", cl_sizes[help_clusters == outgroup_label][1], ")", sep="")
  cl_sizes = factor(cl_sizes, 
                    levels = unique(cl_sizes)[order(nchar(unique(cl_sizes)), unique(cl_sizes))]) 
  plot = ggplot2::ggplot() + 
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y, size = w, 
                                               colour = cl_sizes)) + 
    ggplot2::scale_size(range = c(2, 7), guide = "none") + 
    ggplot2::scale_color_manual(values = c_col[1:(k+1)], name = "Cluster sizes:") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) + 
    ggplot2::labs(x = "x", y = "y", title = NULL, subtitle = NULL) + 
    ggplot2::theme(legend.position = "right", axis.text.x = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + 
    ggplot2::geom_point(mapping = ggplot2::aes(x = centers[, 1], y = centers[, 2]), size = 3, 
                        col = "black", show.legend = FALSE, shape = 4, stroke = 3) + 
    ggplot2::ggtitle(subtitle) + 
    ggplot2::theme(
      plot.title = element_text(size=30, face="bold", margin = margin(10, 0, 10, 0)),
      axis.text.x = element_text(angle=0, size = 12),
      axis.text.y = element_text(angle = 0, size = 12),
      axis.title = element_text(size=14),
      legend.text = element_text(size=14), 
      legend.title = element_text(size = 16)
    ) +
    coord_fixed()
  
  return(plot)
  
}


#' Determinates the squared distance between two points given their longitudes and latitudes
#'
#' @param x1 Longitude and latitude of point 1.
#' @param x2 Longitude and latitude of point 2.
#'
#' @return Squared euclidean distance between two points.
#' @export
#'
#' @examples
hav.dist2 <- function(x1, x2) {
  long1 <- x1[1]
  lat1 <- x1[2]
  long2 <- x2[1]
  lat2 <- x2[2]
  R <- 6371
  diff.long <- (long2 - long1)
  diff.lat <- (lat2 - lat1)
  a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
  b <- 2 * asin(pmin(1, sqrt(a))) 
  d = R * b
  return(d^2)
}

#' Determinates the squared distance between two points given their longitudes and latitudes and their parameters
#'
#' @param x1 Longitude and latitude of point 1 (coordinates and parameters).
#' @param x2 Longitude and latitude of point 2 (coordinates and parameters).
#' @return Squared euclidean distance between two points.
#' @export
#'
#' @examples
hav.dist2_par <- function(x1, x2, lambda) {

  #lambda <- 0.99

  n_coords <- 2
  x1_coord <- x1[1:n_coords]
  x2_coord <- x2[1:n_coords]
  x1_par <- x1[-(1:n_coords)]
  x2_par <- x2[-(1:n_coords)]


  long1 <- x1_coord[1]
  lat1 <- x1_coord[2]
  long2 <- x2_coord[1]
  lat2 <- x2_coord[2]
  R <- 6371
  diff.long <- (long2 - long1)
  diff.lat <- (lat2 - lat1)
  a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
  b <- 2 * asin(pmin(1, sqrt(a)))

  d_dist <- (R * b)

  d_par <- matrix(c(x1_par, x2_par), nrow = 2, byrow = TRUE) |>
    dist()

  return(lambda * d_dist^2 + (1 - lambda) * d_par^2)
}

plot_point_gradient <- function(dat, par = "par1", high = "#21D6DD", low = "#000000"){
  dat |> 
    ggplot() +
    geom_point(aes_string(x = "x", y = "y", color = par, size = "w")) +
    scale_size(range = c(2, 6)) +
    scale_colour_gradient(
      low = low,
      high = high
    ) +
    guides(size = "none") +
    labs(x = "x", y = "y", color = par) +
    theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )+
    coord_fixed(ratio = 1)
}

#' #' Determinates the squared distance between two points given their longitudes and latitudes and their parameters
#' #'
#' #' @param x1 Longitude and latitude of point 1 (coordinates and parameters).
#' #' @param x2 Longitude and latitude of point 2 (coordinates and parameters).
#' #' @return Squared euclidean distance between two points.
#' #' @export
#' #'
#' #' @examples
# hav.dist2_par <- function(x1, x2, lambda) {
# 
#   #lambda <- 0.99
# 
#   n_coords <- 2
#   x1_coord <- x1[1:n_coords]
#   x2_coord <- x2[1:n_coords]
#   x1_par <- x1[-(1:n_coords)]
#   x2_par <- x2[-(1:n_coords)]
# 
# 
#   long1 <- x1_coord[1]
#   lat1 <- x1_coord[2]
#   long2 <- x2_coord[1]
#   lat2 <- x2_coord[2]
#   R <- 6371
#   diff.long <- (long2 - long1)
#   diff.lat <- (lat2 - lat1)
#   a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
#   b <- 2 * asin(pmin(1, sqrt(a)))
#   d = R * b
#   return(lambda*d^2 + (1-lambda)*(x1_par - x2_par)^2)
# }


plot_nonspat_attributes <- function(dat){
  dat$Y |> 
    ggplot() +
    geom_point(aes(x = x, 
                   y = y, 
                   #size = w, 
                   color = orig_group,
                   fill = orig_group,
                   shape = is_outlier,
                   label = par1), 
               size = 3,
               stroke = 1.5) +
    geom_vline(xintercept = dat$div[1]) + 
    geom_hline(yintercept = dat$div[2]) + 
    guides(size = "none", fill = "none") +
    scale_color_discrete(name = "Original clusters") +
    scale_shape_manual(values = c(21,4)) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           shape = "none")
}

create_data_example <- function(seed = 56117){
  
  set.seed(seed)
  
  data <- simulate_gamma_mixture_n(
    N = 1,
    n_total = 500,
    k = 10,
    n_out = 20,
    out_scale = 5,
    scale_between_range = c(0, 1),
    outgroup_alpha = 0.4,
    place_on_grid = TRUE,
    overlap_scale = 0.5
  )
  
  return(data[[1]])
}

plot_clust_list <- function(list_dat, list_clust, data_id = 1){
  # plot_clusters(
  #   coords = list_dat[[1]]$Y |> dplyr::select(x,y),
  #   weights = list_dat[[1]]$Y |> pull(w),
  #   clusters = list_clust[[1]][[1]]$clusters,
  #   centers = list_clust[[1]][[1]]$centers
  # ) + 
  #   geom_vline(xintercept = list_dat[[1]]$div[1]) + 
  #   geom_hline(yintercept = list_dat[[1]]$div[2]) +
  #   
  
  
    clust_id <- c(3,5,7,9)
    
    plot_clusters(
      coords = list_dat[[1]]$Y |> dplyr::select(x,y),
      weights = list_dat[[1]]$Y |> pull(w),
      clusters = list_clust[[1]][[clust_id[1]]]$clusters,
      centers = list_clust[[1]][[clust_id[1]]]$centers
    ) + 
    geom_vline(xintercept = list_dat[[1]]$div[1]) + 
    geom_hline(yintercept = list_dat[[1]]$div[2]) +
    ggtitle(paste("lambda_par =", list_clust[[1]][[clust_id[1]]]$lambda_par) ) +


    plot_clusters(
      coords = list_dat[[1]]$Y |> dplyr::select(x,y),
      weights = list_dat[[1]]$Y |> pull(w),
      clusters = list_clust[[1]][[clust_id[2]]]$clusters,
      centers = list_clust[[1]][[clust_id[2]]]$centers
    ) +
    geom_vline(xintercept = list_dat[[1]]$div[1]) +
    geom_hline(yintercept = list_dat[[1]]$div[2]) +
    ggtitle(paste("lambda_par =", list_clust[[1]][[clust_id[2]]]$lambda_par) ) +

    plot_clusters(
      coords = list_dat[[1]]$Y |> dplyr::select(x,y),
      weights = list_dat[[1]]$Y |> pull(w),
      clusters = list_clust[[1]][[clust_id[3]]]$clusters,
      centers = list_clust[[1]][[clust_id[3]]]$centers
    ) +
    geom_vline(xintercept = list_dat[[1]]$div[1]) +
    geom_hline(yintercept = list_dat[[1]]$div[2]) +
    ggtitle(paste("lambda_par =", list_clust[[1]][[clust_id[3]]]$lambda_par)) +

      plot_clusters(
        coords = list_dat[[1]]$Y |> dplyr::select(x,y),
        weights = list_dat[[1]]$Y |> pull(w),
        clusters = list_clust[[1]][[clust_id[4]]]$clusters,
        centers = list_clust[[1]][[clust_id[4]]]$centers
      ) +
      geom_vline(xintercept = list_dat[[1]]$div[1]) +
      geom_hline(yintercept = list_dat[[1]]$div[2]) +
    ggtitle(paste("lambda_par =", list_clust[[1]][[clust_id[4]]]$lambda_par) )
}

plot_division_example <- function(dat){
  
  true_mu <- dat$mu_true
  # Ids in interactive plot
  id <-  1:nrow(dat$Y)
  
  div_x <- dat$div[1]
  div_y <- dat$div[2]
  
  dat_division <- tibble(x = c(div_x/2, (1+div_x)/2, div_x/2, (1+div_x)/2),
                         y = c((1+div_y)/2, (1+div_y)/2, div_y/2, div_y/2), 
                         label = c("1", "2", "3", "4"))
  
  plot_sim <- ggplot() +
    geom_point(data = dat$Y, aes(x = x, y = y, size = w, label = id)) +
    scale_size(range = c(2, 6)) +  # Scale objects sizes
    guides(
      color = guide_legend(        # Point size in legend
        override.aes = list(size=5)
      )
    ) +
    labs(x = "x", y = "y", title = "Unclustered data") +
    theme(
      legend.position = "right",            # Legend position and removing ticks from axis
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  plot_sim + 
    geom_hline(data = NULL, yintercept = div_y, color = "red") +
    geom_vline(data = NULL, xintercept = div_x, color = "red") +
    geom_label(data = dat_division, 
               aes(x = x, y = y, label = label),
               size = 10) +
    geom_label(data = NULL, aes(x = c(div_x, 0.025 ),
                                y = c(0, div_y),
                                label = round(c(div_x, div_y), 2)),
               size = 4) +
    guides(color = "none",
           size = "none")+
    ggtitle("") +
    coord_fixed(ratio = 1)
}
