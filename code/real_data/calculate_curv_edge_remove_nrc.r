### Code for calculating NRC and removing edges on real data ### 

# Load required libraries
library(ggplot2)
library(igraph)
library(R.matlab)
source("/code/linear_curv.R")


datasets <- c("caltech", "simmons", "polblog")

n <- 0
for (dataname in datasets){
  graph_path <- paste0("/data/", dataname, "_curvs.gml")
  output_path <- paste0("/result/", dataname, "/nrc")

  # Create subfolders for each metric
  adj_dir <- file.path(output_path, metric, "adj") # Folder to save adjacency matrix
  fig_dir <- file.path(output_path, metric, "fig") # Folder to save edge removal plot
  dir.create(adj_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    
  # Load graph data (one with BFC calculated)
  g <- read_graph(graph_path, format = "gml")
  V(g)$membership <- V(g)$community # Put groundtruth membership data in V(g)$membership
  g <- cal_linear_curv(g)
  
  # Construct data frame on edge attributes
  edge_list <- as.data.frame(get.edgelist(g))
  data <- cbind(
    edge_list,
    nrc = E(g)$nrc,
    same_group = E(g)$membership,
    minni = E(g)$minni
  )
  colnames(data) <- c("node_i", "node_j", "nrc", "same_group", "minni")
  data$sqrtmin <- 1 / sqrt(E(g)$minni)
  data$same_group <- as.factor(data$same_group)
  
  # Add sqrtmin as edge attribute
  E(g)$sqrtmin <- 1 / sqrt(E(g)$minni)
  
  # Keep track of edge number
  total_red <- sum(data$same_group == "Different")
  total_blue <- sum(data$same_group == "Same")
  total_edges <- nrow(edge_list)
  
  # Generate weighted matrix for nrc
  nrc_matrix <- as_adjacency_matrix(g, attr="nrc", sparse=FALSE)
  nrc_matrix[nrc_matrix == 0] <- NA  
  zero_weight_edges <- which(E(g)$nrc == 0)  
  edge_list_zero <- ends(g, E(g)[zero_weight_edges])  
  nrc_matrix[cbind(edge_list_zero[,1], edge_list_zero[,2])] <- 0
  nrc_matrix[cbind(edge_list_zero[,2], edge_list_zero[,1])] <- 0 
  
  # Generate weighted matrix for sqrtmin
  sqrtmin_matrix <- as_adjacency_matrix(g, attr="sqrtmin", sparse=FALSE)
  sqrtmin_matrix[sqrtmin_matrix == 0] <- NA  
  sqrtmin_matrix[cbind(edge_list_zero[,1], edge_list_zero[,2])] <- 0
  sqrtmin_matrix[cbind(edge_list_zero[,2], edge_list_zero[,1])] <- 0  

  # Membership info in matrix
  node_membership <- V(g)$membership 
  same_group_matrix <- outer(node_membership, node_membership, FUN = "==") 

  # Get original adj matrix
  original_adj_matrix <- as_adjacency_matrix(g, sparse = TRUE)
    
  # Define grid search parameters using theta and x-intercept
  thetas <- seq(88, 96, length.out = 50)
  avg_degree <- mean(degree(g))
  xints <- seq(-1/mean(degree(g))*0.25, 1/mean(degree(g)), length.out = 50)      # x-intercept values
  
  
  results <- list()
  original_adj_matrix <- as_adjacency_matrix(g, sparse = TRUE)
  
  for (theta in thetas) {
    slope <- tan(theta * pi / 180)
    
    for (x_int in xints) {
      print(paste0("Theta is ", theta, " degrees and x-intercept is ", x_int))
      n <- n + 1
      x_boundary <- function(y) y / slope + x_int
      
      
      mask <- nrc_matrix > x_boundary(sqrtmin_matrix)
      filtered_nrc <- as.matrix(nrc_matrix)
      filtered_nrc[!mask] <- NA 
      filtered_adj_matrix <- as.matrix(original_adj_matrix)  
      filtered_adj_matrix[!mask] <- 0  
      
      # Sanity check: if locations of NA in filtered_nrc match 0 in filtered_adj_matrix
      sanity_check <- which(is.na(filtered_nrc)) == which(filtered_adj_matrix == 0)
      if (all(sanity_check, na.rm = TRUE)) {
        print("Sanity check passed: NA locations in filtered_nrc match 0 locations in filtered_adj_matrix!")
      } else {
        print("Sanity check failed: Some locations do not match!")
      }
      
      # Membership info in matrix
      node_membership <- V(g)$membership 
      same_group_matrix <- outer(node_membership, node_membership, FUN = "==") 
      
      # Count total red (Different) and blue (Same) edges before filtering
      total_red <- sum(same_group_matrix == 0 & original_adj_matrix > 0, na.rm = TRUE)
      total_blue <- sum(same_group_matrix == 1 & original_adj_matrix > 0, na.rm = TRUE)
      total_edges <- sum(original_adj_matrix > 0, na.rm = TRUE)  # Total edges before filtering
      
      # Count remaining red and blue edges after filtering
      remaining_red <- sum(same_group_matrix == 0 & filtered_adj_matrix > 0, na.rm = TRUE)
      remaining_blue <- sum(same_group_matrix == 1 & filtered_adj_matrix > 0, na.rm = TRUE)
      
      # Compute removal ratios
      removed_red_ratio <- 1 - (remaining_red / total_red)
      removed_blue_ratio <- 1 - (remaining_blue / total_blue)
      edges_removed <- total_edges - sum(filtered_adj_matrix > 0, na.rm = TRUE)
      
      # Print results
      print(paste0("Removed Red Ratio: ", removed_red_ratio))
      print(paste0("Removed Blue Ratio: ", removed_blue_ratio))
      print(paste0("Total Edges Removed: ", edges_removed))
      
      results <- append(results, list(
        list(theta = theta, x_intercept = x_int, removed_diff_ratio = removed_red_ratio, removed_same_ratio = removed_blue_ratio, edges_removed = edges_removed, total_diff = total_red, total_same = total_blue,
             total_edges = total_edges)
      ))
      
      # Save adj matrix
      mat_file <- file.path(adj_dir, paste0("adj__matrix_theta_", round(theta, 2), "_xint_", round(x_int, 4), ".mat"))
      writeMat(mat_file, adj_matrix = filtered_adj_matrix)
      

      # Prepare edge removal data
      edge_indices <- which(!is.na(nrc_matrix), arr.ind = TRUE)  
      plot_data <- data.frame(
        node_i = edge_indices[, 1],  
        node_j = edge_indices[, 2],  
        nrc = nrc_matrix[edge_indices], 
        sqrtmin = sqrtmin_matrix[edge_indices],  
        same_group = factor(ifelse(same_group_matrix[edge_indices] == 1, "Same", "Different"))  
      )
      
      removed_mask <- filtered_adj_matrix[edge_indices] == 0  
      all_removed <- all(removed_mask)  
      plot_data$shape_type <- ifelse(removed_mask, "Removed", "Kept")  
      if (all_removed) {
        plot_data$shape_type <- "Removed"
      }
      
      # Plot edge removal data & save the plot
      plot <- ggplot(plot_data, aes(x = nrc, y = sqrtmin, color = same_group)) +
        geom_point(aes(shape = factor(shape_type)), alpha = 0.7, stroke = 1.1, size = 3) +  
        geom_abline(slope = slope, intercept = -slope * x_int, color = "blue", linewidth = 1.2) +
        labs(
          title = paste("Decision Boundary - Theta:", round(theta, 2), "°", "x-intercept:", round(x_int, 4)),
          subtitle = paste("Removed Diff Ratio:", round(removed_red_ratio, 3),
                           "| Removed Same Ratio:", round(removed_blue_ratio, 3)),
          x = expression(NRC),
          y = expression(SD),
          shape = "Edge Status", color = "Community Membership"
        ) +
        scale_shape_manual(values = c("Kept" = 16, "Removed" = 21)) +
        theme_minimal()

      fig_file <- file.path(fig_dir, paste0("plot_theta_", round(theta, 2), "_xint_", round(x_int, 4), ".pdf"))
      ggsave(fig_file, plot, width = 8, height = 6)
    }
  }
  
  # Save edge removal result in csv
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  write.csv(results_df, file.path(output_path, metric, paste0("grid_results_", metric, ".csv")), row.names = FALSE)
}

