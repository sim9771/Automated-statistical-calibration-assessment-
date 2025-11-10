bi_k_means <- function(data, target_power, eps, alpha) {
  
  # Run k-means with k=2
  K <- 2
  km <- kmeans(data, centers = K, nstart = 25, iter.max = 1000)
  clusters <- km$cluster
  centers <- km$centers
  
  # Compute expected powers for the first two clusters
  cluster_powers <- sapply(1:K, function(g) {
    expected_power(p = centers[g, ], k = sum(clusters == g),
                   m = 3000, eps = eps, alp = alpha)[[1]]
  })
  names(cluster_powers) <- 1:K
  
  results <- list()
  
  # Recursive bi-k-means
  repeat {
    
    avg_power <- mean(cluster_powers)
    
    # Save current iteration
    results[[K-1]] <- list(
      clusters = clusters,
      centers = centers,
      cluster_powers = cluster_powers,
      avg_power = avg_power,
      num_clusters = K
    )
    
    if (max(cluster_powers) < 0.8) break
    
    # Identify the cluster with the highest expected power
    to_split_id <- as.integer(names(which.max(cluster_powers)))
    idx_to_split <- which(clusters == to_split_id)
    subset_data <- data[idx_to_split,]
    if (nrow(subset_data) < 2) break
    
    # run k-means to split the selected cluster
    bi_split <- kmeans(subset_data, centers = 2, nstart = 25, iter.max = 1000)
    sub_clusters <- bi_split$cluster
    sub_centers <- bi_split$centers
    
    # Re-index
    clusters[clusters > to_split_id] <- clusters[clusters > to_split_id] + 1
    clusters[idx_to_split] <- ifelse(sub_clusters == 1, to_split_id, to_split_id + 1)
    
    # Update centers
    centers <- rbind(
      if (to_split_id > 1) centers[1:(to_split_id - 1),],
      sub_centers,
      if (to_split_id < K) centers[(to_split_id + 1):K,]
    )
    
    # Compute expected power for the two new clusters
    new_EP <- sapply(0:1, function(i) {
      id <- to_split_id + i
      expected_power(p = centers[id, ], k = sum(clusters == id),
                     m = 3000, eps = eps, alp = alpha)[[1]]
    })
    
    # Update new expected powers
    cluster_powers[to_split_id] <- new_EP[1]
    cluster_powers <- append(cluster_powers, new_EP[2], after = to_split_id)
    names(cluster_powers) <- seq_along(cluster_powers)
    
    # Increment cluster count and iteration
    K <- K + 1
  }
  
  avg_powers <- sapply(results, function(x) x$avg_power)
  best_idx <- which.min(abs(avg_powers - target_power))
  best_result <- results[[best_idx]]
}
