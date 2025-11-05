library(ompr)
library(ompr.roi)
library(ROI.plugin.symphony)
library(Matrix)
library(dplyr)

SC_clust <- function(data, centers, target_power = target_power, tol = tol, eps = eps, alpha = alpha) {
  
  N <- nrow(data)
  current_means <- centers
  K <- nrow(current_means)
  
  repeat {
    
    # Calculate desired size for each cluster
    new_optimal_sizes <- new_power <- numeric(K)
    for (k in 1:K) {
      expected_power_values <- sapply(m_list, function(m_val) {
        out <- expected_power(p = current_means[k, ], k = m_val, m = 3000, eps = eps, alp = alpha)
        return(out[[1]])
      })
      opt_idx <- which.min(abs(expected_power_values - target_power))
      new_optimal_sizes[k] <- m_list[opt_idx]
      new_power[k] <- expected_power_values[opt_idx]
    }
    adj_cur <- adjust_optimal_size(new_optimal_sizes, N)
    
    # Calculate EP
    
    finalpower <- numeric(K)
    
    for (k in 1:K) {
      out1 <- expected_power(p = current_means[k, ], k = adj_cur[k], m = 3000, eps = eps , alp = alpha)
      finalpower[k] <- out1[[1]]
    }
    
    # Compute Euclidean distance matrix
    eudi <- vector("list", K)
    for (k in 1:K) {
      dist <- data - matrix(current_means[k, ], nrow = N, ncol = ncol(data), byrow = TRUE)
      eudi[[k]] <- sqrt(rowSums(dist^2))
    }
    D_sparse <- Matrix(do.call(cbind, eudi), sparse = TRUE)
    
    # Formulate MILP model
    model <- MIPModel() %>%
      add_variable(x[i, k], i = 1:N, k = 1:K, type = "binary") %>%
      # Each data point must be assigned to exactly one cluster
      add_constraint(sum_expr(x[i, k], k = 1:K) == 1, i = 1:N) %>%
      # Each cluster must have exactly desired number of cluster
      add_constraint(sum_expr(x[i, k], i = 1:N) == adj_cur[k], k = 1:K) %>%
      # Objective is to minimize the total distance
      set_objective(sum_expr(D_sparse[i, k] * x[i, k], i = 1:N, k = 1:K), sense = "min")
    
    # Solve MILP
    solution <- solve_model(model, with_ROI(solver = "symphony"))
    assignments <- get_solution(solution, x[i, k])
    
    # Get assignment results
    assignment_matrix <- matrix(assignments$value, nrow = N, ncol = K)
    assigned_points <- apply(assignment_matrix, 1, which.max)
    
    # Update centers
    new_means <- matrix(0, nrow = K, ncol = ncol(data))
    for (k in 1:K) {
      idx <- which(assigned_points == k)
      new_means[k, ] <- colMeans(data[idx,])
    }
    
    # Check stopping criterion
    Delta <- max(apply(abs(new_means - current_means), 1, sum))
    
    if (Delta < tol) {
      break
    }
    
    # Update centers
    current_means <- new_means
  }
  
  list(assignments = assigned_points, centroids = current_means, power = finalpower)
}