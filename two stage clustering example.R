# Performs automated calibration assessment for multi-class probabilistic classifiers.
# 
# The function evaluates model calibration by computing statistical test results, 
# expected powers, and clustering outcomes based on class probability estimates.
#
# Parameters:
#   data   : A matrix or data frame of predicted class probabilities 
#   true   : A matrix of true class labels corresponding to 'data'.
#   target_power : Desired statistical power (1 - beta)
#   eps    : L1 distance between p and a sampled probability vector q.
#   alpha  : Significance level in hypothesis testing.
#   tol    : Stopping criterion for the stage-2 clustering procedure.
#
# Returns:
#   A list containing:
#     - test_results : Statistical test outcomes for assessing calibration.
#     - expected_powers : Estimated expected powers for each cluster.
#     - clusters : Final clustering assignments and centroids in the probability simplex
#

auto_assessment <- function(data, true, m_list, target_power = target_power, tol = tol, eps = eps, alpha = alpha) {
  
  # Stage 1
  stage1 <- bi_k_means(data, target_power = target_power, eps = eps, alpha = alpha)
  
  # Stage 2
  
  stage2 <- SC_clust(data, centers = stage1$centers, target_power = target_power, tol = tol, eps = eps, alpha = alpha)
  
  # Assessment
  clustering <- data.frame(data, cluster = stage2$assignments)
  
  c <- length(stage1$num_clusters)
  pval <- numeric(c)
  for (g in 1:c) {
    idx <- which(clustering$cluster == g)
    pval[g] <- round(pearson_test(pred_all = clustering[idx, ][, -4], subs = true[idx, ]), 4)
  }
  
  
  return(list(
        power = stage2$power,
        avg_clust_size = mean(table(clustering$cluster)),
        clust_num = c,
        p_value = pval,
        clutering = clustering))
}


Result <- auto_assessment(pred5, true5, m_list, target_power = 0.8,tol = 0.01, eps = 0.1, alpha = 0.05)

