adjust_optimal_size <- function(desired_size, sum) {
  d <- sum - sum(desired_size)
  adjustment <- sign(d) * (abs(d) %/% length(desired_size))
  remainder <- abs(d) %% length(optimal_size)
  new_size <- desired_size + adjustment + sign(d) * (seq_along(desired_size) <= remainder)
  new_size
}