yearly_qaly <- function(util_matrix, yearly_discount = 0) {
  # Assume monthly cycles to calculate number of years
  years <- (ncol(util_matrix) - 1) / 12
  
  # Initialize empty matrix
  util_year <- matrix(NA, nrow = nrow(util_matrix), ncol = years)
  colnames(util_year) <- paste0("util_y", 1:years)
  
  for(i in 1:years) {
    discount_vec <- (1 / (1 + (yearly_discount / 12)) ^ (0:(i * 12)))
    util_year[, i] <- util_matrix[, 0:((i * 12) + 1)] %*% discount_vec
  }
  
  return(util_year)
}

life_expectancy <- function(state_matrix, start_age = 70) {
  
  # Remove baseline state
  state_matrix <- state_matrix[,-1]
  
  # Number of individuals
  n <- nrow(state_matrix)
  max_life <- ncol(state_matrix)
  
  # Index all death
  ind <- which(state_matrix == "Death", arr.ind = TRUE)
  
  # Make sure cases are in order by individual
  ord <- order(ind[, 1])
  sm_ord <- ind[ord, 2] - 1
  pt_ord <- as.double(names(sm_ord))
  
  occur1 <- match(1:n, pt_ord)
  life_expect <- unname(sm_ord[occur1])
  
  # Correct for those that did not die
  life_expect[is.na(life_expect)] <- max_life
  
  return((life_expect / 12) + start_age)
}
