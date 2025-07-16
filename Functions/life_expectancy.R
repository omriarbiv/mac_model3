life_expectancy <- function(state_matrix) {

  # Remove baseline state
  state_matrix <- state_matrix[,-1]
  
  # Number of individuals
  n <- dim(state_matrix)[1]
  max_life <- dim(state_matrix)[2]
  
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
  
  return(life_expect)
}