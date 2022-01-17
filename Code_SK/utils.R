average_over_epsilon <- function(x){
  # x is a list
  # Each element in the list has a scenario
  require(dplyr)
  y <- x
  len <- length(x)
  for(ii in 1:len){
    xx <- x[[ii]]
    y[[ii]] <- xx %>% group_by(Cardinality) %>% summarize(CV_mean_regret = mean(CV_mean_regret, na.rm = TRUE), CV_mean_entropy = mean(CV_mean_entropy, na.rm = TRUE), CV_std_regret = mean(CV_std_regret, na.rm = TRUE), CV_std_entropy = mean(CV_std_entropy, na.rm = TRUE) )
  }
  return(y)
}
