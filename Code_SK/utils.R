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


normalize_for_AUC <- function(x, opt = 1){
  # x has table to be normalized including cardinality
  # opt = 1 -> divide by maximum
  # opt = 2 -> divide by the mean
  # opt = 3 or other number, do not divide
  maxyval <- max(x[ ,-1], na.rm = TRUE)
  meanyval <- mean(as.matrix(x[ ,-1]), na.rm = TRUE)
  maxxval <- dim(x)[1]
  if(opt == 1){
    sc <- maxyval
  }else if(opt == 2){
    sc <- meanyval
  }else{
    sc <- 1
  }
  xout <- x
  xout[ ,-1] <- x[ ,-1]/sc
  xout[ ,1] <- x[ ,1]/maxxval
  xout
}
