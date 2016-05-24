# ----------------------------------------------------------------------------------------
# Test for Cauchy distributions

cauchy.test <- function(x, N = 10^3){
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x)) stop(paste(DNAME, "must be a numeric vector"))
  
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   # adjusted sample size without NA values
  if ( n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  
  #test statistic
  stat <- function(x){
    estim <- mledist(x, distr = "cauchy")$estimate
    theta.hat <- as.numeric(estim[1])
    l.hat    <- as.numeric(estim[2])
    t <- l.hat / mean(abs(x - theta.hat))  # ratio of two scale estimators
    return(t)  
  }
  stat_c <- stat(x)  # observed value of the test statistic
  p.value <- sum(replicate(N, stat(rcauchy(n))) > stat_c) / N  # p-value is approximated by Monte Carlo simulation
  method  <- "Test for the Cauchy distribution based on the ratio of two scale estimators "
  results <- list(statistic = c("T" = stat_c), p.value = p.value, method = method, data.name = DNAME)
  class(results) = "htest"
  return(results)
}