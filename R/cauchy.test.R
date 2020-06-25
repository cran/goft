# ----------------------------------------------------------------------------------------
# Test for Cauchy distributions

cauchy_test <- function(x, N = 10^3,  method = "transf"){
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x)) stop(paste(DNAME, "must be a numeric vector"))
  
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   # adjusted sample size without NA values
  if ( n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if(any(c("ratio","transf") == method)== FALSE) stop("Invalid method. \nValid methods are 'transf'  and 'ratio'. ")  
  
  if(method == "transf"){
    t_c     <- .ad.stat.exp(x) 
    p.value <- sum(replicate(N, .ad.stat.exp(rcauchy(n))) > t_c) / N
    method  <- "Test for the Cauchy distribution based on a  transformation to exponential data"
    results <- list(statistic = c("T" = t_c), p.value = p.value, method = method, data.name = DNAME)
    class(results) = "htest"
    return(results)
    }
  
  # Test based on the ratio of two variance estimators.
  if(method == "ratio"){
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
  
  
}


.ad.stat.exp <- function(x){
  estim    <- mledist(x,distr="cauchy")$estimate
  teta.hat <- estim[1]
  l.hat    <- estim[2]
  Fx       <- pcauchy(x, teta.hat, l.hat)
  ex       <- -log(Fx)   # transformation to approx exponential data
  n        <- length(ex)
  b        <- mean(ex)
  s        <- sort(ex)
  theop    <-  pexp(s,1/b)
  # Anderson-Darling statistic
  ad       <- -n-sum((2*(1:n)-1)*log(theop) + (2*n+1-2*(1:n))*log(1-theop))/n
  return(ad)
}
