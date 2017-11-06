# -------------------------------------------------------------------------------------------------------
# Correlation and variance ratio tests for the extreme value distributions (Gumbel, Frechet and Weibull)
# N is the number of Monte Carlo samples used to approximate critical values of the variance ratio test
# -------------------------------------------------------------------------------------------------------
ev_test <- function(x, dist = "gumbel", method = "cor", N = 1000){  
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x)) stop(paste(DNAME, "must be a numeric vector"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   # adjusted sample size without NA values
  if (n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  
  if (any(c("gumbel", "frechet", "weibull") == dist) == FALSE) stop("Invalid distribution. \nValid distributions are 'gumbel', 'frechet' and 'weibull'. ")  
  
  if (dist == "frechet" && min(x) < 0) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers when testing the Frechet distribution hypothesis.")
  if (dist == "weibull" && max(x) > 0) stop("The dataset contains positive observations. \nAll data must be negative real numbers when testing the Weibull extreme value distribution hypothesis.")
  
  if (any(c("cor", "ratio") == method) == FALSE) stop("Invalid method. \nValid methods are 'cor'  and 'ratio'. ")  
  
  
  if( dist == "frechet")  x <- log(x)
  if( dist == "weibull")  x <- -log(-x)
  
  if(method == "cor"){
      x <- -x
      if (n < 20 || n > 250) stop("The correlation test requires a sample size between 20 and 250.")
    
    cor.stat <- function(x){
      z   <- combn(x,2,max)  
      nz  <- length(z)
      Fn  <- rep(NA,nz)
      Fn  <- ecdf(z)
      FnZ <- Fn(z)
      y   <- rep(NA,nz)
      y   <- log(qweibull(FnZ, 1, lower.tail = TRUE))  
      r   <- NA
      logic  <- (z!= max(z))
      r      <- cor(z[logic], y[logic], method = "pearson")
      rl     <- log(1-r)
      return(rl)
    }
      rl <- cor.stat(x)
      median_n <- -3.02921 - .03846 * n + .00023427 * n**2 - .000000471091 * n**3
      if (n <= 60 ){
        s_n <- 0.7588 - 0.01697 * n + 0.000399 * n**2 - 0.000003 * n**3
      }
      if (n > 60) s_n <- .53
      
      p.value <- pnorm(rl, mean = median_n, sd = s_n, lower.tail = FALSE)
      r <- 1 - exp(rl)
      results <- list(statistic = c(R = r), p.value = p.value, 
                      method = paste("Correlation test of fit for the ", dist, " distribution"), 
                      data.name = DNAME)
  }
  
  if(method == "ratio"){
    ratio.stat <- function(x){
      m    <- mean(x)
      frac <- 1 / (1:n)
      summ <- rep(NA, n)
      for(i in 1:n){
        summ[i] <- sum(frac[i:n])
      }
      y <- sort(x)
      s.kim <- m - mean(y * summ)             # Kimball's estimator for the scale parameter
      t <- (pi * s.kim)^2 / (6 * var(x))      # variance ratio test statistic  
      return(t)
    } 
    ratio.stat_c <- ratio.stat(x)    # observed value of the test statistic
    null.dist <- replicate(N, ratio.stat(-log(rexp(n))))  
    p1 <- sum(ratio.stat_c < null.dist) / N
    p2 <- sum(ratio.stat_c > null.dist) / N
    pvalue <- 2 * min(p1, p2)    # approximated p-value by Monte Carlo simulation
    results <- list(statistic = c(T = ratio.stat_c), p.value = pvalue, 
                    method = paste("Variance ratio test of fit for the ", dist, " distribution"), 
                    data.name = DNAME)
  }
  class(results) = "htest"
  return(results)
}

