
# Tests for exponentiality

exp_test <- function(x, method = "transf", N = 10^3){  
  DNAME <- deparse(substitute(x))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  if (!is.numeric(x) || length(x) <= 1) warning(paste(DNAME, "must be a numeric vector containing more than 1 observation"))
  x <- as.vector(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  n <- length(x)   # adjusted sample size without NA values
  if (min(x) < 0) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  if(any(c("ratio", "transf") == method) == FALSE) stop("Invalid method. \nValid methods are 'transf'  and 'ratio'. ")  
  alternative = paste(DNAME," does not follow an Exponential distribution.")
  # Test based on a transformation to approximately uniform r.v.
  if(method == "transf"){
    stat <- function(x){
      n <- length(x)
      b.check  <-   cov(x, log(x))
      u  <- exp(-x / b.check)
      variance <- (1 + trigamma(1)) / 16 + 1 / 12 + (log(2) - 1) / 8
      t  <- sqrt( n / variance ) * (mean(u) - 1 / 2)
      return(t)
    }  
    stat_c <- stat(x)
    if (n >= 200){
      p1 <- pnorm(stat_c)
      p2 <- 1 - p1
      p.value <- 2*min(p1, p2)
    }
    if( n < 200){
      null.distr  <- replicate(N, stat(rexp(n, rate = 1 ))) 
      p1          <- sum(null.distr < stat_c) / N
      p2          <- 1 - p1
      p.value     <- 2 * min(p1, p2)
    }
    results <- list(statistic = c("T" = stat_c), p.value = p.value, data.name = DNAME,
                    method = "Test for exponentiality based on a transformation to uniformity ")
  }

  # Cox-Oakes test (based on the ratio of two scale estimators).
  if (method == "ratio"){
    co_stat <- function(x){
      m <- mean(x)
      lo <- log(x / m)
      l <- log(x)
      v <- (n + (1/m) * sum(x * (lo)^2) - (sum(x * lo))^2 / (n * m^2))^(-1)
      u <- n + sum(l) - sum(x * l) / m
      obs_stat <- sqrt(v) * u
      return(obs_stat)
    }
    stat_c <- co_stat(x)
    p1 <- pnorm(stat_c)
    p2 <- 1 - p1
    p.value <- 2 * min(p1, p2)
    results <- list(statistic = c("CO" = stat_c), p.value = p.value, 
                    method = "Cox-Oakes test for exponentiality", data.name = DNAME)
  }
  class(results) <- "htest"
  return(results)
}

  
  
