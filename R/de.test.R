# ------------------------------------------------------------------------------------
# Tests for the Laplace or double exponential distribution

laplace_test <- function(x, method = "transf", N = 10^5){
  dname <- deparse(substitute(x))  
  if (!is.numeric(x)) stop(paste(dname, "must be a numeric vector"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   # adjusted sample size without NA values
  if (n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  
  if(any(c("ratio", "transf") == method) == FALSE) stop("Invalid method. \nValid methods are 'transf'  and 'ratio'. ")  
  
  # Test based on a data transformation to approximately iid exponential random variables
  # Anderson-Darling test for exponentiality is used
  if (method == "transf"){
    z <- abs(x - mean(x))   # transformed data
    n <- length(z)
    theop <-  pexp(sort(z), 1 / mean(z)) # theoretical CDF
    # Anderson-Darling statistic for testing exponentiality
    ad <-  -n - sum((2 * (1:n) - 1) * log(theop) + (2 * n + 1 - 2 * (1:n)) * log(1 - theop)) / n
    m <- .6
    l <- 1.6
    # Approximated p-value (the upper tail of AD's null distribution is approximated by an IG(m,l) distribution)
    p.value <- 1 - (pnorm(sqrt( l / ad) * (ad / m - 1)) + exp(2 * l / m) * pnorm( -sqrt( l / ad) * (ad / m + 1)))
    results <- list(statistic = c("AD" = ad), p.value = p.value, method = "Test for the Laplace distribution based on a transformation to exponentiality", data.name = dname)                  
  }
  
  # Test based on the ratio of two estimators for the scale parameter of the Laplace distribution
  if (method == "ratio"){
    stat  <- function(x){
      theta.tilde <- mean(x)
      beta.check  <- sum(abs(x - theta.tilde)) / n                # modified MLE of the scale parameter b (absolute mean deviation)
      beta.tilde  <- sqrt(sum((x - theta.tilde)^2) / (2 * n))     # MME of the scale parameter b
      r           <- beta.tilde / beta.check              # ratio of scale estimators
      t           <- sqrt(4 * n) * (r-1)                  # test statistic
      return(t)
    }
    t <- stat(x)
    if ( n < 500){
      null.distr <- replicate(N, stat(rexp(n) - rexp(n)))  
      p1 <-  sum(null.distr < t) / N
      p.value <- 2 * min(p1, 1 - p1)  
    } 
    if (n >= 500) p.value <- 2 * pnorm(abs(t), lower.tail = FALSE)
    results <- list(statistic = c(R = t), p.value = p.value, method = "Ratio test of fit for the Laplace distribution", data.name = dname)
  }
  
  class(results) <- "htest"
  return(results)
}



