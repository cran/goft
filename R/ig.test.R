# Tests for inverse Gaussian distributions

ig.test <- function(x, method = "transf"){  
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x) & length(x) <= 1) stop(paste(DNAME, "must be a numeric vector containing more than 1 observation"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  n <- length(x)   # adjusted sample size without NA values
  if (min(x) < 0) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  if(any(c("ratio","transf") == method)== FALSE) stop("Invalid method. \nValid methods are 'transf'  and 'ratio'. ")  
  alternative = paste(DNAME," does not follow an Inverse Gaussian distribution.")

    
  if(method == "transf"){
    # Test based on a transformation to normality
    u <- rbinom(n, 1, 0.5)
    y <- abs((x - mean(x)) / sqrt(x))
    result <- shapiro.test(y * (1 - u) + y * (-u))
    results1 <- list(statistic = result$statistic, p.value = result$p.value, data.name = DNAME,
                     method = "Test for Inverse Gaussian distributions using a transformation to normality") 
    class(results1) = "htest"
    
    # Test based on a transformation to gamma variables
    p.value   <- NA
    statistic <- NA
    z <- ((x - mean(x))**2)/x  # transformation given in equation (6) of Villasenor and Gonzalez-Estrada (2015)
      res <- .gammadist.test2(z)
      ad <- res$statistic
      p.value <- res$p.value
      results2 <- list(statistic =  ad, p.value = p.value, data.name = DNAME,
                       method = "Test for Inverse Gaussian distributions using a transformation to gamma variables")
      class(results2) = "htest"
      return(list(results1, results2))  
    }      
    
  # Test based on the ratio of two variance estimators.
  if(method == "ratio"){
      m  <- mean(x)
      v  <- mean(1 / x - 1 / m)
      s2 <- var(x)
      r  <- s2 / (m^3 * v)
      t  <- sqrt(length(x) / (6 * m * v)) * (r - 1)       # test statistic (see Corollary 1 of Villasenor and Gonzalez-Estrada (2015))
      p.value <- 2 * pnorm(abs(t), lower.tail = FALSE)
      results <- list(statistic = c(T1 = t), p.value = p.value, 
                      method = " Variance ratio test for Inverse Gaussian distributions", data.name = DNAME)
      class(results) = "htest"
      return(results)
    }
}

# Anderson-Darling test for the gamma distribution with shape parameter equal to 0.5
# and unknown scale parameter
.gammadist.test2 <- function(z){      
  n  <- length(z)
  b. <- 2 * mean(z)
  s  <- sort(z) / b.  
  theop <-  pgamma(s, .5)
  # Anderson-Darling statistic
  ad <- - n - sum(( 2 * (1 : n) - 1) * log( theop ) + (2 * n + 1 - 2 * (1:n))*log(1 - theop)) / n
  m <- .6655
  l <- 1.6
  # Approximated p-value
  p.value  <- 1-(pnorm(sqrt(l / ad) * (ad / m - 1)) + exp( 2 * l / m ) * pnorm( -sqrt( l / ad) * ( ad / m + 1)))
  results <- list(statistic = c("AD" = ad), p.value = p.value)                  
  class(results) = "htest"
  return(results)  
}
