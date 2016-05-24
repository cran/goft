
# Test of fit for the gamma distribution based on the ratio of two variance estimators

gamma_test <- function(x){  
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x) & length(x) <= 1) stop(paste(DNAME, "must be a numeric vector containing more than 1 observation"))
  x <- x[!is.na(x)]
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- as.vector(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if ( min(x) < 0 ) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  n <- length(x)   # adjusted sample size without NA values
  z <-  log(x)
  x.bar   <- mean(x)
  s2.x    <- var(x)
  b.check <- cov(x, z)  
  a.check <- x.bar / b.check
  v       <- sqrt(n * a.check) * (s2.x / (x.bar * b.check) - 1)
  p.value <- 2 * pnorm(abs(v), mean = 0, sd = sqrt(2), lower.tail = FALSE)
  alternative = c( paste(DNAME, " does not follow a Gamma distribution."))
  results <- list(statistic = c("V" = v), p.value = p.value, data.name = DNAME,
                  method = "Test of fit for the Gamma distribution")
  class(results) = "htest"
  return(results)
}
