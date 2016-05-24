
# Test for the Weibull distribution based on a transformation to Gumbel observations

weibull.test <- function(x, method = "transf", N = 1000){
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x) & length(x) <= 1) stop(paste(DNAME, "must be a numeric vector containing more than 1 observation"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  n <- length(x)   # adjusted sample size without NA values
  if (min(x) < 0) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  result  <- ev.test(-log(x), dist = "gumbel", method = "ratio", N = N)
  method  <-  "Test for the Weibull distribution"
  alternative <-  paste(DNAME," does not follow a Weibull distribution.")
  results  <- list( p.value = result$p.value, data.name = DNAME,
                    method = method) 
  class(results) = "htest"
  return(results)
  
}