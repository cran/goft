
# Transformation test for the lognormal distribution
lnorm_test <- function(x){
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x) & length(x) <= 1) stop(paste(DNAME, "must be a numeric vector containing more than 1 observation"))
  x <- x[!is.na(x)]
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- as.vector(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if (min(x) < 0) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  result  <- shapiro.test(log(x))
  method  <-  "Test for the lognormal distribution based on a transformation to normality"
  alternative <-  paste(DNAME," does not follow a Lognormal distribution.")
  results  <- list( p.value = result$p.value, data.name = DNAME,
         method = method) 
  class(results) = "htest"
  return(results)
}


