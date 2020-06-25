# Transformation test for univariate skew-normal distributions
sn_test <- function(x, method = "transf"){  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted.")
  x <- x[is.na(x) == FALSE]
  n <- length(x)   # adjusted sample size without NA values
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  estim <- sn.mple(y = x, penalty = "Qpenalty", opt.method = "nlminb")$cp
  xi.hat <- cp2dp(estim, family = "SN")[1] 
  y      <- abs(x - xi.hat) * sign(rnorm(n))
  p.value <- shapiro.test(y)$p.value
  results <- list(p.value = p.value, 
                  method = "Shapiro-Wilk test for skew-normal distributions", 
                  data.name = DNAME, alternative = paste(DNAME, "does not follow a skew-normal distribution."))
  class(results) <- "htest"
  return(results)
}

