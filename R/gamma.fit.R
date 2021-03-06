gamma_fit <-  function(x){
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x)) stop(paste(DNAME, "must be a numeric vector"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   # adjusted sample size without NA values
  if (n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if (min(x) < 0) stop("There are negative observations. \nAll data must be positive real numbers.")
  
    b.check <- cov(x, log(x))
    a.check <- mean(x) / b.check
    fit     <- as.matrix(c(a.check, b.check))
    colnames(fit) <- c("Parameter estimates")
    rownames(fit) <- c("shape", "scale")
    return(fit)
}
