# Parameter estimators for the generalized Pareto distribution

gp.fit <-  function(x, method){
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
  if (!(method %in% c("amle", "combined")))  stop("Unknown method. Please check the 'method' argument in the help files.")
  fit <- switch(method, 'amle' = .amle_method(x, ceiling(.2 * n)), 'combined' = .combined_method(x))
  fit <- as.matrix(round(fit,4), ncol = 1)
  colnames(fit) <- c("Parameter estimates")
  rownames(fit) <- c("shape", "scale")
  return(fit)
}
