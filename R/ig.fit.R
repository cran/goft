# MLE for the Inverse Gaussian distribution

ig_fit <- function(x){
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if (min(x) < 0 ) stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  n <- length(x)
  mu.hat <- mean(x)
  lambda.hat <- n / (sum(1 / x - 1 / mu.hat))
  fit <- as.matrix(c(mu.hat, lambda.hat))
  colnames(fit) <- c("Inverse Gaussian MLE ")
  rownames(fit) <- c("mu", "lambda")
  return(fit)
}