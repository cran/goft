# Correlation test for univariate normality
normal_test <- function(x, method = "cor"){  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted.")
  x <- x[is.na(x) == FALSE]
  n <- length(x)   # adjusted sample size without NA values
  if (n < 10 || n > 400) stop("sample size must be between 10 and 400")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  z   <- combn(x, 2, sum)  
  nz  <- length(z)
  Fn  <- rep(NA, nz)
  Fn  <- ecdf(z)
  FnZ <- Fn(z)
  y   <- rep(NA, nz)
  y   <- qnorm(FnZ, lower.tail = TRUE)
  r   <- NA
  logic  <- (z != max(z))
  r      <- cor(z[logic], y[logic], method = "pearson")
  s_n <- NA
  if(n <= 30 && n >= 10) s_n = 0.3994 + 0.6394*log10(n) - 0.2033*log10(n)**2 
  if(n <= 400 && n >= 31) s_n = 0.89994   
  median_n <- - 1.928 - 3.553 * log10(n) + 0.7265 * log10(n)**2 - 0.1243 * log10(n)**3
  rl <- log(1 - r)
  p.value <- pnorm(rl, mean = median_n, sd = s_n, lower.tail = FALSE)
  results <- list(statistic = c(R = r), p.value = p.value, 
                  method = "Correlation test for normality ", 
                  data.name = DNAME, alternative = paste(DNAME, "does not follow a normal distribution."))
  class(results) <- "htest"
  return(results)
}


