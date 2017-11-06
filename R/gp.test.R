# Bootstrap goodness-ofi-fit test for the Generalized Pareto distribution

gp_test <- function(x, B = 999){
  dname <- deparse(substitute(x))
  if (!is.numeric(x)) stop(paste(dname, "must be a numeric vector"))
  if (sum(is.na(x)) > 0) warning("NA values have been deleted")
  x <- x[!is.na(x)]
  x <- as.vector(x)
  n <- length(x)   #  sample size without NA values
  if (n <= 1) stop("sample size must be larger than 1")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) stop("all observations are identical")
  if (min(x) < 0) stop("There are negative observations. \nAll data must be positive real numbers.")
  gammap <- .amle_method(x, k = ceiling(.2 * n))[1]    
  gamman <- .combined_method(x)[1]
  r1 <- .R1(x)     # observed value of R^-
  r2 <- .R2(x)     # observed value of R^+
  p.value1 <- sum(replicate(B, .R1(.rgp(n, shape = gamman))) < r1) / B  # bootstrap p-value for H_0^- 
  p.value2 <- sum(replicate(B, .R2(.rgp(n, shape = gammap))) < r2) / B  # bootstrap p-value for H_0^+ 
  p.value  <- max(p.value1, p.value2)    # p-value of the intersection-union test
  pvalues <- matrix(round(c(c(p.value1, p.value2),c(r1, r2)),4), ncol = 2)
  colnames(pvalues) <- c("p.value","R")
  name1  <- paste("H_0^-: ", dname, " follows a gPd with NEGATIVE shape parameter", sep="")
  name2  <- paste("H_0^+: ", dname, " follows a gPd with POSITIVE shape parameter", sep="")
  rownames(pvalues) <- c(name1, name2) 
  results  <- list("p.value" = p.value, method = "Bootstrap test of fit for the generalized Pareto distribution", data.name=dname, pvalues = pvalues)
  class(results) = "htest"
  return(results)
}


# INTERNAL FUNCTIONS

# Asymptotic maximum likelihood estimators 

.amle_method <- function(x, k){
  x  <- sort(x)
  n  <- length(x)
  nk <- n - k
  x1 <- x[(nk+1):n]
  w  <- log(x1)
  g  <-  - (w[1] - sum(w) / k)        # equation (4)
  sigma <- g * exp(w[1] + g * log(k / n))
  return(c(g, sigma))
}

# Combined estimators 
.combined_method <- function(x){
  m     <- mean(x)
  maxi  <- max(x)
  g     <- m / (m - maxi)    # equation (7) 
  sigma <- - g * maxi        # equation (6)
  return(c(g, sigma))
}

# Test statistic for H_0^-

.R1 <- function(x){
  gamma_neg <- .combined_method(x)[1]
  Fn        <- ecdf(x)
  x1        <- x[x != max(x)]
  z1        <- (1 - Fn(x1))^( - gamma_neg) 
  return(abs(cor(x1, z1)))
}

# Test statistic for H_0^+
.R2  <- function(x){
  n              <- length(x)
  Fn             <- ecdf(x)
  gamma_positive <- .amle_method(x, ceiling(.2 * n))[1]
  x1             <- x[x != max(x)]
  y1             <- (1 - Fn(x1))^( - gamma_positive) 
  x.star         <- log(x1)
  y.star         <- log( y1 -1 )
  if (gamma_positive <= 0.5)	return(cor(x1, y1))  
  if (gamma_positive >  0.5)  return((cor(x.star, y.star)))
}

# Simulation of random numbers from the gPd
.rgp  <- function (n, shape){
  if (shape != 0) 
    return((1 / shape) * (runif(n)^(-shape) - 1))
  else return(rexp(n, 1))
}

