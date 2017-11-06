# A generalization of Shapiro-Wilk test for multivariate normality 
# Villasenor-Alva and Gonzalez-Estrada (2009).

mvshapiro_test <- function(X){
	dname <- deparse(substitute(X))
	if (is.vector(X) == TRUE) X <- cbind(X)
	stopifnot(is.matrix(X))
	n     <-  nrow(X)
	p     <-  ncol(X)
  if (n < 12 || n > 5000) stop("Sample size must be between 12 and 5000.")
	if (n <= p) 	stop("Sample size must be larger than vector dimension.")
	if (n > p)	{
    x       <- scale(X, scale = FALSE)
    eigenv  <- eigen(var(X), symmetric = TRUE)
	  e.vec   <- as.matrix(eigenv$vectors)
	  sqrS	  <- e.vec %*% diag(1 / sqrt(eigenv$values), ncol = p) %*% t(e.vec)
	  z       <- t(sqrS %*% t(x))
		w       <- rep(NA, p)
		w       <- apply(z, 2, .shapiro)
		wast    <- mean(w)
		y       <- log(n)
		w1      <- log(1 - wast)
		m       <- -1.5861 - 0.31082 * y - 0.083751 * y^2 + 0.0038915 * y^3
		s       <- exp(-0.4803 - 0.082676 * y + 0.0030302 * y^2)
		s2      <- s^2
	  sigma2  <- log((p - 1 + exp(s2)) / p)
		mu1     <- m + s2 / 2 - sigma2 / 2
		p.value <- pnorm(w1, mean = mu1, sd = sqrt(sigma2), lower.tail = FALSE) 
		results <- list(statistic = c(MVW = wast), "p.value" = p.value, 
                    method = "Generalized Shapiro-Wilk test for Multivariate Normality", 
                    data.name = dname)
		class(results) = "htest"
		return(results)
	}
}

# Internal function  
.shapiro <- function(x) shapiro.test(x)$statistic
