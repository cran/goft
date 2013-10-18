gumbel.test = function(x)
{  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  n <- length(x)
  if (n < 20 || n > 250) 
    stop("sample size must be between 20 and 250")
  samplerange <- max(x) - min(x)
  if (samplerange == 0) 
    stop("all observations are identical")
  
  z   <- combn(x,2,max)  
  nz  <- length(z)
  Fn  <- rep(NA,nz)
  Fn  <- ecdf(z)
  FnZ <- Fn(z)
  y   <- rep(NA,nz)
  y   <- log(qweibull(FnZ,1,lower.tail=TRUE))   # y es la inversa de la fn. de distn. Gumbel evaluada en FnZ(z) 
 #y   <- qnorm(FnZ,lower.tail=TRUE)
  r   <- NA
  logic  <- (z!= max(z))
  r      <- cor(z[logic],y[logic],method="pearson")
  
  median_n = -3.02921 -.03846*n +.00023427*n**2-.000000471091*n**3
  if(n <= 60 )
  {
    s_n = 0.7588 - 0.01697*n + 0.000399*n**2 - 0.000003*n**3
  }
  if( n >60 ) s_n = .53
  
  rl = log(1-r)
  p.value = pnorm(rl,mean=median_n,sd=s_n,lower.tail=FALSE)
  results <- list(statistic = c(R = r), p.value = p.value, 
                  method = "Test of fit for the Gumbel distribution", 
                  data.name = DNAME)
  class(results) = "htest"
  return(results)
}

