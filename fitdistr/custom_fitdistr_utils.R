library("modi")
start.arg.default <- function(x, distr, weights=NULL)
{
  if (distr == "norm") {
    n <- length(x)
    sd0 <- sqrt((n - 1)/n) * sd(x)
    mx <- mean(x)
    start <- list(mean=mx, sd=sd0)
  }else if (distr == "lnorm") {
    if (any(x <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    start <- list(meanlog=ml, sdlog=sd0)
  }else if (distr == "pois") {
    start <- list(lambda=mean(x))
  }else if (distr == "exp") {
    if (any(x < 0)) 
      stop("values must be positive to fit an exponential  distribution")
    start <- list(rate=1/mean(x))
  }else if (distr == "gamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(shape=m^2/v, rate=m/v)
  }else if (distr == "nbinom") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    size <- ifelse(v > m, m^2/(v - m), 100)
    start <- list(size = size, mu = m) 
  }else if (distr == "geom" ) {
    m <- mean(x)
    prob <- ifelse(m>0, 1/(1+m), 1)
    start <- list(prob=prob)        
  }else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    if (is.null(weights)) {m <- mean(x)} else {m <- weighted.mean(x, weights)}
    if (is.null(weights)) {v <- var(x)} else {v <- weighted.var(x, weights)}
    v <- (n - 1)/n*v
    aux <- m*(1-m)/v - 1
    start <- list(shape1=m*aux, shape2=(1-m)*aux)
  }else if (distr == "weibull") {
    if (any(x < 0)) 
      stop("values must be positive to fit an Weibull  distribution")
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    start <- list(shape = shape, scale = scale)
  }else if (distr == "logis") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(location=m, scale=sqrt(3*v)/pi)
  }else if (distr == "cauchy") {
    start <- list(location=median(x), scale=IQR(x)/2)
  }else if (distr == "unif"){
    start <- list(min=0, max=1)
  }else if (distr == "invgamma")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2*m2-m1^2)/(m2-m1^2)
    scale <- m1*m2/(m2-m1^2)
    start <- list(shape=shape, scale=scale)
  }else if (distr == "llogis")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- 2*log(3)/(log(p75)-log(p25))
    scale <- exp(log(p75)+log(p25))/2
    start <- list(shape=shape, scale=scale)
  }else if (distr == "invweibull")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g*log(p75)-log(p25))/(g-1))
    scale <-log(log(4))/(log(shape)-log(p25))
    start <- list(shape=shape, scale=max(scale, 1e-9))
  }else if (distr == "pareto1")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    #http://www.math.umt.edu/gideon/pareto.pdf
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n*m1-x1)/(n*(m1-x1))
    min <- x1*(n*shape - 1)/(n*shape)
    start <- list(shape=shape, min=min)
  }else if (distr == "pareto")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1*m2)/(m2-2*m1^2)
    shape <- 2*(m2-m1^2)/(m2-2*m1^2)
    start <- list(shape=shape, scale=scale)
  }else if (distr == "lgamma")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-gamma distribution")
    #p228 of Klugmann and Hogg (1984)
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2-m1^2)
    lambda <- m1/(m2-m1^2)
    start <- list(shapelog=alpha, ratelog=lambda)
  }else if (distr == "trgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(shape1=m^2/v, shape2=1, rate=m/v)
  }else if (distr == "invtrgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(1/x)
    m <- mean(1/x)
    v <- (n - 1)/n*var(1/x)
    start <- list(shape1=m^2/v, shape2=1, rate=m/v)
  }else
    stop(paste0("Unknown starting values for distribution ", distr, "."))
  
  return(start)
} 

manageparam <- function(start.arg, fix.arg, obs, distname, weights=NULL)
{
  #if clause with 3 different cases:
  #start.arg : NULL | named list | a function
  
  if(is.null(start.arg))
  {
    trystart <- try(start.arg.default(obs, distname, weights), silent = TRUE)
    if(class(trystart) == "try-error")
    {
      cat("Error in computing default starting values.\n")
      stop(trystart)
    }
    lstart <- trystart
    #lstart should be a named list but check it
    hasnoname <- is.null(names(lstart)) || !is.list(lstart)
    if(hasnoname)
      stop("Starting values must be a named list, error in default starting value.")
    
  }else if(is.list(start.arg))
  {
    hasnoname <- is.null(names(start.arg))
    if(hasnoname)
      stop("Starting values must be a named list (or a function returning a named list).")
    lstart <- start.arg
  }else if(is.function(start.arg))
  {
    trystart <- try(start.arg(obs), silent = TRUE)
    if(class(trystart) == "try-error")
    {
      cat("Error in computing starting values with your function.\n")
      stop(trystart)
    }
    lstart <- trystart
    hasnoname <- is.null(names(lstart)) || !is.list(lstart)
    if(hasnoname)
      stop("Starting values must be a named list, your function does not return that.")
  }else
    stop("Wrong type of argument for start")
  
  #if clause with 3 different cases:
  #fix.arg : NULL | named list | a function
  if(is.null(fix.arg))
  {
    lfix <- NULL
  }else if(is.list(fix.arg))
  {
    hasnoname <- is.null(names(fix.arg))
    if(hasnoname)
      stop("Fixed parameter values must be a named list (or a function returning a named list).")
    lfix <- fix.arg
  }else if(is.function(fix.arg))
  {
    tryfix <- try(fix.arg(obs), silent = TRUE)
    if(class(tryfix) == "try-error")
    {
      cat("Error in computing fixed parameter values with your function.\n")
      stop(tryfix)
    }
    lfix <- tryfix
    hasnoname <- is.null(names(lfix)) || !is.list(lfix)
    if(hasnoname)
      stop("Fixed parameter values must be a named list, your function does not return that.")
  }else
    stop("Wrong type of argument for fix.arg")
  
  #eliminate arguments both in lstart and lfix (when start.arg was NULL)
  if(is.null(start.arg) && !is.null(lfix))
  {
    lstart <- lstart[!names(lstart) %in% names(lfix)]
    if(length(lstart) == 0)
      stop("Don't need to use fitdist() if all parameters have fixed values")
  }
  
  list("start.arg"=lstart, "fix.arg"=lfix)
}

computegetparam <- function(argdistname)
{
  #remove first argument, that should be "x", "p", "q", or "n", see ?dgamma, pgamma, qgamma
  argdistname <- argdistname[-1]
  nonparaminR <- c("x", "p", "q", "n") #defensive programming
  #remove other arguments, see ?dgamma, pgamma, qgamma, dbeta
  nonparaminR <- c(nonparaminR, "log", "log.p", "lower.tail", "ncp")
  nonparaminActuar <- c("limit", "order", "t")
  nonparaminGamlssdist <- "fast"
  nonparamspecial <- c("...", "..1", "..2")
  #see ?dnig, dhyperb, dskewlap, dgig,...
  nonparaminGenHyperbolic <- c("param", "KOmega", "ibfTol", "nmax", "method", "intTol",
                               "valueOnly", "nInterpol", "uniTol", "subdivisions", "logPars")
  #see ?dsn
  nonparamsn <- "dp"
  
  plist <- setdiff(argdistname, nonparaminR)
  plist <- setdiff(plist, nonparaminActuar)
  plist <- setdiff(plist, nonparaminGamlssdist)
  plist <- setdiff(plist, nonparamspecial)
  plist <- setdiff(plist, nonparaminGenHyperbolic)
  plist <- setdiff(plist, nonparamsn)
  
  plist
}

checkparamlist <- function(start.arg, fix.arg, argdistname, hasnodefaultval)
{
  errtxt <- list(t3="'start' must specify names which are arguments to 'distr'.",
                 t4="'fix.arg' must specify names which are arguments to 'distr'.",
                 t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.",
                 t6="'start' should not have NA or NaN values.",
                 t7="'fix.arg' should not have NA or NaN values.",
                 t8="Some parameter names have no starting/fixed value.",
                 t9="Some parameter names have no starting/fixed value but have a default value.")
  
  vstart <- unlist(start.arg)
  #check unexpected names
  m <- match(names(vstart), argdistname)
  if (any(is.na(m))) 
    stop(errtxt$t3)
  #check NA/NaN values
  cat(vstart)
  if(any(is.na(vstart) || is.nan(vstart)))
    stop(errtxt$t6)
  if(!is.null(fix.arg))
  {
    vfix <- unlist(fix.arg)
    #check unexpected names
    mfix <- match(names(vfix), argdistname)
    if (any(is.na(mfix))) 
      stop(errtxt$t4)
    
    # check that some parameters are not both in fix.arg and start
    minter <- match(names(vstart), names(vfix))
    if (any(!is.na(minter)))
      stop(errtxt$t5)
    
    #check NA/NaN values
    if(any(is.na(vfix) || is.nan(vfix)))
      stop(errtxt$t7)
    allparname <- names(c(vstart, vfix))
  }else
    allparname <- names(vstart)
  
  theoparam <- computegetparam(argdistname)
  #special case where both scale and rate are allowed, see ?dgamma
  if("scale" %in% theoparam && "rate" %in% theoparam)
  {
    errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
    #special case where both prob and mu are allowed, see ?dnbinom
  }else if(length(theoparam) == 3 && all(c("size", "prob", "mu") %in% theoparam))
  {
    errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
  }else
    errt8 <- any(!theoparam %in% allparname)
  #only make a warning if unset arguments have a default value
  if(errt8)
  {
    unsetarg <- theoparam[!theoparam %in% allparname] 
    if(any(hasnodefaultval[unsetarg]))
      stop(errtxt$t8)
    else
      warning(errtxt$t9)
  }
  
  list("start.arg"=start.arg, "fix.arg"=fix.arg)
}

testdpqfun <- function(distr, fun=c("d","p","q"), start.arg, 
                       fix.arg=NULL, discrete=FALSE)
{
  stopifnot(all(is.character(fun)))
  fun <- fun[fun %in% c("d","p","q")]
  stopifnot(length(fun) > 0)
  
  if(is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  if(is.function(fix.arg))
    stop("fix.arg should be either a named list or NULL but not a function")
  
  op <- options() #get current options
  #print(getOption("warn"))
  options(warn=-1)
  
  res <- NULL
  if("d" %in% fun)
    res <- rbind(res, test1fun(paste0("d", distr), start.arg, fix.arg))
  if("p" %in% fun)
    res <- rbind(res, test1fun(paste0("p", distr), start.arg, fix.arg))
  if("q" %in% fun)
    res <- rbind(res, test1fun(paste0("q", distr), start.arg, fix.arg))
  
  options(op)     # reset (all) initial options
  res
} 

test1fun <- function(fn, start.arg, fix.arg, dpqr)
{
  res <- data.frame(ok=FALSE, txt="")
  stopifnot(is.list(start.arg))
  if(!is.null(fix.arg))
    stopifnot(is.list(fix.arg))
  
  #does the function exist?
  if(!exists(fn, mode="function"))
  {
    res$txt <- paste("The", fn, "function must be defined")
    return(res)
  }
  
  #naming convention
  if(missing(dpqr))
    dpqr <- substr(fn, 1, 1)
  firstarg_theo <- switch(dpqr, "d"="x", "p"="q", "q"="p", "r"="n")
  firstarg_found <- names(formals(fn))[1]
  if(firstarg_found != firstarg_theo)
  {    
    t0 <- paste("The", fn, "function should have its first argument named:", firstarg_theo)
    res$txt <- paste(t0, "as in base R")
    return(res)
  }
  
  #zero-component vector
  res0 <- try(do.call(fn, c(list(numeric(0)), start.arg, fix.arg)), silent=TRUE)
  t0 <- paste("The", fn, "function should return a zero-length vector when input has length zero and not raise an error")
  t1 <- paste("The", fn, "function should return a zero-length vector when input has length zero")
  if(class(res0) == "try-error")
  {
    res$txt <- t0
    return(res)
  }
  if(length(res0) != 0)
  {
    res$txt <- t1
    return(res)
  }
  
  #inconsistent value
  x <- c(0, 1, Inf, NaN, -1)
  res1 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t2 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent values and not raise an error")
  if(class(res1) == "try-error")
  {
    res$txt <- t2
    return(res)
  }
  
  #missing value 
  x <- c(0, 1, NA)
  res2 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t4 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not raise an error")
  t5 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not remove missing values")
  if(class(res2) == "try-error")
  {
    res$txt <- t4
    return(res)
  }
  if(length(res2) != length(x))
  {
    res$txt <- t5
    return(res)
  }
  
  #inconsistent parameter
  x <- 0:1
  start.arg <- lapply(start.arg, function(x) -x)
  res3 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t6 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent parameters and not raise an error")
  if(class(res3) == "try-error")
  {
    res$txt <- t6
    return(res)
  }
  
  #wrong parameter name
  x <- 0:1
  names(start.arg) <- paste0(names(start.arg), "_")
  res4 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t8 <- paste("The", fn, "function should raise an error when names are incorrectly named")
  if(class(res4) != "try-error")
  {
    res$txt <- t8
    return(res)
  }  
  return(data.frame(ok=TRUE, txt=""))
}
