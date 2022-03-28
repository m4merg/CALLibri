source("custom_fitdistr_utils.R")

mgedist <- function (data, distr, gof = "CvM", start = NULL, fix.arg = NULL, 
                     optim.method = "default", lower = -Inf, upper = Inf, custom.optim = NULL, 
                     silent = TRUE, gradient = NULL, checkstartfix = FALSE, weights, ...) 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else distname <- distr
  if (is.element(distname, c("binom", "nbinom", "geom", "hyper", 
                             "pois"))) 
    stop("Maximum goodness-of-fit estimation method is not intended to fit discrete distributions")
  pdistname <- paste("p", distname, sep = "")
  if (!exists(pdistname, mode = "function")) 
    stop(paste("The ", pdistname, " function must be defined"))
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  argddistname <- names(formals(ddistname))
  if (is.null(custom.optim)) 
    optim.method <- match.arg(optim.method, c("default", 
                                              "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                                              "Brent"))
  gof <- match.arg(gof, c("CvM", "KS", "AD", "ADR", "ADL", 
                          "AD2R", "AD2L", "AD2", "ADW"))
  start.arg <- start
  if (is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  my3dots <- list(...)
  if ("weights" %in% names(my3dots)) 
    stop("Weights is not allowed for maximum GOF estimation")
  if (is.vector(data)) {
    cens <- FALSE
    if (!(is.numeric(data) & length(data) > 1)) 
      stop("data must be a numeric vector of length greater than 1 for non censored data\n            or a dataframe with two columns named left and right and more than one line for censored data")
  }
  else {
    cens <- TRUE
    censdata <- data
    if (!(is.vector(censdata$left) & is.vector(censdata$right) & 
          length(censdata[, 1]) > 1)) 
      stop("data must be a numeric vector of length greater than 1 for non censored data\n        or a dataframe with two columns named left and right and more than one line for censored data")
    pdistname <- paste("p", distname, sep = "")
    if (!exists(pdistname, mode = "function")) 
      stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))
  }
  if (cens) {
    lcens <- censdata[is.na(censdata$left), ]$right
    if (any(is.na(lcens))) 
      stop("An observation cannot be both right and left censored, coded with two NA values")
    rcens <- censdata[is.na(censdata$right), ]$left
    ncens <- censdata[censdata$left == censdata$right & 
                        !is.na(censdata$left) & !is.na(censdata$right), 
                      ]$left
    icens <- censdata[censdata$left != censdata$right & 
                        !is.na(censdata$left) & !is.na(censdata$right), 
                      ]
    data <- c(rcens, lcens, ncens, (icens$left + icens$right)/2)
  }
  #if (DEBUG) {cat("START\n")}
  #print(start)
  if (!checkstartfix) {
    arg_startfix <- manageparam(start.arg = start, fix.arg = fix.arg, 
                                obs = data, distname = distname)
    hasnodefaultval <- sapply(formals(ddistname), is.name)
    arg_startfix <- checkparamlist(arg_startfix$start.arg, 
                                   arg_startfix$fix.arg, argddistname, hasnodefaultval)
    if (is.function(fix.arg)) 
      fix.arg.fun <- fix.arg
    else fix.arg.fun <- NULL
  }
  else {
    arg_startfix <- list(start.arg = start, fix.arg = fix.arg)
    fix.arg.fun <- NULL
  }
  vstart <- unlist(arg_startfix$start.arg)
  if (is.null(vstart)) 
    stop("Starting values could not be NULL with checkstartfix=TRUE")
  fix.arg <- arg_startfix$fix.arg
  if (!cens) {
    if (gof == "CvM") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        1/(12 * n) + sum((theop - (2 * 1:n - 1)/(2 * 
                                                   n))^2)
      }
    else if (gof == "KS") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        obspu <- seq(1, n)/n
        obspl <- seq(0, n - 1)/n
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        max(pmax(abs(theop - obspu), abs(theop - obspl)))
      }
    else if (gof == "AD")
      fnobj <- function(par, fix.arg, obs, pdistnam, DEBUG = FALSE)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
        
        err <- - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
	sum <- 0
	for(i in 1:n) {
		i1 <- i + 0
		n1 <- n + 0
		element <- (2*i1 - 1)*log(theop[i])/(2*n1) + (1-(2*i1-1)/(2*n1))*log(1-theop[i])
		sum <- sum + element
		if (DEBUG) {
			#cat(sum,":",element," ")
			}
		}
	if(DEBUG) {
		#cat("\n")
		}
	err1 <- -n - 2*sum
        if (DEBUG) {
	  #print(par)
	  #print(theop)
          #print(err)
	  #print(err1)
        }
        err1
      }
    else if (gof == "ADW") 
      fnobj <- function(par, fix.arg, obs, pdistnam, DEBUG = FALSE) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        #err <- - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) )
	sum <- 0
	sumw <- 0
	for(i in 1:n) {
		sumw <- sumw + weights[i]
		if((!is.na(theop[i]))&(theop[i] > 0.9999999999999999)) {theop[i] <- 0.9999999999999999}
		element <- (2*n*sumw - 1)*log(theop[i])/(2*n) + (1-(2*n*sumw-1)/(2*n))*log(1-theop[i])
		sum <- sum + weights[i]*element
		if ((DEBUG)&(i > 68)) {
			cat(sum,":",weights[i],":",element,"|","\n")
			}
		}
	if(DEBUG) {
		#cat("\n")
		}
	err <- -n - 2*sum
        if (DEBUG) {
	  #print(par)
	  #print(tail(theop))
          #print(err)
        }
        err
      }
    else if (gof == "ADR") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        n/2 - 2 * sum(theop) - mean((2 * 1:n - 1) * 
                                      log(1 - rev(theop)))
      }
    else if (gof == "ADL") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        -3 * n/2 + 2 * sum(theop) - mean((2 * 1:n - 
                                            1) * log(theop))
      }
    else if (gof == "AD2R") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        2 * sum(log(1 - theop)) + mean((2 * 1:n - 1)/(1 - 
                                                        rev(theop)))
      }
    else if (gof == "AD2L") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        2 * sum(log(theop)) + mean((2 * 1:n - 1)/theop)
      }
    else if (gof == "AD2") 
      fnobj <- function(par, fix.arg, obs, pdistnam) {
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), 
                                     as.list(fix.arg)))
        2 * sum(log(theop) + log(1 - theop)) + mean(((2 * 
                                                        1:n - 1)/theop) + ((2 * 1:n - 1)/(1 - rev(theop))))
      }
  }
  else stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
  loglik <- function(par, fix.arg, obs, ddistnam) {
    sum(log(do.call(ddistnam, c(list(obs), as.list(par), 
                                as.list(fix.arg)))))
  }
  owarn <- getOption("warn")
  if (is.null(custom.optim)) {
    hasbound <- any(is.finite(lower) | is.finite(upper))
    if (optim.method == "default") {
      meth <- ifelse(length(vstart) > 1, "Nelder-Mead", 
                     "BFGS")
    }
    else meth <- optim.method
    if (meth == "BFGS" && hasbound && is.null(gradient)) {
      meth <- "L-BFGS-B"
      txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
      txt2 <- "The method is changed to L-BFGS-B."
      warning(paste(txt1, txt2))
    }
    options(warn = ifelse(silent, -1, 0))
    if (hasbound) {
      if (!is.null(gradient)) {
        opt.fun <- "constrOptim"
      }
      else {
        if (meth == "Nelder-Mead") 
          opt.fun <- "constrOptim"
        else if (meth %in% c("L-BFGS-B", "Brent")) 
          opt.fun <- "optim"
        else {
          txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
          txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
          stop(paste(txt1, txt2))
        }
      }
      if (opt.fun == "constrOptim") {
        npar <- length(vstart)
        lower <- as.double(rep_len(lower, npar))
        upper <- as.double(rep_len(upper, npar))
        haslow <- is.finite(lower)
        Mat <- diag(npar)[haslow, ]
        hasupp <- is.finite(upper)
        Mat <- rbind(Mat, -diag(npar)[hasupp, ])
        colnames(Mat) <- names(vstart)
        rownames(Mat) <- paste0("constr", 1:NROW(Mat))
        Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
        names(Bnd) <- paste0("constr", 1:length(Bnd))
        initconstr <- Mat %*% vstart - Bnd
        if (any(initconstr < 0)) 
          stop("Starting values must be in the feasible region.")
        opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                              f = fnobj, ui = Mat, ci = Bnd, grad = gradient, 
                                              fix.arg = fix.arg, obs = data, pdistnam = pdistname, 
                                              hessian = !is.null(gradient), method = meth, 
                                              ...), silent = TRUE)
        if (!inherits(opttryerror, "try-error")) 
          if (length(opt$counts) == 1) 
            opt$counts <- c(opt$counts, NA)
      }
      else {
        opttryerror <- try(opt <- optim(par = vstart, 
                                        fn = fnobj, fix.arg = fix.arg, obs = data, 
                                        gr = gradient, pdistnam = pdistname, hessian = TRUE, 
                                        method = meth, lower = lower, upper = upper, 
                                        ...), silent = TRUE)
      }
    }
    else {
      opt.fun <- "optim"
      opttryerror <- try(opt <- optim(par = vstart, fn = fnobj, 
                                      fix.arg = fix.arg, obs = data, gr = gradient, 
                                      pdistnam = pdistname, hessian = TRUE, method = meth, 
                                      lower = lower, upper = upper, ...), silent = TRUE)
    }
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The function optim encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        #print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, loglik = NA, hessian = NA))
    }
    if (opt$convergence > 0) {
      warnings("The function optim failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, 
                value = opt$value, hessian = opt$hessian, optim.function = opt.fun, 
                optim.method = meth, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                weights = NULL, counts = opt$counts, optim.message = opt$message, 
                loglik = loglik(opt$par, fix.arg, data, ddistname), 
                gof = gof)
  }
  else {
    options(warn = ifelse(silent, -1, 0))
    if (!cens) 
      opttryerror <- try(opt <- custom.optim(fn = fnobj, 
                                             fix.arg = fix.arg, obs = data, pdistnam = pdistname, 
                                             par = vstart, ...), silent = TRUE)
    else stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The customized optimization function encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        #print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, value = NA, hessian = NA))
    }
    if (opt$convergence > 0) {
      warnings("The customized optimization function failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    argdot <- list(...)
    method.cust <- argdot$method
    res <- list(estimate = opt$par, convergence = opt$convergence, 
                value = opt$value, hessian = opt$hessian, optim.function = custom.optim, 
                optim.method = method.cust, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                weights = NULL, counts = opt$counts, optim.message = opt$message, 
                loglik = loglik(opt$par, fix.arg, data, ddistname), 
                gof = gof)
  }
  
  #pars <- list(shape1 = 2.0, shape2 = 1000.0)
  #print(fnobj(pars, fix.arg, data, 'pbeta', DEBUG = TRUE))
  return(res)
}


fitdistC <- function (data, distr, method = c("mge"), start=NULL, 
                     fix.arg=NULL, discrete, keepdata = TRUE, keepdata.nb=100, weights = NULL, ...) 
{
  #print(weights)
  if (!is.null(weights)) {
	  weights <- weights/sum(weights)
  	}
  #check argument distr
  if (!is.character(distr)) 
    distname <- substring(as.character(match.call()$distr), 2)
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  #pdistname <- paste("p", distname, sep="")
  #if (!exists(pdistname, mode="function"))
  #    stop(paste("The ", pdistname, " function must be defined"))
  #check argument discrete
  if(missing(discrete))
  {
    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
      discrete <- TRUE
    else
      discrete <- FALSE
  }
  if(!is.logical(discrete))
    stop("wrong argument 'discrete'.")
  if(!is.logical(keepdata) || !is.numeric(keepdata.nb) || keepdata.nb < 2)
    stop("wrong arguments 'keepdata' and 'keepdata.nb'")
  #check argument method
  if(any(method == "mom"))
    warning("the name \"mom\" for matching moments is NO MORE used and is replaced by \"mme\"")
  method <- match.arg(method, c("mle", "mme", "qme", "mge"))
  if(method %in% c("mle", "mme", "mge"))
    dpq2test <- c("d", "p")
  else
    dpq2test <- c("d", "p", "q")
  #check argument data
  if (!(is.vector(data) & is.numeric(data) & length(data)>1))
    stop("data must be a numeric vector of length greater than 1")
  #encapsulate three dots arguments
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  n <- length(data)
  # manage starting/fixed values: may raise errors or return two named list
  arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data, 
                              distname=distname, weights=weights)
  #check inconsistent parameters
  argddistname <- names(formals(ddistname))
  hasnodefaultval <- sapply(formals(ddistname), is.name)
  arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg, 
                                 argddistname, hasnodefaultval)
  #arg_startfix contains two names list (no longer NULL nor function)
  #store fix.arg.fun if supplied by the user
  if(is.function(fix.arg))
    fix.arg.fun <- fix.arg
  else
    fix.arg.fun <- NULL
  # check d, p, q, functions of distname
  resdpq <- testdpqfun(distname, dpq2test, start.arg=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, discrete=discrete)
  if(any(!resdpq$ok))
  {
    for(x in resdpq[!resdpq$ok, "txt"])
      warning(x)
  }
  
  
  # Fit with mledist, qmedist, mgedist or mmedist
  if (method == "mme")
  {
    stop("Method is unavailable")
  }else if (method == "mle")
  {
    stop("Method is unavailable")
  }else if (method == "qme")
  {
    stop("Method is unavailable")
  }else if (method == "mge")
  {
    if (!"gof" %in% names(my3dots))
      warning("maximum GOF estimation has a default 'gof' argument set to 'CvM'")    
    
    if (is.null(weights)) {
      weights <- rep(1, length(data))
    }
    
    mge <- mgedist(data, distname, start=arg_startfix$start.arg, 
                   fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, weights = weights, ...)
    
    estimate <- mge$estimate
    sd <- NULL
    loglik <- mge$loglik
    npar <- length(estimate)
    aic <- -2*loglik+2*npar
    bic <- -2*loglik+log(n)*npar
    correl <- varcovar <- NULL
    
    convergence <- mge$convergence
    fix.arg <- mge$fix.arg
    weights <- NULL
  }else
  {
    stop("match.arg() does not work correctly")
  }
  
  #needed for bootstrap
  if (!is.null(fix.arg)) 
    fix.arg <- as.list(fix.arg)
  
  if(keepdata)
  {
    reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl, 
                    vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=data,
                    distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                    dots = my3dots, convergence = convergence, discrete = discrete, 
                    weights = weights)
  }else #just keep a sample set of all observations
  {
    n2keep <- min(keepdata.nb, n)-2
    imin <- which.min(data)
    imax <- which.max(data)
    subdata <- data[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE)]
    subdata <- c(subdata, data[c(imin, imax)])
    
    reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl, 
                    vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=subdata,
                    distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                    dots = my3dots, convergence = convergence, discrete = discrete, 
                    weights = weights)  
  }
  
  return(structure(reslist, class = "fitdist"))
  
}
