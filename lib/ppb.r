options(warn=1)
args <- commandArgs()
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "poolr",
  "scModels",
  "Rmpfr",
  "statmod"
  )

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
      )
    )
}

gq_altered =function (n, kind = "legendre", alpha = 0, beta = 0)
{
    n <- as.integer(n)
    if (n < 0L)
        stop("need non-negative number of nodes")
    if (n == 0L)
        return(list(nodes = numeric(0L), weights = numeric(0L)))
    kind <- match.arg(kind, c("legendre", "chebyshev1", "chebyshev2",
        "hermite", "jacobi", "laguerre"))
    i <- 1L:n
    i1 <- i[-n]
    switch(kind, legendre = {
        lnmuzero <- log(2)
        a <- rep_len(0, n)
        b <- i1/sqrt(4 * i1^2 - 1)
    }, chebyshev1 = {
        lnmuzero <- log(pi)
        a <- rep_len(0, n)
        b <- rep_len(0.5, n - 1L)
        b[1] <- sqrt(0.5)
    }, chebyshev2 = {
        lnmuzero <- log(pi/2)
        a <- rep_len(0, n)
        b <- rep_len(0.5, n - 1L)
    }, hermite = {
        lnmuzero <- log(pi)/2
        a <- rep_len(0, n)
        b <- sqrt(i1/2)
    }, jacobi = {
        ab <- alpha + beta
        lnmuzero <- (ab + 1) * log(2) + lgamma(alpha + 1) + lgamma(beta +
            1) - lgamma(ab + 2)
        a <- i
        a[1] <- (beta - alpha)/(ab + 2)
        i2 <- i[-1]
        abi <- ab + 2 * i2
        a[i2] <- (beta^2 - alpha^2)/(abi - 2)/abi
        b <- i1
        b[1] <- sqrt(4 * (alpha + 1) * (beta + 1)/(ab + 2)^2/(ab +
            3))
        i2 <- i1[-1]
        abi <- ab + 2 * i2
        b[i2] <- sqrt(4 * i2 * (i2 + alpha) * (i2 + beta) * (i2 +
            ab)/(abi^2 - 1)/abi^2)
    }, laguerre = {
        a <- 2 * i - 1 + alpha
        b <- sqrt(i1 * (i1 + alpha))
        lnmuzero <- lgamma(alpha + 1)
    })
    b <- c(b, 0)
    z <- rep_len(0, n)
    z[1] <- 1
    ierr <- 0L
    out <- .Fortran("gausq2", n, as.double(a), as.double(b),
        as.double(z), ierr, PACKAGE = "statmod")
    x <- out[[2]]
    w <- out[[4]]
    w <- exp(mpfr(lnmuzero + 2 * log(abs(w)), 333))
    list(nodes = x, weights = w)
}


pBP = function(x,alp,bet,lam1=1,lam2=1) {
  if (missing(bet)){
    par=alp;
    alp=par[1]; bet=par[2]; lam1=par[3]; lam2=1; if (length(par)> 3) lam2=par[4];
  }

  x=x/lam2
  ff = function(k,m){
    if (max(m) < 100000) res= dpois(k,m) else res = dnorm(k,m,sqrt(m))
    return(res)
  }
  w = gq_altered(10,'jacobi', alpha=bet-1, beta=alp-1)
  gs = sum(mpfr(w$weight, 333)*ff(mpfr(x, 333),m=mpfr(lam1, 333)*(mpfr(1,333)+mpfr(w$node,333))/mpfr(2,333)))
  prob = 1/beta(mpfr(alp, 333),mpfr(bet,333))*2^(-mpfr(alp, 333)-mpfr(bet,333)+mpfr(1,333))*mpfr(gs,333)   # get P(X=x)
  prob=prob
  return(prob)
}


file_in <- args[6]
file_out_total <- args[7]
file_out_detailed <- args[8]
data <- read.table(file = file_in, sep = '\t', header = FALSE, col.names = c('seed', 'AD', 'DP', 'alpha', 'beta', 'mean', 'panel_size'), stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "nA"))
data$pvalue <- 'NA'

round_custom <- function(value) {
        if (value < 1) {
                return(0)
                }
        return(round(value))
        }

get_pval <- function(AD, DP, alpha, beta, mean, panel_size) {
	if (is.na(AD) || is.na(DP) || is.na(mean) || is.na(alpha) || is.na(beta)) {
		return("NA")
		}
	if (AD == 'NA' || DP == 'NA' || mean == 'NA' || alpha == 'NA' || beta == 'NA') {
		return("NA")
		}

	pval_local <- 1
	if (DP == 0) {
		return("NA")
		} else {
		pval_local <- pBP(floor(round_custom(AD) - 1), alpha, beta, DP, 1)
		if (pval_local == 0){
			pval_local <- exp(mpfr(1200,333))
			}
		if ((AD/DP) < mean) {
			return(1)
			} else {
			pval_local <- min(mpfr(1,333), pval_local*panel_size)
			return(pval_local)
			}
		}
	return(pval_local)
	}

fun_c_mpfr <- function(x, df) {
 res <- mpfr(1,333)
 res <- res/mpfr((2**(df/2)),333)
 res <- res/mpfr(gamma(df/2), 333)
 res <- res*(x**mpfr(((df/2) - 1), 333))
 res <- res*(exp(-x/2))
 return(res)
}

pchisq_c_mpfr <- function(stat, df) {
        lower = mpfr(stat,333)
        upper <- mpfr(0,333)
        if (stat < 100) {
                upper <- mpfr(1000,333)
                } else {upper <- stat*mpfr(10,333)}
        suppressWarnings(a <- integrateR(fun_c_mpfr, lower, upper, mpfr(df, 333))$value)
        suppressWarnings(b <- integrateR(fun_c_mpfr, mpfr(0,333), upper, mpfr(df, 333))$value)
        res = (a/b)
        return(res)
        }

fisher_mpfr <- function(p) {
        k <- length(mpfr(p, 333))
        statistic <- mpfr(-2,333) * sum(log(mpfr(p, 333)))
        pval <- pchisq_c_mpfr(statistic, df = 2 * k)
        return(pval)
        }


phred_mpfr <- function(value) {
	if ((typeof(value) == 'character')&&(value == 'NA')) {
		return('NA')
		} else {
		output <- (-10)*log(mpfr(value, 333))/log(10)
		return(output)
		}
        }

unphred_mpfr <- function(value) {
	if ((typeof(value) == 'character')&&(value == 'NA')) {
		return('NA')
		} else {
		output <- mpfr(mpfr(10,333)**(-1*mpfr(value, 333)/mpfr(10, 333)),333)
		return(output)
		}
	}

mpfr_to_str <- function(value) {
	if ((typeof(value) == 'character')&&(value == 'NA')) {
		return('NA')
		} else {
		pval.output <- capture.output(value)[2]
		return(substr(pval.output,5,nchar(pval.output)))
		}
	}

make_mpfr_vector <- function(vector.inp) {
        res <- c(unphred_mpfr(vector.inp[1]))
        if (length(vector.inp) > 1) {
                for(i in 2:length(vector.inp)){
                        res <- append(res, unphred_mpfr(vector.inp[i]))
                        }
                }
        return(res)
        }


for (i in 1:nrow(data)) {
	pval <- get_pval(data[i,]$AD, data[i,]$DP, data[i,]$alpha, data[i,]$beta, data[i,]$mean, data[i,]$panel_size)
	if (typeof(pval) == 'list') {
		data[i,8] <- round(as.numeric(mpfr_to_str(phred_mpfr(pval))))
		#data[i, 8] <- 1
		} else {
		data[i,8] <- pval
		}
	#suppressWarnings(data[i,8] <- as.numeric(pval))
	}

unphred_numeric <- function(value) {
        if ((typeof(value) == 'character')&&(value == 'NA')) {
                return('NA')
                } else {
                output <- 10**(-1*value/10)
                return(output)
                }
	}

seeds <- unique(data$seed)
y <- df <- data.frame(seed = character(),
                 pval = integer(),
                 stringsAsFactors = FALSE)
for (i in 1:(length(seeds))) {
	p_vector <- data[data$seed == seeds[i],]$pval
	p_vector <- as.numeric(p_vector[!p_vector %in% "NA"])
	p_vector <- na.omit(p_vector)
	pval <- 'NA'
	if (length(p_vector) > 0) {
		print(p_vector)
		#pval <- fisher(unlist(lapply(p_vector, unphred_numeric)))$p
		#if (pval == 0) {
		#	if (min(p_vector) > 5000) {
		#		pval <- min(p_vector)
		#		} else {
				pval <- round(as.numeric(phred_mpfr(mpfr_to_str(fisher_mpfr(make_mpfr_vector(p_vector))))))
		#		}
		#	}
		}
	y[nrow(y) + 1,] = c(seeds[i],pval)
	}

write.table(y, file = file_out_total, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data[,c("seed", "AD", "DP", "pvalue", "alpha", "beta")], file = file_out_detailed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)








