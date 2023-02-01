list.of.packages <- c(
  "Rmpfr",
  "statmod",
  "poolr"
  )

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
      )
    )
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

unphred_numeric <- function(value) {
        if ((typeof(value) == 'character')&&(value == 'NA')) {
                return('NA')
                } else {
                output <- 10**(-1*as.numeric(value)/10)
                return(output)
                }
        }


args <- commandArgs()
p_vector <- args[6:length(args)]

pval <- fisher(unlist(lapply(p_vector, unphred_numeric)))$p
if (pval == 0) {
       if (min(p_vector) > 5000) {
               pval <- min(p_vector)
               } else {
		p_vector <- make_mpfr_vector(p_vector)
                pval <- round(as.numeric(phred_mpfr(mpfr_to_str(fisher_mpfr(p_vector)))))
               }
       }

#input.v <- make_mpfr_vector(input.v)
#f <- round(as.numeric(phred_mpfr(mpfr_to_str(fisher_mpfr(input.v)))))
#f <- fisher(as.numeric(args[6:length(args)]))
cat(pval)













