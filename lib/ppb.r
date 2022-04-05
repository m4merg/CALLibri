options(warn=1)
args <- commandArgs()
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "poolr",
  "scModels"
  )

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
      )
    )
}

file_in <- args[6]
file_out_total <- args[7]
file_out_detailed <- args[8]
data <- read.table(file = file_in, sep = '\t', header = FALSE, col.names = c('seed', 'AD', 'DP', 'alpha', 'beta', 'mean', 'panel_size'), stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "nA"))
#data <- cbind(data, pvalue)
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
		pval_local <- max(0.00000000001, (1 - ppb(floor(round_custom(AD) - 1), alpha, beta, c = DP)))
		if ((AD/DP) < mean) {
			return(1)
			} else {
			pval_local <- min(1, pval_local*panel_size)
			return(pval_local)
			}
		}
	return(pval_local)
	}

#my.cluster <- parallel::makeCluster(13, type = "PSOCK")
#print(my.cluster)
#doParallel::registerDoParallel(cl = my.cluster)
#foreach::getDoParRegistered()
#foreach::getDoParWorkers()


for (i in 1:nrow(data)) {
	pval <- get_pval(data[i,]$AD, data[i,]$DP, data[i,]$alpha, data[i,]$beta, data[i,]$mean, data[i,]$panel_size)
	data[i,8] <- pval
	}

seeds <- unique(data$seed)
y <- df <- data.frame(seed = character(),
                 pval = integer(),
                 stringsAsFactors = FALSE)
#y <- foreach (i = 1:(length(seeds)), .combine = rbind, .packages='poolr') %do% {
for (i in 1:(length(seeds))) {
	p_vector <- data[data$seed == seeds[i],]$pval
	p_vector <- as.numeric(p_vector[!p_vector %in% "NA"])
	p_vector <- na.omit(p_vector)
	pval <- 'NA'
	if (length(p_vector) > 0) {
		pval <- fisher(p_vector)$p
		}
	#data.frame(seed = seeds[i], pval = pval)
	y[nrow(y) + 1,] = c(seeds[i],pval)
	}

write.table(y, file = file_out_total, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data[,c("seed", "AD", "DP", "pvalue", "alpha", "beta")], file = file_out_detailed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#warnings()
#print(y[y$seed == 'P8ahHtkQeLvfp32QoAbFxDGbHcqpgEep4ua',])
#print(x[x$seed == 'P8ahHtkQeLvfp32QoAbFxDGbHcqpgEep4ua',])








