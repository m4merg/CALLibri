# execusion:
# R --slave -f <path to script> --args <path to fitdistr> <path to input file>
# example: 
# $ R --slave -f /home/onco-admin/RnD/UEBAcall/YU_beta_error_approx.R --args /home/onco-admin/RnD/UEBAcall/fitdistr/ /home/onco-admin/RnD/UEBAcall/test/indexdata_N397206687986727744
library("modi")
args <- commandArgs()
setwd(args[6])
source("custom_fitdistr.R")
suppressMessages(library("modi"))
suppressMessages(library("scModels"))
#suppressMessages(library("fitdistrplus"))

pval_threshold_for_distr <- 0.1
baselineDataFraction <- 0.3
baselineDataCount <- 15

panelSize <- args[8]
ADobs <- args[9]
DPobs <- args[10]

data <- read.table(file = args[7], sep = '\t', header = FALSE, col.names = c('AD', 'DP', 'weight'))
data <- transform(data, AF= AD/DP)
data <- na.omit(data)
data <- data[order(data$AF),]
print(data)

fitCustom <- function(data, weights) {
	gof <- "ADW"
	customFit <- fitdistC(data = data, "beta", method="mge", gof=gof, weights = weights)
	#if (is.na(customFit$estimate[[1]])) {
	#	gof <- "CvM"
	#	customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof)
	#	}
	#if (is.na(customFit$estimate[[2]])) {
	#	gof <- "CvM"
	#	customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof)
	#	}
	#if (customFit$estimate[[1]] < 1) {
	#	customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof, fix.arg=list(shape1=1))
	#	}
	return(customFit)
	}

getNextProbe <- function(AD, DP, alpa, beta, threshold) {
	#pval_local <- dpb(AD, alpha, beta, c = DP)
	#if (pval_local < threshold) {
	#	return 0
	#	} else {
	#	return 1
	#	}
	}

baselineDataCount <- min(baselineDataCount, round(baselineDataFraction*nrow(data)))
dataCurrent <- data[1:baselineDataCount,]
print(dataCurrent)
distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
print(distrCurrent$estimate[[1]])
print(distrCurrent$estimate[[2]])

for (i in ((baselineDataCount + 1):nrow(data))) {
	pval_local <- dpb(round(data[i,]$AD), distrCurrent$estimate[[1]], distrCurrent$estimate[[2]], c = data[i,]$DP)
	if (pval_local < pval_threshold_for_distr) {

		} else {
		dataCurrent <- data[1:i,]
		distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
		}
	cat(i,"\t",data[i,]$AD,'-',data[i,]$DP,':',data[i,]$AF,"\t",pval_local,"\t",pval_threshold_for_distr,"\t",distrCurrent$estimate[[1]],"\t",distrCurrent$estimate[[2]],"\n")
	}
q()
if (length(data[data>0]) < 10) {
	cat("NA\n");
	} else {
	
	dataCurrent <- data[data>0][1:round(baselineDataFraction*length(data[data>0]))]
	distrCurrent <- fitCustom(dataCurrent)
	dataBaseline <- dataCurrent
	
	for (i in (round(baselineDataFraction*length(data[data>0])) + 1):length(data[data>0]) ) {
		prob <- 1-pbeta(data[data>0][i], distrCurrent$estimate[[1]], distrCurrent$estimate[[2]])
		if (prob > pval) {
			dataCurrent <- c(dataCurrent, data[data>0][i])
			distrCurrent <- fitCustom(dataCurrent)
			#prob <- 1-pbeta(data[data>0][i], distrCurrent$estimate[[1]], distrCurrent$estimate[[2]])
			}
		cat (i,data[data>0][i],prob,pval,distrCurrent$estimate[[1]],distrCurrent$estimate[[2]],"\n", sep="\t")
		}

	prob <- 1-pbeta(data_exp, distrCurrent$estimate[[1]], distrCurrent$estimate[[2]])
	cat(prob,"\t", distrCurrent$estimate[[1]],"\t", distrCurrent$estimate[[2]],"\t", distrCurrent$aic,"\t", distrCurrent$bic,"\t", length(dataCurrent),"/",length(data[data>0]),"\n")
	}
