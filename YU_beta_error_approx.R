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

pval_threshold_for_distr <- 0.00001
baselineDataFraction <- 0.3
baselineDataCount <- 15

panel_size <- as.numeric(args[8])
ADobs <- as.numeric(args[9])
DPobs <- as.numeric(args[10])

data <- read.table(file = args[7], sep = '\t', header = FALSE, col.names = c('AD', 'DP', 'weight'))
data <- transform(data, AF= AD/DP)
data <- na.omit(data)
data <- data[order(data$AF),]

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

isSuitable <- function(distribution) {
	a <- distribution$estimate[[1]]
	b <- distribution$estimate[[2]]
	variance <- (a*b)/((a+b)*(a+b)*(a+b+1))
	if (a < 1 && b < 1) {
		return(0)
		}
	if (variance > 0.05) {
		return(0)
		} else {
		return(1)
		}
	}

baselineDataCount <- min(baselineDataCount, round(baselineDataFraction*nrow(data)))
dataCurrent <- data[1:baselineDataCount,]
distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
print(dataCurrent)

############### MAKING PRIMER DATA SHIFT ##################
base_shift <- 0
skip <- 0
while (1) {
	if (is.na(distrCurrent$estimate[[1]])) {
		if ((baselineDataCount + base_shift) >= nrow(data)) {
			cat("NA\n")
			q()
			}
		base_shift <- base_shift + 1
		cat("SHIFT_s1: base shift: ", base_shift, " / skip: ", skip,"\n")
		dataCurrent <- data[(1+base_shift):(baselineDataCount + base_shift),]
		distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
		} else if (!isSuitable(distrCurrent)) {
		base_shift <- base_shift + 1
		cat("SHIFT_s2: base shift: ", base_shift, " / skip: ", skip,"\t",distrCurrent$estimate[[1]],"-",distrCurrent$estimate[[2]],"\n")
		dataCurrent <- data[(1+base_shift):(baselineDataCount + base_shift),]
		distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
		} else{
		break
		}
	}

for (i in (1:round(baselineDataCount/4))) {
	if (pbeta(dataCurrent[i,]$AF, distrCurrent$estimate[[1]], distrCurrent$estimate[[2]]) < pval_threshold_for_distr) {
		skip <- i
		cat("SHIFT_s3: base shift: ", base_shift, " / skip: ", skip,"\n")
		} else {
		break
		}
	}

while (1) {
	if ((baselineDataCount + base_shift) >= nrow(data)) {
		cat("NA\n")
		q()
		}
	dataCurrent <- data[(1 + base_shift + skip):(baselineDataCount + base_shift),]
	distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
	if (is.na(distrCurrent$estimate[[1]])) {
		base_shift <- base_shift + 1
		cat("SHIFT_s4: base shift: ", base_shift, " / skip: ", skip,"\n")
		next
		}
	if (!isSuitable(distrCurrent)) {
		base_shift <- base_shift + 1
		cat("SHIFT_s5: base shift: ", base_shift, " / skip: ", skip,"\t",distrCurrent$estimate[[1]],"-",distrCurrent$estimate[[2]],"\n")
		next
		}
	if (pbeta(dataCurrent[1,]$AF, distrCurrent$estimate[[1]], distrCurrent$estimate[[2]]) < pval_threshold_for_distr) {
		base_shift <- base_shift + 1
		cat("SHIFT_s6: base shift: ", base_shift, " / skip: ", skip,"\n")
		next
		}
	break
	}

###########################################

cat("SHIFT_final beta distribution base shift: ", base_shift, " / skip: ", skip,"\n")
print(dataCurrent)
cat("Starting beta distribution est: ",distrCurrent$estimate[[1]], "-", distrCurrent$estimate[[2]], "\n")

for (i in ((baselineDataCount + 1 + base_shift):nrow(data))) {
	pval_local <- 1 - ppb(round(data[i,]$AD), distrCurrent$estimate[[1]], distrCurrent$estimate[[2]], c = data[i,]$DP)
	if (pval_local < pval_threshold_for_distr) {
		break
		} else {
		dataCurrent <- data[(1+base_shift+skip):i,]
		#print(dataCurrent)
		distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
		#print(distrCurrent$estimate[[1]])
		#print(distrCurrent$estimate[[2]])
		}
	#print(dataCurrent)
	cat(i,"\t",data[i,]$AD,'-',data[i,]$DP,':',data[i,]$AF,"\t",pval_local,"\t",pval_threshold_for_distr,"\t",distrCurrent$estimate[[1]],"\t",distrCurrent$estimate[[2]],"\n")
	}
cat("Final beta distribution est: ",distrCurrent$estimate[[1]], "-", distrCurrent$estimate[[2]], "\n")
pval_local <- max(0.00000000001, (1 - ppb(round(ADobs), distrCurrent$estimate[[1]], distrCurrent$estimate[[2]], c = DPobs)))
#cat(round(ADobs),"\t",distrCurrent$estimate[[1]],"\t",distrCurrent$estimate[[2]], "\t",DPobs,"\n")
#cat(ppb(round(ADobs), distrCurrent$estimate[[1]], distrCurrent$estimate[[2]], c = DPobs),"\n")
cat(ADobs, "\t", DPobs, "\t", ADobs/DPobs, "\t pval_local:", pval_local,"\n")
if ((ADobs/DPobs) < mean(dataCurrent$AF)) {
	cat("1\n")
	} else {
	pval_local <- min(1, pval_local*panel_size)
	cat(pval_local,"\n")
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
