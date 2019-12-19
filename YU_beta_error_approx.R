pval <- 0.000000001
baselineDataFraction <- 0.6

suppressMessages(library("fitdistrplus"))
args <- commandArgs()
inputFile <-args[6];
data <- as.numeric(scan(inputFile, what="", sep="\n", quiet=TRUE))
data_exp <- as.numeric(args[7])
data <- sort(data[0:length(data)])

fitCustom <- function(data) {
	gof <- "CvM"
	customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof)
	if (is.na(customFit$estimate[[1]])) {
		gof <- "CvM"
		customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof)
		}
	if (is.na(customFit$estimate[[2]])) {
		gof <- "CvM"
		customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof)
		}
	if (customFit$estimate[[1]] < 1) {
		customFit <- fitdistrplus::fitdist(data, "beta", method="mge", gof=gof, fix.arg=list(shape1=1))
		}
	return(customFit)
	}

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













