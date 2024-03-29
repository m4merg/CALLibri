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
#ADobs <- as.numeric(args[9])
#DPobs <- as.numeric(args[10])
ADobs <- 1;
DPobs <- 1;

data <- read.table(file = args[7], sep = '\t', header = FALSE, col.names = c('AD', 'DP', 'weight'))
data <- transform(data, AF= AD/DP)
data <- na.omit(data)
data <- data[order(data$AF),]

gap_measure <- function(data) {
	measure <- (max(data$AF)*length(data$AF)-sum(data$AF))/length(data$AF)
	return(measure)
	}

round_custom <- function(value) {
	if (value < 1) {
		return(0)
		}
	return(round(value))
	}

fitCustom <- function(data, weights) {
	gof <- "ADW"
	customFit <- fitdistC(data = data, "beta", method="mge", gof=gof, weights = weights)
	return(customFit)
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

if (nrow(data) < baselineDataCount) {
	cat("NA\n")
	q()
	}

baselineDataCount <- min(baselineDataCount, round(baselineDataFraction*nrow(data)))
dataCurrent <- data[1:baselineDataCount,]
distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)

cat("WHOLE DATASET:\n")
print(data)
cat("\nDATA PRIMER:\n")
print(dataCurrent)
cat("SHIFTING DATA PRIMER FOR SUITABLE DISTRIBUTION\n")

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
cat("WORKING DATA PRIMER\n")
print(dataCurrent)
cat("Starting beta distribution est: ",distrCurrent$estimate[[1]], "-", distrCurrent$estimate[[2]], " / gap: ",gap_measure(dataCurrent),"\n")
cat("ITERATIONAL SIGNAL SEARCH\n")

A <- distrCurrent$estimate[[1]]
B <- distrCurrent$estimate[[2]]
for (i in ((baselineDataCount + 1 + base_shift):nrow(data))) {
	ad_actual <- floor(round_custom(data[i,]$AD) - 1)
	#cat("ad_actual: ",ad_actual," current 1: ",A," current 2: ",B, " c: ",data[i,]$DP,"\n")
	pval_local <- 1 - ppb(ad_actual, A, B, c = data[i,]$DP)
	left_count <- nrow(data) - i
	#cat(pval_local,"\n")
	if (pval_local < pval_threshold_for_distr) {
		if (((left_count/(nrow(data) + 1)) > 0.1)&&((ppois(ad_actual, data[i,]$DP/10000) > pval_threshold_for_distr) || (ad_actual < 4))&&((pval_local > pval_threshold_for_distr/1000) || (ad_actual < 4))) {
			data[i:nrow(data),]$weight <- 5
			} else {
			cat("FAILED: ",i," (",nrow(data)-i," left)\tAD: ",data[i,]$AD,' (actual: ',ad_actual,') /DP:',data[i,]$DP,' /AF:',data[i,]$AF,"\t /P:",pval_local,"\t /P:",pval_threshold_for_distr,"\t /A",distrCurrent$estimate[[1]],"\t /B",distrCurrent$estimate[[2]],"\n")
			break
			}
		}
	dataCurrent <- data[(1+base_shift+skip):i,]
	distrCurrent <- fitCustom(data = dataCurrent$AF, weights = dataCurrent$weight)
	if ((!(is.na(distrCurrent$estimate[[1]]))) && (!(is.na(distrCurrent$estimate[[2]])))) {
		A <- distrCurrent$estimate[[1]]
		B <- distrCurrent$estimate[[2]]
		}
	#print((distrCurrent$estimate[[1]]))
	#print((distrCurrent$estimate[[2]]))
	#print(dataCurrent)
	cat("PASSED: ",i," (",nrow(data)-i," left)\tAD: ",data[i,]$AD,' (actual: ',ad_actual,') /DP:',data[i,]$DP,' /AF:',data[i,]$AF,"\t /P:",pval_local,"\t /P:",pval_threshold_for_distr,"\t /A",distrCurrent$estimate[[1]],"\t /B",distrCurrent$estimate[[2]],"\n")
	}
cat("Final beta distribution est: ",distrCurrent$estimate[[1]], "-", distrCurrent$estimate[[2]], "\n")
cat(distrCurrent$estimate[[1]],"\t",distrCurrent$estimate[[2]],"\t",mean(dataCurrent$AF));
q()
pval_local <- '1'
if (DPobs > 0) {
	pval_local <- max(0.00000000001, (1 - ppb(floor(round_custom(ADobs) - 1), distrCurrent$estimate[[1]], distrCurrent$estimate[[2]], c = DPobs)))
	}

cat("AD sample: ",ADobs, "\t DP sample:", DPobs, "\t AF sample:", ADobs/DPobs, "\t Pval sample:", pval_local,"\n")
if (DPobs == 0) {
	cat("NA\n")
	} else {
	if ((ADobs/DPobs) < mean(dataCurrent$AF)) {
		cat("1\n")
		} else {
		pval_local <- min(1, pval_local*panel_size)
		cat(pval_local,"\n")
		}
	}

q()


