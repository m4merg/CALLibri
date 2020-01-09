library(parallel)
suppressMessages(library(hash))
options(warn=-1)
args <- commandArgs()

inputFolder <-args[1];
no_cores  <-args[2];
sampleFile <- paste(inputFolder,sample_N1,sep="/")

for (controlFile in list.files(path = inputFolder, pattern = "control")) {
	controlData <- scan(controlFile, what="", sep="\n", quiet=TRUE)
	for (i in (1:length(controlData))) {
		line <- as.character(unlist(strsplit(controlData[i], "\t")))
		str <- paste()
		}
	}
