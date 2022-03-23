library("poolr")
args <- commandArgs()
f <- fisher(as.numeric(args[6:length(args)]))
cat(f$p)
