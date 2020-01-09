# unit tests for custom_fitdist.R


datafilename <- './fitdistr/custom_fitdistr_units_data.txt'
data <- as.numeric(scan(datafilename, what="", sep="\n", quiet=TRUE))
data_exp <- as.numeric(args[7])
data <- sort(data[0:length(data)])

data <- data[data > 0]

print(paste("Data loaded from ", datafilename))
print(paste("Points in dataset: ", length(data)))

source("./fitdistr/custom_fitdistr.R")

print("If shape equal for: AD[0:76] vs. ADW[0:76]xrep(1, 76)")
res1 = fitdist(data[0:76], 'beta', 'mge', gof = "AD")
#print(res1$estimate)
res2 = fitdist(data, 'beta', 'mge', gof = "ADW", weights = rep(1, 76))
#print(res2$estimate)
print(res1$estimate == res2$estimate)
      
print("If shape equal for: AD[0:75] vs. ADW[0:176]xrep(1, 76) + rep(0, 101)")
res1 = fitdist(data[0:75], 'beta', 'mge', gof = "AD")
print(res1$estimate)
res2 = fitdist(c(data, rep(0.013, 100)), 'beta', 'mge', gof = "ADW", weights = c(rep(1, 75), rep(0, 101)))
print(res2$estimate)
print(res1$estimate == res2$estimate)
