args <- commandArgs()
countAlt <- as.numeric(args[6]);
countAll <- as.numeric(args[7]);
alpha <- as.numeric(args[8]);

get_weigth <- function(countAlt, countAll, alpha) {
	if (countAlt < 1) {
		countAlt = 1;
		}
	f_step <- (countAlt/countAll)/100;
	step <- 0;
	lower <- -1;
	upper <- -1;
	while (1) {
		f <- step * f_step;
		if ((lower < 0)&(ppois(countAlt, f*countAll) < (1-alpha))) {
			lower <- f;
			}
		if ((upper < 0)&(ppois(countAlt, f*countAll) < (alpha))) {
			upper <- f;
			break;
			}
		step <- step + 1;
		}
	return((log(upper - lower))**2)
	}

if (countAll > 1) {
	cat(get_weigth(countAlt, countAll, alpha), "\n");
	} else {
	cat(1,"\n")
	}
