args <- commandArgs()
countAlt <- as.integer(args[6]);
countAll <- as.integer(args[7]);
alpha <- as.double(args[8]);
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

get_weight <- function(lower, upper) {
	return((log(upper - lower))**2)
        }

cat(get_weight(lower, upper), "\n");
