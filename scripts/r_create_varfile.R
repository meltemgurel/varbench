library(ggplot2)
library(Cairo)

args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','T','G')

mutate <- function(df){
  idx <- seq(10, nrow(df), 101)
  df[idx, 'vaf'] <- round(((rbeta(length(idx), 1, 5, ncp = 0)*9.9)+0.1)/100, 3)
  df[idx, 'alt'] <- sapply(df[idx, 'ref'], function(r) sample(nucs[nucs!=r], 1))
  return(df[idx,])
}

#section1-create and save the varfile
data <- droplevels(subset(read.delim(args[1], header = FALSE,
                   col.names = c('chr', 'pos', 'ref', 'cov')),
                   cov >= mean(cov)))

mutated <- do.call("rbind", by(data, data$chr, mutate))

write.table(x = mutated[,c('chr','pos','pos','vaf','alt')], file = args[2], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')

#section2-create and save the region info for checksum
write.table(x = mutated[,c('chr','pos','ref','alt','vaf')], file = args[3], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')

#section3-create and save vaf plot for reporting
ggplot(mutated, aes(vaf)) + geom_histogram(binwidth = 0.005)
ggsave(file=args[4], type = "cairo-png")
