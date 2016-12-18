library(data.table)
library(reshape2)
library(ggplot2)
library(Cairo)

poss <- function(v){
  m <- c()
  p <- v[1]
  while(is.finite(p)){
    m <- append(m, p)
    p <- min(v[v >= p + 101])
  }
  return(m)
}

args = commandArgs(trailingOnly=TRUE)
args
#section1-create and save the varfile
data <- droplevels(subset(read.delim(args[1], header = FALSE,
                   col.names = c('chr', 'pos', 'nuc', 'cov')),
                   cov >= mean(cov)))
data$mut <- 0

muts <- melt(sapply(tapply(data$pos, data$chr, poss), function(m) m))[,2:3]
colnames(muts) <- keys <- c('chr', 'pos')
muts$pos2 <- muts$pos
muts$vaf <- round(((rbeta(nrow(muts), 1, 5, ncp = 0)*9.9)+0.1)/100, 3)

write.table(x = muts, file = args[2], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')

#section2-create and save plots for reporting
tData <- data.table(data, key=keys)
tMutd <- data.table(muts, key=keys)
tData[tMutd, mut := 1L]

p <- ggplot(tData, aes(pos, as.factor(mut)))
p + geom_point(aes(colour = as.factor(mut))) + facet_wrap(~chr)
ggsave(file=args[3], type = "cairo-png")

ggplot(muts, aes(vaf)) + geom_histogram(binwidth = 0.005)
ggsave(file=args[4], type = "cairo-png")

#section3-create and save the region info for checksum
tData[tData$mut > 0, mut := muts$vaf]
write.table(x = tData[mut > 0], file = args[5], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')
