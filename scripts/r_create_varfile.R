args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','T','G')

mutate <- function(df){
  #picking interior bases to prevent low coverage within 6bp window
  idx <- seq(10, nrow(df), 101)
  df[idx, 'vaf'] <- ((rbeta(length(idx), 1, 5, ncp = 0)*9.9)+0.1)/100
  df[idx, 'alt'] <- sapply(df[idx, 'ref'], function(r) sample(nucs[nucs!=r], 1))
  return(df[idx,])
}

#create and save the varfile
data <- droplevels(subset(read.delim(args[1], header = FALSE,
                          col.names = c('chr', 'pos', 'ref', 'cov')),
                    subset=(cov >= mean(cov) & ! ref %in% c('a','c','g','t'))))

mutated <- do.call("rbind", by(data, data$chr, mutate))

write.table(x = mutated[,c('chr','pos','pos','vaf','alt')],
            file = args[2], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')

#create and save the region info for checksum
write.table(x = mutated[,c('chr','pos','ref','alt','vaf')],
            file = args[3], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE, dec = '.')

#create and save variant allele frequency plot for reporting
pdf(args[4])
hist(mutated$vaf, breaks=seq(0.001,0.1,by=0.005), xlab='allele frequency',
     main='allele frequency distribution')
dev.off()
