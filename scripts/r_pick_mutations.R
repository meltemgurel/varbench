args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','G','T')

source("utils.R")

calculate_af_interval <- function(r)
{
  af   <- as.numeric(r[-c(1:4)])
  sd   <- sd(af, na.rm = TRUE)
  med  <- median(af, na.rm = TRUE)
  cov  <- if (med == 0) NA else sd/med
  lbnd <- quantile(af, 0.025, na.rm = TRUE)
  ubnd <- quantile(af, 0.975, na.rm = TRUE)
  return(signif(c(cov, sd, med, lbnd, ubnd),3))
}

plot_af_densities <- function(r)
{
  af   <- as.numeric(r[-c(1:4)])
  med  <- median(af, na.rm = TRUE)
  lbnd <- quantile(af, 0.025, na.rm = TRUE)
  ubnd <- quantile(af, 0.975, na.rm = TRUE)

  plot(density(af, na.rm = TRUE),
       main = 'Allele Frequency', sub = paste(paste(r[1:3], collapse = ':'),
                                              r[4], sep = "->"))
  abline(v = med, col = 'red', lty = 2)
  abline(v = lbnd, col = 'blue', lty = 2)
  abline(v = ubnd, col = 'blue', lty = 2)
}

vcfs   <- grep('.vcf', args, value=TRUE)
args   <- args[-match(vcfs, args)]
caller <- args[1]
nmuts  <- as.numeric(args[2])

pred.mut.list <- lapply(vcfs, attach_af_to_dataframe, caller)
pred.mut.list <- Map(function(x, i) setNames(x, c('CHR', 'POS', 'REF', 'ALT',
                                             paste('VAF', i, sep = '.'))),
                     pred.mut.list, seq_along(pred.mut.list))

afmatrix <- Reduce(function(...) merge(..., by = c('CHR', 'POS', 'REF', 'ALT'),
                                     all = TRUE), pred.mut.list)
afmatrix[is.na(afmatrix)] <- 0

afsumm <- apply(afmatrix, 1, calculate_af_interval)
afsumm <- cbind(afmatrix[,1:4],t(afsumm))

colnames(afsumm) <- c('CHR', 'POS', 'REF', 'ALT',
                      'coef.variation', 'standart.dev', 'median.allele.freq',
                      'lower.bound', 'upper.bound')
                      
afsumm <- afsumm[order(afsumm[,5]),]

#filtering out mutations with median=0 (median=0 => coef.variation is NA)
afsumm <- afsumm[!is.na(afsumm$coef.variation),]
afsumm <- afsumm[round(seq(1, nrow(afsumm), by = nrow(afsumm)/nmuts)),]

pdf(args[3])
par(mfrow = c(4,2))
apply(afmatrix[rownames(afsumm),], 1, plot_af_densities)
dev.off()

write.table(afsumm, file = args[4], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE, dec = '.')
