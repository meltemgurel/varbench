args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','G','T')

source("scripts/utils.R")

calculate_af_interval <- function(r)
{
  af <- as.numeric(r[-c(1:4)])
  med <- median(af, na.rm = TRUE)
  lbnd <- quantile(af, 0.025, na.rm = TRUE)
  ubnd <- quantile(af, 0.975, na.rm = TRUE)

  plot(density(af, na.rm = TRUE),
       main = 'Allele Frequency', sub = paste(paste(r[1:3], collapse = ':'),
                                          r[4], sep = "->"))

  abline(v = med, col = 'red', lty = 2)
  abline(v = lbnd, col = 'blue', lty = 2)
  abline(v = ubnd, col = 'blue', lty = 2)
  return(signif(c(med, lbnd, ubnd),3))
}

caller <- args[2]
vcfs   <- list.files(args[3], pattern = "mutations.*.vcf", full.names = TRUE)

pred.mut.list <- lapply(vcfs, attach_af_to_dataframe, caller)
pred.mut.list <- Map(function(x, i) setNames(x, c('CHR', 'POS', 'REF', 'ALT', paste('VAF', i, sep = '.'))),
              pred.mut.list, seq_along(pred.mut.list))
afmatrix <- Reduce(function(...) merge(..., by = c('CHR', 'POS', 'REF', 'ALT'),
                                  all = TRUE), pred.mut.list)
afmatrix[is.na(afmatrix)] <- 0

pdf(args[4])
par(mfrow = c(4,2))
afintervals <- apply(afmatrix, 1, calculate_af_interval)
dev.off()

afmatrix <- cbind(afmatrix[,1:4],t(afintervals))
colnames(afmatrix) <- c('CHR', 'POS', 'REF', 'ALT',
                      'median.allele.freq', 'lower.bound', 'upper.bound')

write.table(afmatrix, file = args[5], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE, dec = '.')
