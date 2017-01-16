args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','T','G')

vcfprocVardict <- function(clist){
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL')
  clist <- clist[sapply(strsplit(clist$INFO, ";"),
                        function(x) return(if(x[1] == "STATUS=LikelySomatic" & x[3] == "TYPE=SNV") TRUE else FALSE)),
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]
  clist$VAF <- round(sapply(strsplit(clist$TUMOR, ":"),
                      function(x) as.integer(x[3])/as.integer(x[2])), 3)
  return(clist)
}

vcfprocSomaticSniper <- function(clist){
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')
  clist <- clist[sapply(strsplit(clist$TUMOR, ":"), "[")[12,] == 2 &
                   sapply(strsplit(clist$NORMAL, ":"), "[")[12,] == 0,
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]
  af <- t(sapply(sapply(sapply(strsplit(clist$TUMOR, ":"), "[")[5,], strsplit, ","),
                 function(x) round(as.integer(x)/sum(as.integer(x)),3)))
  colnames(af) <- nucs
  clist$VAF <- af[cbind(seq_along(clist$ALT), match(clist$ALT,colnames(af)))]
  return(clist)
}

vcfprocMutect <- function(clist){
  #TODO
  return(clist)
}

calc_af <- function(vcf, caller, rlist, rsoms){
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)

  if(caller == 'vardict') clist <- vcfprocVardict(clist)
  else if(caller == 'somaticsniper') clist <- vcfprocSomaticSniper(clist)
  else if(caller == 'mutect') clist <- vcfprocMutect(clist)

  return(clist[,c('CHR', 'POS', 'REF', 'ALT', 'VAF')])
}

calc_metrics <- function(r) {
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
  return(c(med, lbnd, ubnd))
}

caller <- args[2]
vcfs   <- list.files(args[3], pattern = "mutations.*.vcf", full.names = TRUE)
aflist <- lapply(vcfs, calc_af, caller)
aflist <- Map(function(x, i) setNames(x, c('CHR', 'POS', 'REF', 'ALT', paste('VAF', i, sep = '.'))),
              aflist, seq_along(aflist))
afmtrx <- Reduce(function(...) merge(..., by = c('CHR', 'POS', 'REF', 'ALT'),
                                  all = TRUE), aflist)
afmtrx <- afmtrx[!(rowSums(is.na(afmtrx)) > sum(grepl("VAF", names(afmtrx))) - 6),]

pdf(args[4])
par(mfrow = c(4,2))
afdens <- apply(afmtrx, 1, calc_metrics)
dev.off()

afmtrx <- cbind(afmtrx[,1:4],t(afdens))
colnames(afmtrx) <- c('CHR', 'POS', 'REF', 'ALT',
                      'median.allele.freq', 'lower.bound', 'upper.bound')

write.table(afmtrx, file = args[5], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE, dec = '.')
