args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','G','T')

vcfprocVardict <- function(clist)
{
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL')

  clist <- clist[sapply(strsplit(clist$INFO, ";"), function(x)
    return(if((x[1] == "STATUS=LikelySomatic" | x[1] == "STATUS=StrongSomatic") & x[3] == "TYPE=SNV") TRUE else FALSE)),
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]

  clist$VAF <- signif(sapply(strsplit(clist$TUMOR, ":"),
                             function(x) as.integer(x[3])/as.integer(x[2])), 3)
  return(clist)
}

vcfprocSomaticSniper <- function(clist)
{
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')

  clist <- clist[sapply(strsplit(clist$TUMOR, ":"), "[")[12,] == 2,
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]

  af <- t(sapply(sapply(sapply(strsplit(clist$TUMOR, ":"), "[")[5,], strsplit, ","),
                 function(x) signif(as.integer(x)/sum(as.integer(x)),3)))

  colnames(af) <- nucs
  clist$VAF <- af[cbind(seq_along(clist$ALT), match(clist$ALT,colnames(af)))]
  return(clist)
}

vcfprocMutect <- function(clist)
{
  #TODO
  return(clist)
}

calc_af <- function(vcf, caller)
{
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)

  if(caller == 'vardict') clist <- vcfprocVardict(clist)
  else if(caller == 'somaticsniper') clist <- vcfprocSomaticSniper(clist)
  else if(caller == 'mutect') clist <- vcfprocMutect(clist)

  return(clist[,c('CHR', 'POS', 'REF', 'ALT', 'VAF')])
}

calc_bootstrap_summary <- function(r)
{
  af <- as.numeric(r[-c(1:4)])
  sd <- sd(af, na.rm = TRUE)
  med <- median(af, na.rm = TRUE)
  cov <- if (med == 0) NA else sd/med
  lbnd <- quantile(af, 0.025, na.rm = TRUE)
  ubnd <- quantile(af, 0.975, na.rm = TRUE)
  return(signif(c(cov, sd, med, lbnd, ubnd),3))
}

draw_allelefreq_densities <- function(r)
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
}

vcfs <- grep('.vcf', args, value=TRUE)
args <- args[-match(vcfs, args)]
caller <- args[1]
nmuts  <- as.numeric(args[2])
aflist <- lapply(vcfs, calc_af, caller)
aflist <- Map(function(x, i) setNames(x, c('CHR', 'POS', 'REF', 'ALT', paste('VAF', i, sep = '.'))),
              aflist, seq_along(aflist))
afmtrx <- Reduce(function(...) merge(..., by = c('CHR', 'POS', 'REF', 'ALT'),
                                     all = TRUE), aflist)
afmtrx[is.na(afmtrx)] <- 0
afsumm <- apply(afmtrx, 1, calc_bootstrap_summary)
afsumm <- cbind(afmtrx[,1:4],t(afsumm))
colnames(afsumm) <- c('CHR', 'POS', 'REF', 'ALT',
                      'coef.variation', 'standart.dev', 'median.allele.freq', 'lower.bound', 'upper.bound')
afsumm <- afsumm[order(afsumm[,5]),]
#filtering out mutations with median=0 (median=0 => coef.variation is NA)
afsumm <- afsumm[!is.na(afsumm$coef.variation),]
afsumm <- afsumm[round(seq(1, nrow(afsumm), by = nrow(afsumm)/nmuts)),]

pdf(args[3])
par(mfrow = c(4,2))
apply(afmtrx[rownames(afsumm),], 1, draw_allelefreq_densities)
dev.off()

write.table(afsumm, file = args[4], quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE, dec = '.')