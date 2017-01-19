args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','T','G')

vcfprocVardict <- function(clist){
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL')
  clist <- clist[sapply(strsplit(clist$INFO, ";"),
                        function(x) return(if(x[1] == "STATUS=LikelySomatic" & x[3] == "TYPE=SNV") TRUE else FALSE)),
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]
  clist$VAF <- round(sapply(strsplit(clist$TUMOR, ":"),
                      function(x) as.integer(x[3])/as.integer(x[2])), 6)
  return(clist)
}

vcfprocSomaticSniper <- function(clist){
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')
  clist <- clist[sapply(strsplit(clist$TUMOR, ":"), "[")[12,] == 2 &
                   sapply(strsplit(clist$NORMAL, ":"), "[")[12,] == 0,
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]
  af <- t(sapply(sapply(sapply(strsplit(clist$TUMOR, ":"), "[")[5,], strsplit, ","),
                 function(x) round(as.integer(x)/sum(as.integer(x)),6)))
  colnames(af) <- nucs
  clist$VAF <- af[cbind(seq_along(clist$ALT), match(clist$ALT,colnames(af)))]
  return(clist)
}

vcfprocMutect <- function(clist){
  #TODO
  return(clist)
}

calc_qm <- function(vcf, caller, rlist, rsoms){
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)

  if(caller == 'vardict') clist <- vcfprocVardict(clist)
  else if(caller == 'somaticsniper') clist <- vcfprocSomaticSniper(clist)
  else if(caller == 'mutect') clist <- vcfprocMutect(clist)

  csoms <- apply(clist[,1:4], 1, paste, collapse = "|")

  tp <- length(intersect(rsoms, csoms))
  fp <- length(setdiff(csoms, rsoms))
  fn <- length(setdiff(rsoms, csoms))

  a <- tp / (tp + fp)
  b <- tp / (tp + fn)
  c <-  if (a == 0 & b == 0) NA else 2 * a * b / (a + b)
  u <- union(rsoms, csoms)
  af.a <- rlist$RVAF[match(u, rsoms)]
  af.a[is.na(af.a)] <- 0
  af.p <- clist$VAF[match(u, csoms)]
  af.p[is.na(af.p)] <- 0
  d <- summary(m <- lm(af.p ~ af.a))$r.squared
  return(c(a,b,c,d))
}

caller <- args[2]
rlist  <- droplevels(subset(read.table(args[3], header = TRUE),
                           RVAF > 0))[, c('CHR', 'POS', 'REF', 'ALT', 'RVAF')]
rsoms  <- apply(rlist[,-5], 1, paste, collapse = "|")
vcfs   <- list.files(args[4], pattern = paste0("mutations.", caller, ".*.vcf"),
                     full.names = TRUE)
qmatrx <- t(sapply(vcfs, calc_qm, caller, rlist, rsoms))

summ.mean <- apply(qmatrx, 2, mean)
summ.std.error <- apply(qmatrx, 2, function(x) sd(x)/sqrt(length(x)))
colnames(qmatrx) <- c('precision', 'recall', 'F-score', 'Rsq')

pdf(args[6])
med <- median(qmatrx[,'Rsq'], na.rm = TRUE)
lbnd <- quantile(qmatrx[,'Rsq'], 0.025, na.rm = TRUE)
ubnd <- quantile(qmatrx[,'Rsq'], 0.975, na.rm = TRUE)
plot(density(qmatrx[,'Rsq'], na.rm = TRUE), main = 'R squared', sub = caller)
abline(v = med, col = 'red', lty = 2)
abline(v = lbnd, col = 'blue', lty = 2)
abline(v = ubnd, col = 'blue', lty = 2)
dev.off()

qmatrx <- rbind(summ.mean, summ.std.error)
cat("\n----------------------------------------------------------------------\n\n",
    file = args[5], append = TRUE)
cat(paste("Called with ", caller, "\n\n"), file=args[5], append = TRUE)
write.table(qmatrx, file = args[5], quote = FALSE, sep = '\t',
            row.names = TRUE, col.names = TRUE, dec = '.', append = TRUE)
cat("\n--------------------------------------------------------------------\n\n",
    file = args[5], append = TRUE)
