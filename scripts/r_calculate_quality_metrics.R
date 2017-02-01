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

drawCorr <- function(actual, predicted, r2, vcf, caller)
{
  pos.x <- max(actual)
  pos.y <- min(predicted)
  eqn <- bquote(r^2 == .(r2))
  plot(actual, predicted, main=paste(caller, "-", vcf),
       xlab="Actual allele frequencies", ylab="Called allele frequencies", pch=20,
       cex=1)
  abline(lm(predicted ~ actual), col="red")
  text(pos.x, pos.y, eqn, adj=c(1,0))
}

calc_qm <- function(vcf, caller, rlist, rsoms)
{
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)

  if(caller == 'vardict') clist <- vcfprocVardict(clist)
  else if(caller == 'somaticsniper') clist <- vcfprocSomaticSniper(clist)
  else if(caller == 'mutect') clist <- vcfprocMutect(clist)

  csoms <- gsub("\\s", "", apply(clist[,1:4], 1, paste, collapse = "|"))

  tp <- length(intersect(rsoms, csoms))
  fp <- length(setdiff(csoms, rsoms))
  fn <- length(setdiff(rsoms, csoms))

  a <- tp / (tp + fp)
  b <- tp / (tp + fn)
  c <-  if (a == 0 & b == 0) NA else 2 * a * b / (a + b)

  af.a <- rlist$AVAF
  af.p <- clist$VAF[match(rsoms, csoms)]
  af.p[is.na(af.p)] <- 0
  r2 <- signif(summary(m <- lm(af.p ~ af.a))$r.squared,3)

  drawCorr(af.a, af.p, r2, vcf, caller)

  return(signif(c(a,b,c,r2), 3))
}

draw_qm <- function(val, title, caller)
{
  med <- median(val, na.rm = TRUE)
  lbnd <- quantile(val, 0.025, na.rm = TRUE)
  ubnd <- quantile(val, 0.975, na.rm = TRUE)
  plot(density(val, na.rm = TRUE), main = title, sub = caller)
  abline(v = med, col = 'red', lty = 2)
  abline(v = lbnd, col = 'blue', lty = 2)
  abline(v = ubnd, col = 'blue', lty = 2)
}

get_af <- function(vcf, caller, rlist, rsoms)
{
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)

  if(caller == 'vardict') clist <- vcfprocVardict(clist)
  else if(caller == 'somaticsniper') clist <- vcfprocSomaticSniper(clist)
  else if(caller == 'mutect') clist <- vcfprocMutect(clist)

  csoms <- gsub("\\s", "", apply(clist[,1:4], 1, paste, collapse = "|"))

  af.p <- clist$VAF[match(rsoms, csoms)]
  af.p[is.na(af.p)] <- 0

  return(af.p)
}

plotqm <- function(avaf, data)
{
  col <- data[,as.character(avaf)]
  points(col~rep(avaf,length(col)), pch='.')
  points(rep(avaf,3), c(quantile(col, 0.75),
                        quantile(col, 0.50),
                        quantile(col, 0.25)),
         pch='_')
}

caller <- args[2]
rlist  <- droplevels(subset(read.table(args[3], header = TRUE),
                           AVAF > 0))[, c('CHR', 'POS', 'REF', 'ALT', 'AVAF')]
rsoms  <- apply(rlist[,-5], 1, paste, collapse = "|")
vcfs   <- list.files(args[4], pattern = paste0("mutations.", caller, ".*.vcf"),
                     full.names = TRUE)

afmatrx <- t(sapply(vcfs, get_af, caller, rlist, rsoms))
colnames(afmatrx) <- rlist$AVAF

pdf(args[5])
plot(apply(afmatrx, 2, median)~rlist$AVAF, pch=20, col='darkred',
     main="Predicted allele frequency distributions", sub=caller,
  xlab="Actual allele frequencies", ylab="Predicted allele frequencies")
sapply(rlist$AVAF, plotqm, afmatrx)
dev.off()

pdf(args[6])
par(mfrow = c(3,2))
qmatrx <- t(sapply(vcfs, calc_qm, caller, rlist, rsoms))
dev.off()

summ.mean <- signif(apply(qmatrx, 2, mean),3)
summ.sd   <- signif(apply(qmatrx, 2, sd),3)
colnames(qmatrx) <- c('precision', 'recall', 'F-score', 'Rsq')

pdf(args[7])
draw_qm(qmatrx[,'precision'], 'Precision', caller)
dev.off()

pdf(args[8])
draw_qm(qmatrx[,'recall'], 'Recall', caller)
dev.off()

pdf(args[9])
draw_qm(qmatrx[,'F-score'], 'F-score', caller)
dev.off()

pdf(args[10])
draw_qm(qmatrx[,'Rsq'], 'R squared', caller)
dev.off()

qmsumm <- rbind(summ.mean, summ.sd)
colnames(qmsumm) <- c('precision', 'recall', 'F-score', 'Rsq')
cat("\n----------------------------------------------------------------------\n\n",
    file = args[11], append = TRUE)
cat(paste("Called with ", caller, "\n\n"), file=args[5], append = TRUE)
write.table(qmsumm, file = args[11], quote = FALSE, sep = '\t',
            row.names = TRUE, col.names = TRUE, dec = '.', append = TRUE)
cat("\n--------------------------------------------------------------------\n\n",
    file = args[11], append = TRUE)
