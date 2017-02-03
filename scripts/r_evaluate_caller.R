args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','G','T')

source("utils.R")

plot_af_corr <- function(actual, predicted, r2, vcf, caller)
{
  pos.x <- max(actual)
  pos.y <- min(predicted)
  eqn <- bquote(r^2 == .(r2))
  plot(actual, predicted, main=paste(caller, "-", vcf),
       xlab="Actual allele frequencies", ylab="Called allele frequencies",
       pch=20, cex=1)
  abline(lm(predicted ~ actual), col="red")
  text(pos.x, pos.y, eqn, adj=c(1,0))
}

plot_af_distributions <- function(af, data)
{
  col <- data[,as.character(af)]
  points(col~rep(af,length(col)), pch='.')
  points(rep(af,3), c(quantile(col, 0.75),
                      quantile(col, 0.50),
                      quantile(col, 0.25)),
         pch='_')
}

plot_quality_metrics <- function(metric, title, caller)
{
  med <- median(metric, na.rm = TRUE)
  lbnd <- quantile(metric, 0.025, na.rm = TRUE)
  ubnd <- quantile(metric, 0.975, na.rm = TRUE)
  plot(density(metric, na.rm = TRUE), main = title, sub = caller)
  abline(v = med, col = 'red', lty = 2)
  abline(v = lbnd, col = 'blue', lty = 2)
  abline(v = ubnd, col = 'blue', lty = 2)
}

calculate_quality_metrics <- function(vcf, caller, real.mut.df, real.mut.id)
{
  pred.mut.df <- attach_af_to_dataframe(vcf, caller)
  pred.mut.id <- gsub("\\s", "", apply(pred.mut.df[,1:4], 1, paste, collapse = "|"))

  tp <- length(intersect(real.mut.id, pred.mut.id))
  fp <- length(setdiff(pred.mut.id, real.mut.id))
  fn <- length(setdiff(real.mut.id, pred.mut.id))

  a <- tp / (tp + fp)
  b <- tp / (tp + fn)
  c <-  if (a == 0 & b == 0) NA else 2 * a * b / (a + b)

  af.a <- real.mut.df$AVAF
  af.p <- pred.mut.df$VAF[match(real.mut.id, pred.mut.id)]
  af.p[is.na(af.p)] <- 0
  r2 <- signif(summary(m <- lm(af.p ~ af.a))$r.squared,3)

  plot_af_corr(af.a, af.p, r2, vcf, caller)

  return(signif(c(a,b,c,r2), 3))
}

get_predicted_af <- function(vcf, caller, real.mut.df, real.mut.id)
{
  pred.mut.df <- attach_af_to_dataframe(vcf, caller)
  pred.mut.id <- gsub("\\s", "", apply(pred.mut.df[,1:4], 1, paste, collapse = "|"))

  af.p <- pred.mut.df$VAF[match(real.mut.id, pred.mut.id)]
  af.p[is.na(af.p)] <- 0

  return(af.p)
}

caller <- args[2]
real.mut.df <- droplevels(subset(read.table(args[3], header = TRUE),
                          AVAF > 0))[, c('CHR', 'POS', 'REF', 'ALT', 'AVAF')]
real.mut.id <- apply(real.mut.df[,-5], 1, paste, collapse = "|")

vcfs <- list.files(args[4], pattern = paste0("mutations.", caller, ".*.vcf"),
                     full.names = TRUE)
afmatrix <- t(sapply(vcfs, get_predicted_af, caller, real.mut.df, real.mut.id))
colnames(afmatrix) <- real.mut.df$AVAF

pdf(args[5])
plot(apply(afmatrix, 2, median)~real.mut.df$AVAF, pch=20, col='darkred',
     main="Predicted allele frequency distributions", sub=caller,
  xlab="Actual allele frequencies", ylab="Predicted allele frequencies")
sapply(real.mut.df$AVAF, plot_af_distributions, afmatrix)
dev.off()

pdf(args[6])
par(mfrow = c(3,2))
qmmatrix <- t(sapply(vcfs, calculate_quality_metrics,
                     caller, real.mut.df, real.mut.id))
dev.off()

summ.mean <- signif(apply(qmmatrix, 2, mean),3)
summ.sd   <- signif(apply(qmmatrix, 2, sd),3)
colnames(qmmatrix) <- c('precision', 'recall', 'F-score', 'Rsq')

pdf(args[7])
plot_quality_metrics(qmmatrix[,'precision'], 'Precision', caller)
dev.off()

pdf(args[8])
plot_quality_metrics(qmmatrix[,'recall'], 'Recall', caller)
dev.off()

pdf(args[9])
plot_quality_metrics(qmmatrix[,'F-score'], 'F-score', caller)
dev.off()

pdf(args[10])
plot_quality_metrics(qmmatrix[,'Rsq'], 'R squared', caller)
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
