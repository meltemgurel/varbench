args = commandArgs(trailingOnly=TRUE)
nucs = c('A','C','T','G')

calc_qm <- function(vcf, rlist, rsoms){
  clist<-read.table(vcf, header=FALSE, stringsAsFactors = FALSE)
  colnames(clist) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                       'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')
  clist <- clist[sapply(strsplit(clist$TUMOR, ":"), "[")[12,] == 2 &
                   sapply(strsplit(clist$NORMAL, ":"), "[")[12,] == 0,
                 c('CHR', 'POS', 'REF', 'ALT', 'TUMOR', 'NORMAL')]

  af <- t(sapply(sapply(sapply(strsplit(clist$TUMOR, ":"), "[")[5,], strsplit, ","),
                 function(x) round(as.integer(x)/sum(as.integer(x)),3)))
  colnames(af) <- nucs

  clist$VAF <- af[cbind(seq_along(clist$ALT), match(clist$ALT,colnames(af)))]

  csoms <- apply(clist[,1:4], 1, paste, collapse = "|")

  tp <- length(intersect(rsoms, csoms))
  fp <- length(setdiff(csoms, rsoms))
  fn <- length(setdiff(rsoms, csoms))

  a <- tp / (tp + fp)
  b <- tp / (tp + fn)
  c <-  if (a == 0 & b == 0) 0 else 2 * a * b / (a + b)
  u <- union(rsoms, csoms)
  af.a <- rlist$RVAF[match(u, rsoms)]
  af.a[is.na(af.a)] <- 0
  af.p <- clist$VAF[match(u, csoms)]
  af.p[is.na(af.p)] <- 0
  d = round(summary(m <- lm(af.p ~ af.a))$r.squared, 4)
}

rlist <- droplevels(subset(read.table(args[2], header=TRUE),
                           RVAF > 0))[, c('CHR', 'POS', 'REF', 'ALT', 'RVAF')]
rsoms <- apply(rlist[,-5], 1, paste, collapse = "|")
vcfs <- list.files(args[3], pattern = "mutations.*.vcf")
qms <- sapply(vcfs, calc_qm, rlist, rsoms)
