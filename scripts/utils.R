# Contains some widely used utility functions.

# Processes vcf files generated with vardict:
# Keeps only somatic mutations and calculates allele frequency
process_vcf_Vardict <- function(clist)
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

# Processes vcf files generated with SomaticSniper:
# Keeps only somatic mutations and calculates allele frequency
process_vcf_SomaticSniper <- function(clist)
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

# Processes vcf files generated with Mutect:
# Keeps only somatic mutations and calculates allele frequency
process_vcf_Mutect <- function(clist)
{
  #TODO
  return(clist)
}

# Attaches the variant allele frequency column to the dataframe
attach_af_to_dataframe <- function(vcf, caller)
{
  clist<-read.table(vcf, header = FALSE, stringsAsFactors = FALSE)
  if(caller == 'vardict')
  {
    clist <- process_vcf_Vardict(clist)
  }
  else if(caller == 'somaticsniper')
  {
    clist <- process_vcf_SomaticSniper(clist)
  }
  else if(caller == 'mutect')
  {
    clist <- process_vcf_Mutect(clist)
  }
  return(clist[,c('CHR', 'POS', 'REF', 'ALT', 'VAF')])
}
