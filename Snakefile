'''
varbench pipeline

write smth about the pipeline here.. what it does etc.
-----------------------------------------------------------------------

Requirements:
  See environment.yml

Usage:
  snakemake \
  	--snakefile Snakefile \
  	--configfile config.yml
'''

import json
import sys
import subprocess
import wget
from os.path import join, basename, dirname
from setuptools import setup, find_packages

include: "src/functions.py", "src/functions.R"

#------------------------------------------------------------------------------
#--------------------------------------------------------------------- Globals-
#------------------------------------------------------------------------------

# Config
configfile: 'config.yml'

# Path to the reference genome.
REFERENCE = config['REFERENCE'] if config['REFERENCE'] else sys.exit('ERROR: You must provide a reference file')

# Path to the sample reads.
SAMPLE = config['SAMPLE'] if config['SAMPLE'] else sys.exit('ERROR: You must provide a sample file')

# Directory where intermediate files will be written.
OUT_DIR = config['OUT_DIR'] if config['OUT_DIR'] else 'output/'

# For now threads, will scale to cores
NTHREADS = config['NTHREADS'] if config['NTHREADS'] else 4

#------------------------------------------------------------------------------
#----------------------------------------------------------------------- Rules-
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Step 0. Set-up: set up environment
#
# TODO: add the final rule
#------------------------------------------------------------------------------

def get_name(x):
    return basename(os.path.splitext(x)[0])

rule all:
    input:
        "report.txt"

#------------------------------------------------------------------------------
# Step 1. Alignment: Align the SAMPLE to the REFERENCE and sort
#------------------------------------------------------------------------------

# Align the sample to reference genome
rule bwa_map:
    """Run bwa mem"""
    input:
        reference = REFERENCE,
        reads = SAMPLE
    output:
        temp(join(OUT_DIR, "{prefix}.aligned.bam"))
    params:
        rg="@RG\tID:"+get_name(SAMPLE)+"\tSM:"+get_name(SAMPLE),
        tc=8
    log:
        join(OUT_DIR, "logs/bwa_map/"+get_name(SAMPLE)+".log")
    message:
        "Running alignment with {params.tc} threads on {input.reads}."
    run:
        #create reference index if it doesn't exist
        if not os.path.isfile(input.reference+'.bwt'):
            shell("bwa index "+input.reference)
        #run bwa mem to align the reads
        shell("(bwa mem -R '{params.rg}' -t {params.tc} {input} | "
              "samtools view -Sb - > {output}) 2> {log}")

# Sort the aligned reads
rule samtools_sort:
    input:
        join(OUT_DIR, "{prefix}.aligned.bam")
    output:
        protected(join(OUT_DIR, "{prefix}.sorted.bam"))
    shell:
        "samtools sort -T "+OUT_DIR+"{wildcards.prefix} "
        "-O bam {input} > {output}"

# Index the aligned reads
rule samtools_index:
    input:
        join(OUT_DIR, "{prefix}.sorted.bam")
    output:
        join(OUT_DIR, "{prefix}.sorted.bam.bai")
    shell:
        "samtools index {input}"

#------------------------------------------------------------------------------
# Step 2. Coverage: Identify coverage using bedtools
#------------------------------------------------------------------------------

# Identify coverage
## Generates per-base-coverage table (chr pos base coverage)
rule samtools_mpileup:
    input:
        reference = REFERENCE,
        reads = join(OUT_DIR, "{prefix}.sorted.bam")
    output:
        temp(join(OUT_DIR, "{prefix}.regions"))
    params:
        maxdepth=10000, # At a position, read maximally INT reads per input file (default is 8000)
        mindepth=10
    shell:
        "samtools mpileup -f {input.reference} {input.reads} -d {params.maxdepth} | "
        "awk '$4 > {params.mindepth} {printf \"%s\t%s\t%s\t%s\n\",$1,$2,$3,$4}'"

# Select bases to create the varfile
rule r_create_varfile:
    input:
        join(OUT_DIR, "{prefix}.regions")
    output:
        varfile=join(OUT_DIR, "{prefix}.varfile"),
        mutpos=join(OUT_DIR, "{prefix}.mutations.pos.png"),
        mutfreq=join(OUT_DIR, "{prefix}.mutations.freq.png")
    run:
        R("""
        #section1-create and save the varfile
        data <- read.delim("{input}", header = FALSE)
        colnames(data) <- c('chr', 'pos', 'nuc', 'cov')
        data$mut <- 0

        pos <- melt(sapply(tapply(data$pos, data$chr, poss), function(m) m))[,2:3]
        vaf <- (rbeta(nrow(pos), 1, 5, ncp = 0)*9.9)+0.1
        pos$vaf <- round(vaf, 3)

        write.table(x = pos, file = "{output.varfile}", quote = FALSE, sep = '\t',
                    row.names = FALSE, col.names = FALSE, dec = '.')

        #section2-create and save plots for reporting
        colnames(pos) <- c('chr', 'pos')
        keys <- colnames(pos)

        tData <- data.table(data, key=keys)
        tMutd <- data.table(pos, key=keys)

        tData[tMutd, mut := 1L]
        tData$mut <- as.factor(tData$mut)

        p <- ggplot(tData, aes(pos, mut))
        p + geom_point(aes(colour = mut)) + facet_wrap(~chr)
        ggsave(file="{output.mutpos}")

        ggplot(pos, aes(vaf)) + geom_histogram(binwidth = 0.1)
        ggsave(file="{output.mutfreq}")
        """)

#------------------------------------------------------------------------------
# Step 3. Mutations: Generating the desired somatic mutations with bamsurgeon
#------------------------------------------------------------------------------

# Generate mutated bam files
rule bamsurgeon_addsnv:
    input:
        coverregs=join(OUT_DIR, "{prefix}.varfile"),
        targetbam=join(OUT_DIR, "{prefix}.sorted.bam"),
        reference=REFERENCE
    output:
        join(OUT_DIR, "{prefix}.mut.bam")
    params:
        "-m 0.1 -n 10 -p 8 --mindepth --maxdepth --single --tagreads "

    shell:
        "awk '{$NF=""; print $0}' {input.coverregs} >  "+join(OUT_DIR, "{prefix}.varfile") +
        "; addsnv.py {params} -v {input.coverregs} -f {input.targetbam} -r {input.reference} "
        "-o {output}"

# Temp
rule finalize:
    input:
        expand(join(OUT_DIR, "{prefix}.mut.bam"), prefix=get_name(SAMPLE)),
        expand(join(OUT_DIR, "{prefix}.mutations.pos.png"), prefix=get_name(SAMPLE)),
        expand(join(OUT_DIR, "{prefix}.mutations.freq.png"), prefix=get_name(SAMPLE))
    output:
        "report.txt"
    shell:
        "echo DONE > {output}"

#Get targeted regions


#------------------------------------------------------------------------------
# notestoself;
# benchmark:
# "benchmarks/somecommand/{sample}.txt"
