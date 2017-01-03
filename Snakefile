'''
varbench pipeline - Simulating somatic mutations

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
from snakemake.utils import R

include: "scripts/functions.py"

#------------------------------------------------------------------------------
#--------------------------------------------------------------------- Globals-
#------------------------------------------------------------------------------

# Config
configfile: 'config.yml'

# Path to the reference genome.
REFERENCE = config['REFERENCE'] if config['REFERENCE'] else sys.exit('ERROR: You must provide a reference file')

# Path to the sample reads.
SAMPLE1 = config['SAMPLE1'] if config['SAMPLE1'] else sys.exit('ERROR: You must provide a sample file')
SAMPLE2 = config['SAMPLE2'] if config['SAMPLE2'] else sys.exit('ERROR: You must provide a sample file')
SAMPLE = SAMPLE1

# Directory where intermediate files will be written.
OUT_DIR = config['OUT_DIR'] if config['OUT_DIR'] else 'output/'

# For now threads, will scale to cores
NTHREADS = config['NTHREADS'] if config['NTHREADS'] else 8

#------------------------------------------------------------------------------
# Step 0. Set-up: set up environment
#
# TODO: update the final rule
#------------------------------------------------------------------------------

# Check for dependencies
setup()

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
        reference=REFERENCE,
        read1=SAMPLE1,
        read2=SAMPLE2
    output:
        temp(join(OUT_DIR, "{prefix}.aligned.bam"))
    params:
        rg="@RG\tID:"+get_name(SAMPLE)+"\tSM:"+get_name(SAMPLE),
        tc=NTHREADS
    log:
        join(OUT_DIR, "logs/bwa_map."+get_name(SAMPLE)+".log")
    message:
        "Running alignment with {params.tc} threads on {input.read1} and {input.read2}."
    run:
        #create reference index if it doesn't exist
        if not os.path.isfile(input.reference+'.bwt'):
            shell("bwa index "+input.reference)
        #run bwa mem to align the reads
        shell("(bwa mem -R '{params.rg}' -t {params.tc} {input.reference} {input.read1} {input.read2} | "
              "samtools view -Sb -F 3844 - > {output}) 2> {log}")

# Sort the aligned reads
rule samtools_org_sort:
    input:
        join(OUT_DIR, "{prefix}.aligned.bam")
    output:
        protected(join(OUT_DIR, "{prefix}.org.bam"))
    log:
        join(OUT_DIR, "logs/samtools_org_sort."+get_name(SAMPLE)+".log")
    shell:
        "(samtools sort -T "+OUT_DIR+"{wildcards.prefix} "
        "-O bam {input} > {output}) 2> {log}"

# Index the aligned reads
rule samtools_org_index:
    input:
        join(OUT_DIR, "{prefix}.org.bam")
    output:
        join(OUT_DIR, "{prefix}.org.bam.bai")
    log:
        join(OUT_DIR, "logs/samtools_org_index."+get_name(SAMPLE)+".log")
    shell:
        "(samtools index {input}) 2> {log}"

#------------------------------------------------------------------------------
# Step 2. Coverage: Identify coverage using bedtools
#------------------------------------------------------------------------------

# Identify coverage
## Generates per-base-coverage table (chr pos base coverage)
rule samtools_mpileup:
    input:
        reference=REFERENCE,
        reads=join(OUT_DIR, "{prefix}.org.bam"),
        readsi=join(OUT_DIR, "{prefix}.org.bam.bai")
    output:
        temp(join(OUT_DIR, "{prefix}.regions"))
    params:
        maxdepth=10000, # At a position, read maximally INT reads per input file (default is 8000)
        mindepth=10
    log:
        join(OUT_DIR, "logs/samtools_mpileup."+get_name(SAMPLE)+".log")
    shell:
        "(samtools mpileup -f {input.reference} {input.reads} -d {params.maxdepth} | "
        "awk '$4 > {params.mindepth} {{printf \"%s\\t%s\\t%s\\t%s\\n\",$1,$2,$3,$4}}' "
        "> {output}) 2> {log}"

#------------------------------------------------------------------------------
# Step 3. Mutations: Generating the desired somatic mutations
#
# TODO: make sure the 6bp window around the nuc is covered
# TODO: maybe specify a 6bp window and let bamsurgeon randomly choose the nuc?
# TODO: check R libs, silence warnings
#------------------------------------------------------------------------------

# Select bases to create the varfile
rule r_create_varfile:
    input:
        join(OUT_DIR, "{prefix}.regions")
    output:
        varfile=join(OUT_DIR, "{prefix}.varfile"),
        alleles=join(OUT_DIR, "{prefix}.alleles"),
        mutfreq=join(OUT_DIR, "{prefix}.mutations.freq.png")
    shell:
        "Rscript --vanilla scripts/r_create_varfile.R {input} {output.varfile} "
        "{output.alleles} {output.mutfreq}"

#------------------------------------------------------------------------------
# Step 4. Mutations: Introducing somatic mutations in the original BAM file
#
# TODO: {CW} depth at site is greater than cutoff, aborting mutation.
#------------------------------------------------------------------------------

# Generate mutated bam files
rule bamsurgeon_addsnv:
    input:
        varfile=join(OUT_DIR, "{prefix}.varfile"),
        targetbam=join(OUT_DIR, "{prefix}.org.bam"),
        targetbami=join(OUT_DIR, "{prefix}.org.bam.bai"),
        reference=REFERENCE
    output:
        temp(join(OUT_DIR, "{prefix}.mut.bam"))
    params:
        "-p 128 --mindepth 10 --maxdepth 10000 --ignoresnps --force --tagreads "
        "--tmpdir "+join(OUT_DIR, "addsnv")
    log:
        join(OUT_DIR, "logs/bamsurgeon_addsnv."+get_name(SAMPLE)+".log")
    shell:
        "(addsnv.py {params} -v {input.varfile} -f {input.targetbam} -r {input.reference} "
        "-o {output}) 2> {log}"

# Sort the mutated reads
rule samtools_mut_sort:
    input:
        join(OUT_DIR, "{prefix}.mut.bam")
    output:
        protected(join(OUT_DIR, "{prefix}.mutant.bam"))
    log:
        join(OUT_DIR, "logs/samtools_mut_sort."+get_name(SAMPLE)+".log")
    shell:
        "(samtools sort -T "+OUT_DIR+"{wildcards.prefix} "
        "-O bam {input} > {output}) 2> {log}"

# Index the mutated reads
rule samtools_mut_index:
    input:
        join(OUT_DIR, "{prefix}.mutant.bam")
    output:
        join(OUT_DIR, "{prefix}.mutant.bam.bai")
    log:
        join(OUT_DIR, "logs/samtools_mut_index."+get_name(SAMPLE)+".log")
    shell:
        "(samtools index {input}) 2> {log}"

#------------------------------------------------------------------------------
# Step 5. Checksum: Calculating the allele frequency of the induced mutations
#------------------------------------------------------------------------------

# Calculate the actual allele frequencies for each position in *.alleles
rule py_allele_freq_cal:
    input:
        alleles=join(OUT_DIR, "{prefix}.alleles"),
        original=join(OUT_DIR, "{prefix}.org.bam"),
        modified=join(OUT_DIR, "{prefix}.mutant.bam"),
        modifiedi=join(OUT_DIR, "{prefix}.mutant.bam.bai"),
    output:
        join(OUT_DIR, "{prefix}.mutations")
    shell:
        "python scripts/py_allele_freq_cal.py {input.alleles} {input.original} {input.modified}"

# Temp
rule finalize:
    input:
        expand(join(OUT_DIR, "{prefix}.mutations"), prefix=get_name(SAMPLE)),
        expand(join(OUT_DIR, "{prefix}.mutations.pos.png"), prefix=get_name(SAMPLE)),
        expand(join(OUT_DIR, "{prefix}.mutations.freq.png"), prefix=get_name(SAMPLE))
    output:
        "report.txt"
    shell:
        "echo DONE > {output}"

#------------------------------------------------------------------------------
# notestoself;
# benchmark:
# "benchmarks/somecommand/{sample}.txt"
