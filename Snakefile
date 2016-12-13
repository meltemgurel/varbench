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
from os.path import join, basename, dirname
from setuptools import setup, find_packages
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

include: "src/dependencies.py"

#------------------------------------------------------------------------------
#--------------------------------------------------------------------- Globals-
#------------------------------------------------------------------------------

# Config
configfile: 'config.yml'

# Path to the reference genome.
REFERENCE = config['REFERENCE']

#FTP connection
FTP = FTPRemoteProvider()

# Link to default reference genome
REFLINK = config['REFLINK']

# Path to the sample reads.
SAMPLE = config['SAMPLE'] if config['SAMPLE'] else sys.exit('ERROR: You must provide a sample file')

# Directory where intermediate files will be written.
OUT_DIR = config['OUT_DIR'] if config['OUT_DIR'] else '.'

# For now threads, will scale to cores
THREADS = config['CORES'] if config['CORES'] else 4

#------------------------------------------------------------------------------
#------------------------------------------------------------------------ Init-
#------------------------------------------------------------------------------

# Check for dependencies
if not check_python(): sys.exit('Dependency problem: python >= 2.7.2 is required')
if not check_bwa(): sys.exit('Dependency problem: bwa >= 0.7.12 not found')
if not check_samtools(): sys.exit('Dependency problem: samtools >= 1.2 not found')
if not check_bcftools(): sys.exit('Dependency problem: bcftools >= 1.2 not found')
if not check_wgsim(): sys.exit('Dependency problem: wgsim not found')
if not check_velvet(): sys.exit('Dependency problem: velvet >= 1.2 not found')
if not check_exonerate():
    url="https://github.com/adamewing/exonerate.git"
    os.system("git clone "+url+"; cd exonerate; git checkout v2.4.0; autoreconf -i;"
    " ./configure && make && make check && make install")
if not check_bamsurgeon(): sys.exit('Dependency problem: bamsurgeon not found')

#------------------------------------------------------------------------------
#----------------------------------------------------------------------- Rules-
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#
# Step 1. Set-up: set up environment
#
# TODO: add the final rule
#
#------------------------------------------------------------------------------

#rule all:
#    input:
#        "report.html"

#------------------------------------------------------------------------------
#
# Step 2. Alignment: Align the SAMPLE to the REFERENCE using BWA
#
#------------------------------------------------------------------------------

# Yes I could skip this intermediate step but I need the index file downstream
rule bwa_index:
    """bwa index a reference"""
    input:
        ref = REFERENCE if REFERENCE else FTP.remote(REFLINK)
    output:
        'reference.bwt'
    message: "Creating an index file for {input}."
    shell: "bwa index {input.ref}"

# Align the sample to reference genome
rule bwa_mem:
    """Run bwa mem"""
    params:
        index = 'reference.bwt'
    input:
        read1 = SAMPLE
    output:
        bam = "sample.aligned.bam"
    log:
        log = "bwa.alignment.log"
    threads:
        THREADS
    message:
        "Running alignment with {threads} threads on {input.read1}."
    shell:
        "bwa mem -t {threads} {params.index} " + \
        "{input.read1} | " + \
        " samtools view -Sb - > {output.bam}"



#------------------------------------------------------------------------------
# notestoself;
# benchmark:
# "benchmarks/somecommand/{sample}.txt"
