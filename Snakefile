'''
varbench.snakefile

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
import wget
from os.path import join, basename, dirname

# Globals ---------------------------------------------------------------------

configfile: 'config.yml'

# Path to the reference genome.
REFERENCE = config['REFERENCE']

# Path to the sample reads.
SAMPLE = config['SAMPLE']

# Directory where intermediate files will be written.
OUT_DIR = config['OUT_DIR']

# Functions -------------------------------------------------------------------

def get_reference():
    url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz'
    ref = wget.download(url)
    os.system("gunzip "+ref)
    return ref

# Rules -----------------------------------------------------------------------

#------------------------------------------------------------------------------
#
# Step 1. Set-up: set up environment
#
# TODO: add package checks before installation
#       add actual resource checks: reference + sample file
#
#------------------------------------------------------------------------------

rule install_exonerate:
    params:
        url="https://github.com/adamewing/exonerate.git",
        dir=OUT_DIR
    log:
        "dependencies.exonerate.log"
    shell:
        "git clone {params.url}; cd exonerate; git checkout v2.4.0; autoreconf -i;"
        " ./configure && make && make check && make install"

rule install_bamsurgeon:
    params:
        url="https://github.com/adamewing/bamsurgeon.git",
        dir=OUT_DIR
    log:
        "dependencies.bamsurgeon.log"
    shell:
        "git clone {params.url}; cd bamsurgeon; python setup.py build; "
        "python setup.py install"

rule all:
    input:
        'out.txt'

#------------------------------------------------------------------------------
#
# Step 2. Alignment: Align the SAMPLE to the REFERENCE using BWA
#
#------------------------------------------------------------------------------

rule bwa_index:
    """bwa index a reference"""
    input:
        indexref = REFERENCE if REFERENCE else get_reference()
    output:
        'reference.bwt'
    shell: "bwa index {input.indexref}"
