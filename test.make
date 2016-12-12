SHELL=/bin/bash -o pipefail

#------------------------------------------------------
#
# Step i. Set-up: set up environment
#
#------------------------------------------------------

# delete failed files
.DELETE_ON_ERROR:

# but keep intermediate files
.SECONDARY:

# export dependencies to the PATH
#export PATH := ./BAMSurgeon/bin/:./samtools/:./bwa/:$(PATH)
export PATH := ./bcftools/:./samtools/:./bwa/:./samtools/misc/wgsim/:./velvetg/:./velveth/:$(PATH)

# Path to script
SCRIPT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# Download links for programs that are not on github
PT_LINK=https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
VV_LINK=https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz

# Execution parameters
CORES=8
THREADS=4
NC_PROCESS=$(CORES)
NP_PROCESS=$(shell expr $(CORES) / $(THREADS))

# Resources and references
SAMPLE=temp.fastq
REFERENCE=GRCh38.primary_assembly.genome.fa

#------------------------------------------------------
#
# Step ii. Installation: install dependencies
#
#------------------------------------------------------

# Install Python dependencies
pythonlibs.version:
	pip install pysam > $@
	pip install cython >> $@

# Install samtools
samtools.version:
	git clone --recursive https://github.com/samtools/htslib.git
	cd htslib; make
	git clone --recursive https://github.com/samtools/samtools.git
	cd samtools; make
	-cd samtools; git log | head -1 > ../$@

#Install BCFTools
bcftools.version:
	git clone https://github.com/samtools/bcftools.git
	cd bcftools; make
	-cd bcftools; git log | head -1 > ../$@

# Install bwa
bwa.version:
	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	-cd bwa; git log | head -1 > ../$@

# Install picard tools
poa.version:
	wget $(PT_LINK)
	unzip picard-tools-1.131.zip
	echo $(PT_LINK) > $@

# Install exonerate
# installing from adamewing fork, for why see:
# https://github.com/adamewing/exonerate/blob/master/README.md
exonerate.version:
	git clone https://github.com/adamewing/exonerate.git
	cd exonerate
	git checkout v2.4.0
	autoreconf -i
	./configure && make && make check && make install
	-cd exonerate; git log | head -1 > ../$@

# Install velvet
velvet.version:
	wget $(VV_LINK)
	tar xvzf velvet_1.2.10.tgz
	make -C velvet_1.2.10
	echo $(VV_LINK) > $@


#------------------------------------------------------
#
# Step 1. Alignment: Align the SAMPLE to the REFERENCE
#										 using BWA
#
#------------------------------------------------------

# Download the primary assembly if no reference genome was specified
$(REFERENCE):
curl -L ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz | gunzip - > $@

# Index the reference for BWA
reference.bwt: $(REFERENCE)
	bwa index $<

# Index the reference for faidx
reference.fasta.fai: $(REFERENCE)
	samtools faidx $<

# Align the reads to the reference
reads.sorted.bam: $(REFERENCE) reference.bwt $(SAMPLE) bwa.version samtools.version
	bwa mem -t $(THREADS) $(REFERENCE) $(SAMPLE) | samtools view -Sb - | samtools sort -f - $@

# Index the sorted reads bam file
reads.sorted.bam.bai: reads.sorted.bam
samtools index $<





#$(error VAR is $(SCRIPT_DIR))
