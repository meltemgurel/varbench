# Varbench

**Varbench** is a suite of tools developed to compare somatic variant callers for targeted deep sequencing data, run the best performing one, and report the allele frequencies of the called variants using estimate intervals.

**Varbench** comprises four pipelines, namely **simulate**, **compare**, **estimate** and **validate**. `doc/report.html` describes each pipeline and reports the results obtained by running them on targeted deep sequencing ctDNA data.

## Prerequisites
To run **Varbench** you need the following software to be installed:
- Python
- Snakemake
- Bamsurgeon
    - Pysam
    - BWA
    - SAMtools
    - BCFtools
    - Velvet
    - Picard tools
    - Exonerate

## Installing dependencies
**Varbench** uses the workflow system Snakemake. Snakemake follows the GNU Make paradigm. With a Python 3 setup, installation of Snakemake can be performed by issuing
```sh
$ easy_install3 snakemake
```
or
```sh
$ pip3 install snakemake
```
in your terminal.

To induce somatic mutations, **Varbench** uses [Bamsurgeon](https://github.com/adamewing/bamsurgeon); a suite of tools for adding mutations to .bam files. See instructions [here](https://github.com/adamewing/bamsurgeon) for installation.

Instead of manually installing Bamsurgeon's dependencies you can create a conda environment with the ```environment.yml``` file provided in the root directory.
```sh
conda create -n myworkflow --file environment.yml
```
This will install all the necessary dependencies in an isolated conda environment. You can then activate this environment by running
```sh
$ source activate myworkflow
```
## Running Varbench
Clone **Varbench** into your working directory
```sh
$ git clone https://github.com/meltemgurel/varbench.git path/to/workdir
$ cd path/to/workdir
```

### simulate: Simulating somatic mutations
After you edit ```config.simulate.yml``` to specify the
- path to the reference genome file (```REFERENCE```),
- path to the paired-end .fastq read files (```SAMPLE1``` and ```SAMPLE2```),
- minimum and maximum coverage depth to filter reads (```MINDEPTH``` and ```MAXDEPTH```, default is 10 and 10K),
- directory to store the pipeline outputs (```OUT_DIR```, will be created if it doesn't exist already),
- number of processes to split into (```NTHREADS```)

run `simulate` with
```sh
snakemake --snakefile simulate --configfile config.simulate.yml
```
### compare: Pick the best variant caller
To run `compare` you will need to install VarDict, SomaticSniper, Strelka, and Varscan. These are already included in the conda environment. If, instead, you manually installed the dependencies please make sure to install these variant callers before you run `compare`. More variant callers will keep on being integrated.

```config.compare.yml``` needs to be edited before running the pipeline. Here the ```NORMAL``` and ```TUMOR``` samples are the .bam files created with the ```simulate``` pipeline, and the ```MUTATIONS``` parameter refers to the list of mutations generated again with the same pipeline.
You also need to submit the
- path to the reference genome file (```REFERENCE```),
- directory to store the pipeline outputs (```OUT_DIR```, will be created if it doesn't exist already),
- number of processes to split into (```NTHREADS```)

To run `compare` issue
```sh
snakemake --snakefile compare --configfile config.compare.yml
```
in your terminal.
### estimate: Calling somatic mutations with allele frequency confidence intervals
Requires at least one of VarDict or SomaticSniper to be installed. And you must specify either 'vardict' or 'somaticsniper' as the preferred ```CALLER``` in the ```config.estimate.yml``` configuration file. Example setup files are provided in the root directory: ```config.p3_vardict.yml``` and ```config.p3_somaticsniper.yml```. The ```NORMAL``` and ```TUMOR``` samples are the .bam files created with the ```simulate``` pipeline.
Also specify the
- path to the reference genome file (```REFERENCE```),
- directory to store the pipeline outputs (```OUT_DIR```, will be created if it doesn't exist already),
- number of processes to split into (```NTHREADS```)

`estimate` can then be executed with:
```sh
snakemake --snakefile estimate --configfile config.estimate.yml
```
### validate: Pick mutations for validation
When you create .VCF files with either ```compare``` or ```estimate``` you can run
```sh
snakemake --snakefile validate --configfile config.validate.yml
```
to receive a list of mutations you should validate for your analysis. The ```VCF_DIR``` parameter in the ```config.validate.yml``` should point to the .VCF files and ```NMUTATIONS``` refers to the number of mutations you wish to validate.
Other required fields are the
- directory to store the pipeline outputs (```OUT_DIR```, will be created if it doesn't exist already),
- variant caller used to generate the .VCF files (```CALLER```)
- number of processes to split into (```NTHREADS```)
