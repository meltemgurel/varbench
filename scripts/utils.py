# Contains some widely used utility functions.

#Checks if Exonerate is installed
def check_exonerate():
    p = subprocess.Popen(['exonerate'], stdout=subprocess.PIPE)
    for line in p.stdout:
        sline = line.decode("utf-8")
        if sline.startswith('exonerate from exonerate'):
            major, minor = sline.strip().split()[-1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False

#Checks if Bamsurgeon is installed
def check_bamsurgeon():
    p = subprocess.Popen(['addsnv.py'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('usage: addsnv.py'):
            return True
    return False

# Checks for dependencies before running simulate
# Attempts to install Exonerate if it's not installed
def check_dependencies():
    if not check_exonerate():
        url="https://github.com/adamewing/exonerate.git"
        os.system("git clone "+url+"; cd exonerate; git checkout v2.4.0; "
        "touch aclocal.m4 configure Makefile.am && Makefile.in "
        "&& ./configure && make && make install")
    if not check_bamsurgeon():
        sys.exit('Dependency problem: bamsurgeon not found')

# Returns file names for labelling
def get_name(x):
    return basename(os.path.splitext(x)[0])

# Returns the execution command for given variant caller
def get_caller_cmd(caller, ref, muts, outdir):
    if caller == 'vardict':
        return cmd_vardict(ref, muts, outdir)
    elif caller == 'mutect':
        return cmd_mutect(ref, muts, outdir)
    elif caller == 'somaticsniper':
        return cmd_somaticsniper(ref)
    elif caller == 'strelka':
        return cmd_strelka(ref, outdir)

# Returns the execution command for Vardict
def cmd_vardict(ref, muts, outdir):
    in_df  = pd.DataFrame.from_csv(muts, sep='\t',
                                   index_col=None)
    out_df = in_df.groupby(['CHR'], sort=False).agg({'POS' : [min, max]})
    out_df.to_csv(join(outdir, "regs.vardict.bed"), sep="\t",
                  index=True, header=False)
    return 'vardict-java -G '+ref+' -f 0.001 -N "{}" -b "{}|{}" ' \
           '-c 1 -S 2 -E 3 '+join(outdir, "regs.vardict.bed")+' | scripts/testsomatic.R | scripts/var2vcf_paired.pl ' \
           '-N "{}|{}" -f 0.001 > {}'

# Returns the execution command for Mutect
def cmd_mutect(ref, muts, outdir):
    in_df  = pd.DataFrame.from_csv(muts, sep='\t',
                                   index_col=None)
    out_df = in_df.groupby(['CHR'], sort=False).agg({'POS' : [min, max]})
    out_df.to_csv(join(outdir, "regs.mutect.txt"), sep="\t",
                  index=True, header=False)
    os.system('awk {print $1":"$2"-"$2} '+join(outdir, "regs.mutect.txt")+
              ' > '+join(outdir, "regs.mutect.txt"))
    return 'java -Xmx2g -jar mutect/mutect.jar --analysis_type MuTect --reference_sequence '+ref+ \
           '--cosmic mutect/b37_cosmic_v54_120711.vcf --dbsnp mutect/dbsnp_132_b37.leftAligned.vcf.gz ' \
           '--intervals '+join(outdir, "regs.mutect.txt")+ \
           '--input_file:normal {} --input_file:tumor {} --out {}'

# Returns the execution command for Somatic Sniper
def cmd_somaticsniper(ref):
    return "bam-somaticsniper -F vcf -f "+ref+" {} {} {}"

# Returns the execution command for Strelka
def cmd_strelka(ref, outdir):
    os.system('')
    outdir = join(outdir, "strelka")
    return 'configureStrelkaWorkflow.pl --ref '+ref+ \
           ' --config config/config.eland.ini ' \
           '--normal {} --tumor {} --out '+outdir+ \
           '; make -j 8 -C '+outdir+'; mv '+join(outdir, "results/passed.somatic.snvs.vcf")+ ' {}; ' \
           'rm -r '+outdir
