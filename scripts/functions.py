def check_bwa():
    p = subprocess.Popen(['bwa'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('Version:'):
            major, minor, sub = sline.strip().split()[1].split('.')
            sub = sub.split('-')[0]
            if int(major) >= 0 and int(minor) >= 7 and int(sub) >= 12:
                return True
    return False

def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('Version:'):
            major, minor = sline.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_bcftools():
    p = subprocess.Popen(['bcftools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('Version:'):
            major, minor = sline.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_wgsim():
    p = subprocess.Popen(['wgsim'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('Version:'):
            major, minor = sline.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 0 and int(minor) >= 2:
                return True
    return False


def check_velvet():
    p = subprocess.Popen(['velvetg'], stdout=subprocess.PIPE)
    for line in p.stdout:
        sline = line.decode("utf-8")
        if sline.startswith('Version'):
            major, minor = sline.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


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

def check_bamsurgeon():
    p = subprocess.Popen(['addsnv.py'], stderr=subprocess.PIPE)
    for line in p.stderr:
        sline = line.decode("utf-8")
        if sline.startswith('usage: addsnv.py'):
            return True
    return False

def check_python():
    return sys.hexversion >= 0x20702f0

def setup():
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

def get_name(x):
    return basename(os.path.splitext(x)[0])

def get_caller_cmd(caller, ref, muts, outdir):
    if caller == 'vardict':
        return cmd_vardict(ref, muts, outdir)
    elif caller == 'mutect':
        return cmd_mutect(ref, muts, outdir)
    elif caller == 'somaticsniper':
        return cmd_somaticsniper(ref)

def cmd_vardict(ref, muts, outdir):
    os.system("awk '{print $1\"\t\"$2\"\t\"$2}' "+muts+ \
              ' | tail -n +2 > '+join(outdir, "regs.vardict.bed"))
    return 'vardict-java -G '+ref+' -f 0.001 -N "{}|{}" -b "{}|{}" ' \
           '-c 1 -S 2 -E 3 '+join(outdir, "regs.vardict.bed")+' | scripts/testsomatic.R | scripts/var2vcf_paired.pl ' \
           '-N "{}|{}" -f 0.001 > {}'

def cmd_mutect(ref, muts, outdir):
    os.system('awk {print $1":"$2"-"$2} '+muts+
              ' | tail -n +2 > '+join(outdir, "regs.mutect.txt"))
    return 'java -Xmx2g -jar mutect/mutect.jar --analysis_type MuTect --reference_sequence '+ref+ \
           '--cosmic mutect/b37_cosmic_v54_120711.vcf --dbsnp mutect/dbsnp_132_b37.leftAligned.vcf.gz ' \
           '--intervals '+join(outdir, "regs.mutect.txt")+ \
           '--input_file:normal {} --input_file:tumor {} --out {}'

def cmd_somaticsniper(ref):
    return "bam-somaticsniper -F vcf -f "+ref+" {} {} {}"
