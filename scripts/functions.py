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

def get_af(bam, chr, pos):
    pos -= 1
    baseFreq = {"A":0,"T":0,"C":0,"G":0}
    bamFile = pysam.AlignmentFile(bam, 'rb')
    for pileupcolumn in bamFile.pileup(chr, pos):
        found = []
        if pileupcolumn.pos == pos:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        found.index(pileupread.alignment.qname)
                    except ValueError :
                        baseFreq[str(pileupread.alignment.seq[pileupread.query_position])] += 1
                        found.append(pileupread.alignment.qname)
    bamFile.close()
    return baseFreq
