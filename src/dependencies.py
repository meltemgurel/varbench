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
