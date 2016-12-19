import sys
import pysam
import pandas

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

df = pandas.DataFrame.from_csv(sys.argv[1], sep='\t', index_col=None, header=None)
for index, row in df.iterrows():
    print(get_af(bam=sys.argv[2], chr=str(row[0]), pos=int(row[1])))
    print(get_af(bam=sys.argv[3], chr=str(row[0]), pos=int(row[1])))
    print("-----------------------------------------------------------")
