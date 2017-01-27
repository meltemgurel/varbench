import sys
import pysam
import pandas

MAXDEPTH = int(sys.argv[4])

def get_allele_frequencies(bam, chr, pos):
    pos -= 1
    base_freq = {"A":0,"T":0,"C":0,"G":0}
    bam_file = pysam.AlignmentFile(bam, 'rb')
    for pileupcolumn in bam_file.pileup(chr, pos, max_depth = MAXDEPTH):
        found = []
        if pileupcolumn.pos == pos:
            for pileupread in pileupcolumn.pileups:
                if pileupread.query_position is not None \
                        and not pileupread.alignment.is_secondary \
                        and bin(pileupread.alignment.flag & 2048) != bin(2048) \
                        and not pileupread.alignment.mate_is_unmapped:
                    try:
                        found.index(pileupread.alignment.qname)
                    except ValueError :
                        base_freq[str(pileupread.alignment.seq[pileupread.query_position])] += 1
                        found.append(pileupread.alignment.qname)
    bam_file.close()
    total = sum(base_freq.itervalues(), 0.0)
    freqs = {k: float('%.3g' % (v / total)) for k, v in base_freq.iteritems()}
    print freqs
    return total, freqs

in_df  = pandas.DataFrame.from_csv(sys.argv[1], sep='\t',
                                   index_col=None, header=None)
out_df = pandas.DataFrame(index=range(0,len(in_df.index)),
                          columns=['CHR','POS','COV','REF','ALT','RVAF','AVAF'])

for index, row in in_df.iterrows():
    total, freqs = get_allele_frequencies(bam = sys.argv[2],
                                          chr = row[0],
                                          pos = row[1])
    out_df.loc[index] = [row[0], row[1], int(total),
                         row[2], row[3], '%f' % float('%.3g' % row[4]),
                         freqs[row[3]]]

out_df.to_csv(sys.argv[3], sep="\t", index=False)
