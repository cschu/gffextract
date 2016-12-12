#!/usr/bin/env python
import sys
import csv

from ktoolu_io import readFasta

def reverseComplement(seq, alphabet='ACGT'):
    """
    Returns the reverse complement of nucleic acid seqence input.
    """
    # compl = dict(zip(alphabet, alphabet[::-1]))
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]

mRNAs = dict()
for row in csv.reader(sys.stdin, delimiter='\t'):
    if len(row) > 1:
        if row[2] == 'mRNA' and 'protein_coding' in row[8]:
            mRNAs[row[0]] = mRNAs.get(row[0], set())
            mRNAs[row[0]].add((int(row[3]), int(row[4]), row[6], row[8].split(';')[0][3:].split('.')[0]))

for _id, _seq in readFasta(sys.argv[1]):
    for start, end, strand, gid in mRNAs.get(_id[1:], set()):
        tid = '>' + gid + '_%i-%i:%c' % (start, end, strand)
        tseq = _seq[start - 1: end]
        if strand == '-':
            tseq = reverseComplement(tseq)
        print tid
        print tseq
        



