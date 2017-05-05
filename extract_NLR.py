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

sequences = dict()
for row in csv.reader(sys.stdin, delimiter='\t'):
    if len(row) > 1:
        if row[2] == 'NBSLRR':
            sequences[row[0]] = sequences.get(row[0], set())
            sequences[row[0]].add((int(row[3]), int(row[4]), row[6], row[8].strip().split('=')[1]))

for key in sequences:
    print(sequences[key])

for _id, _seq in readFasta(sys.argv[1]):
    for start, end, strand, gid in mRNAs.get(_id[1:], set()):
        tid = '>' + gid + '_%i-%i:%c' % (start, end, strand)
        tseq = _seq[start - 1: end]
        if strand == '-':
            tseq = reverseComplement(tseq)
        print(tid + '\n' + tseq)
