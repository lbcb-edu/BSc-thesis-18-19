import sys
import getopt
import os
PSEUDOMONAS_KOREENSIS = 198620
PANTONEA_AGGLOMERANS = 549
KLEBSIELLA_PNEUMONIAE = 573
NOT_CLASS = -10

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

argv = sys.argv[1:]
bact1 = argv[0]
bact2 = argv[1]
bact3 = argv[2]
notbact = argv[3]
percentage = int(argv[4])/float(100)
outfile = argv[5]

outprec = {}

infiles = [bact1, bact2, bact3, notbact]

pantoagglo = 0
pseudokoreens = 0
klebipneum = 0
notclass = 0

outnumseq = []
outlist = []
outlistId = {}
l = 0

for infile in infiles:
    seq = list(SeqIO.parse(infile,"fasta"))
    totseq = len(seq)
    numseqOut = int(percentage*totseq)+1
    index = int(totseq/numseqOut)+1

    i = 1
    m = 0
    l = 0
    
    outseqfile = infile
    
    if (outseqfile == 'pantonea_agglomerans_r9_4_1.fasta') or (outseqfile == 'pantonea_agglomerans_miseq.fasta'):
        outseqid = PANTONEA_AGGLOMERANS
        pantoagglo = numseqOut
    elif (outseqfile == 'pseudomonas_koreensis_r9_5.fasta') or (outseqfile == 'pseudomonas_koreensis_miseq_2_300.fasta'):
        outseqid = PSEUDOMONAS_KOREENSIS
        pseudokoreens = numseqOut
    elif (outseqfile == 'klebsiella_pneumoniae_INF274_r9.fasta') or (outseqfile == 'klebsiella_pneumoniae_INF125_hiseq_2000.fasta'):
        outseqid = KLEBSIELLA_PNEUMONIAE
        klebipneum = numseqOut
    else:
        outseqid = NOT_CLASS
        notclass = numseqOut
    for j in range(numseqOut):
        l += 1
        if j != 0:
            i = j + index
        outseq = seq[i]
        outlist.append(outseq)
        outlistId[outseq.id] = outseqid

    f = open (outfile + "IDs.txt", "w")
    f.write(str(outlistId))
    f.close()
    
    
SeqIO.write(outlist, outfile + '.fasta', "fasta")

prec_pantoagglo = 100 * (pantoagglo/float(len(outlist)))
prec_pseudokoreens = 100 * (pseudokoreens/float(len(outlist)))
prec_klebiprenum = 100 * (klebipneum/float(len(outlist)))
prec_notclass = 100 * (notclass/float(len(outlist)))

outprec[PANTONEA_AGGLOMERANS] = prec_pantoagglo
outprec[PSEUDOMONAS_KOREENSIS] = prec_pseudokoreens
outprec[KLEBSIELLA_PNEUMONIAE] = prec_klebiprenum
outprec[NOT_CLASS] = prec_notclass

f = open (outfile + "Percentage.txt", "w")
f.write(str(outprec))
f.close()

exit()
