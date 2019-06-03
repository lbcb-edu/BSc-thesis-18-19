import sys
import getopt
import os
import getopt

PSEUDOMONAS_KOREENSIS = 198620
PANTONEA_AGGLOMERANS = 549
KLEBSIELLA_PNEUMONIAE = 573
NOT_CLASS = -10
version = "1.0"
numOfRequiredArg = 6


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def help():
    print ("""usage: Script for generating metagenomic data set in FASTA format.
        %s [-f|-s|-t|-n|-p|-o|-v|-h]
        -f, --first_bacteria: first bacteria
        -s, --second_bacteria: second bacteria
        -t, --third_bacteria: third baceteria
        -n, --not_bacteria: sequences that are not from bacteria
        -p, --percentage: percentage of sequences
        -o, --output_name: output name of files
        -v, --version: version
        -h, --help: help
        """)

argv = sys.argv[1:]
opts = []
args = []

try:
    opts, args = getopt.getopt(argv, "f:s:t:n:p:o:vh", ["first_bacteria=", "second_bacteria=", "third_bacteria=", "not_bacteria=", "percentage=", "output_name=", "version", "help"])
except getopt.GetoptError as err:
    print(err)
    help()

if opts == []:
    print ("Wrong input arguments. Number of required arguments: " + str(numOfRequiredArg) + " For more information read help: ")
    help()
    exit()

for opt,arg in opts:
    if opt in ['-f', '--first_bacteria']:
        bact1 = arg
    elif opt in ['-s', '--second_bacteria']:
        bact2 = arg
    elif opt in ['-t', '--third_bacteria']:
        bact3 = arg
    elif opt in ['-n', '--not_bacteria']:
        notbact = arg
    elif opt in ['-p', '--percentage']:
        percentage = int(arg)/float(100)
    elif opt in ['-o', '--output_name']:
        outfile = arg
    elif opt in ['-v', '--version']:
        print ("version: " + str(version))
        exit()
    elif opt in ['-h', '--help']:
        help()
        exit()

if opts != numOfRequiredArg:
    difference = numOfRequiredArg - len(opts)
    print ("Number of missing arguments: " + str(difference) + ". Please input all arguments. For more information read help")
    help()
    exit()

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
