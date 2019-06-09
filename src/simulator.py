import sys
import getopt
import os
import getopt

PSEUDOMONAS_KOREENSIS = 198620
PANTONEA_AGGLOMERANS = 549
KLEBSIELLA_PNEUMONIAE = 573
NOT_CLASS = -10
version = "1.0"


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

f = 0
s = 0
t = 0
n = 0
percentage = 15/float(100)
outfile = "outputSimulator"



def help():
    print ("""\nScript for generating metagenomic data set in FASTA format.
        %s [-f|-s|-t|-n|-p|-o|-v|-h]
        -f, --first_bacteria: first bacteria
        -s, --second_bacteria: second bacteria
        -t, --third_bacteria: third baceteria
        -n, --not_bacteria: sequences that are not from bacteria
        -p, --percentage: percentage of sequences (default = 15)
        -o, --output_name: output name of files (default = "outputSimulator")
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
    print ("\nWrong input arguments. For more information read help: ")
    help()
    exit()

for opt,arg in opts:
    if opt in ['-f', '--first_bacteria']:
        bact1 = arg
        f = 1
    elif opt in ['-s', '--second_bacteria']:
        bact2 = arg
        s = 1
    elif opt in ['-t', '--third_bacteria']:
        bact3 = arg
        t = 1
    elif opt in ['-n', '--not_bacteria']:
        notbact = arg
        n = 1
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

if (not f or not s or not t or not n):
    print ("\nPlease input all required arguments. For more information read help: ")
    help()
    exit()

outpercentage = {}

infiles = [bact1, bact2, bact3, notbact]

pantoneaagglomerans = 0
pseudomonaskoreensis = 0
klebsiellapneumoniae = 0
notclass = 0

outnumseq = []
outlist = []
outlistId = {}

for infile in infiles:
    seq = list(SeqIO.parse(infile,"fasta"))
    totseq = len(seq)
    numseqOut = int(percentage*totseq)+1
    index = int(totseq/numseqOut)+1

    i = 1
    
    outseqfile = os.path.basename(infile)
    
    if (outseqfile == 'pantonea_agglomerans_r9_4_1.fasta') or (outseqfile == 'pantonea_agglomerans_miseq.fasta'):
        outseqid = PANTONEA_AGGLOMERANS
        pantoneaagglomerans = numseqOut
    elif (outseqfile == 'pseudomonas_koreensis_r9_5.fasta') or (outseqfile == 'pseudomonas_koreensis_miseq_2_300.fasta'):
        outseqid = PSEUDOMONAS_KOREENSIS
        pseudomonaskoreensis = numseqOut
    elif (outseqfile == 'klebsiella_pneumoniae_INF274_r9.fasta') or (outseqfile == 'klebsiella_pneumoniae_INF125_hiseq_2000.fasta'):
        outseqid = KLEBSIELLA_PNEUMONIAE
        klebsiellapneumoniae = numseqOut
    else:
        outseqid = NOT_CLASS
        notclass = numseqOut
    for j in range(numseqOut):
        if j != 0:
            i = j + index
        outseq = seq[i]
        outlist.append(outseq)
        outlistId[outseq.id] = outseqid

    f = open (outfile + "IDs.txt", "w")
    f.write(str(outlistId))
    f.close()
    
    
SeqIO.write(outlist, outfile + '.fasta', "fasta")

perc_pantoneaagglomerans = 100 * (pantoneaagglomerans/float(len(outlist)))
perc_pseudomonaskoreensis = 100 * (pseudomonaskoreensis/float(len(outlist)))
perc_klebsiellapneumoniae = 100 * (klebsiellapneumoniae/float(len(outlist)))
perc_notclass = 100 * (notclass/float(len(outlist)))

outpercentage[PANTONEA_AGGLOMERANS] = perc_pantoneaagglomerans
outpercentage[PSEUDOMONAS_KOREENSIS] = perc_pseudomonaskoreensis
outpercentage[KLEBSIELLA_PNEUMONIAE] = perc_klebsiellapneumoniae
outpercentage[NOT_CLASS] = perc_notclass

f = open (outfile + "Percentage.txt", "w")
f.write(str(outpercentage))
f.close()

exit()
