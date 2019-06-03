import sys
import getopt
import os
import csv
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

argv = sys.argv[1:]
krakenfiletitle = argv[0]
clarkfiletitle = argv[1]
percentagefletitle = argv[2]
idsfiletitle = argv[3]
outputfile = argv[4]

percentagedict = {}
seqtaxoutputdict = {}
seqtaxkrakendict = {}
seqtaxclarkdict = {}

numseqkraken = 0
numseqclark = 0

numofclassifiedkraken = 0
numofclassifiedclark = 0

numofpantoaggloout = 0
numofpseudokoreensout = 0
numofklebipneumout = 0

numofpantoagglokraken = 0
numofpseudokoreenskraken = 0
numofklebipneumkraken = 0

numofpantoaggloclark = 0
numofpseudokoreensclark = 0
numofklebipneumclark = 0

numcorrectseqclark = 0
numcorrectseqkraken = 0

percentagefile = open(percentagefletitle, "r")
for line in percentagefile:
    sequences = line.split(", ")
    for seq in sequences:
        parts = seq.split(": ")
        if parts[0].startswith("{"):
            id = parts[0][1:]
        else:
            id = parts[0]
        if parts[1].endswith("}"):
            percentage = parts[1][:-1]
        else:
            percentage = parts[1]
        percentagedict[id] = percentage

idsfile = open(idsfiletitle, "r")
line = idsfile.read()
sequences = line.split(", ")
for seq in sequences:
    parts = seq.split(": ")
    if parts[0].startswith("{"):
        seqid = parts[0][2:-1]
    else:
        seqid = parts[0][1:-1]
    if parts[1].endswith("}"):
        taxid = parts[1][:-1]
    else:
        taxid = parts[1]
    if taxid == '198620':
        numofpseudokoreensout += 1
    elif taxid == '549':
        numofpantoaggloout += 1
    elif taxid == '573':
        numofklebipneumout += 1
    seqtaxoutputdict[seqid] = taxid

taxonomyID = ""
with open(krakenfiletitle) as csv_file:
    csvreader = csv.reader(csv_file, delimiter = '\t')
    linecount = 0
    for row in csvreader:
        if linecount == 0:
            print "Empty row."
        else:
            numseqkraken += 1
            isclassified = row[0]
            sequenceID = row[1]
            taxonomyID = row[2]
            sequencelength = row[3]
            if isclassified == 'C':
                numofclassifiedkraken += 1
            seqtaxkrakendict[sequenceID] = (taxonomyID, isclassified)
        linecount += 1
        if taxonomyID == '198620':
            numofpseudokoreenskraken += 1
        elif taxonomyID == '549':
            numofpantoagglokraken += 1
        elif taxonomyID == '573':
            numofklebipneumkraken += 1

clarkisclasified = 0
with open(clarkfiletitle) as csv_file:
    csvreader = csv.reader(csv_file, delimiter = ',')
    linecount = 0
    for row in csvreader:
        if linecount == 0:
            print "Empty row."
        else:
            numseqclark += 1
            sequenceID = row[0]
            sequencelength = row[1]
            taxonomyID = row[2]
            if taxonomyID != 'NA':
                numofclassifiedclark += 1
                clarkisclassified = 1
            else:
                clarkisclassified = 0
            seqtaxclarkdict[sequenceID] = (taxonomyID, clarkisclassified)
        linecount += 1
        if taxonomyID == '198620':
            numofpseudokoreensclark += 1
        elif taxonomyID == '549':
            numofpantoaggloclark += 1
        elif taxonomyID == '573':
            numofklebipneumclark += 1

numseqclark -= 1

numofoutseq = len(seqtaxoutputdict)
numsameseqkraken = 0
numsameseqclark = 0
#i = 0
for seq in seqtaxoutputdict:
    outtaxid = seqtaxoutputdict[seq]
    if seq in seqtaxkrakendict:
        krakentaxid, isclassifiedkraken = seqtaxkrakendict[seq]
        if outtaxid == krakentaxid:
            numsameseqkraken += 1
            if isclassifiedkraken == 'C':
                numcorrectseqkraken += 1
    if seq in seqtaxclarkdict:
        clarktaxid, isclassifiedclark = seqtaxclarkdict[seq]
        if outtaxid == clarktaxid:
            numsameseqclark += 1
            if isclassifiedclark:
                numcorrectseqclark += 1
#print i
#i += 1

percentageofsameseqkraken = (numsameseqkraken/float(numofoutseq)) * 100
percentageofsameseqclark = (numsameseqclark/float(numofoutseq)) * 100

percentageofclassifiedkraken = (numofclassifiedkraken/float(numofoutseq)) * 100
percentageofclassifiedclark = (numofclassifiedclark/float(numofoutseq)) * 100

numofunclassifiedkraken = numseqkraken - numofclassifiedkraken
numofunclassifiedclark = numseqclark - numofclassifiedclark

numoffalseclassifiedkraken = numofclassifiedkraken - numcorrectseqkraken
numoffalseclassifiedclark = numofclassifiedclark - numcorrectseqclark

percentageoffalseclassifiedkraken = (numoffalseclassifiedkraken/float(numofoutseq)) * 100
percentageoffalseclassifiedclark = (numoffalseclassifiedclark/float(numofoutseq)) * 100

percentageofunclassifiedkraken = (numofunclassifiedkraken/float(numofoutseq)) * 100
percentageofunclassifiedclark = (numofunclassifiedclark/float(numofoutseq)) * 100

percentageofcorrectinclassifiedkraken = (numcorrectseqkraken/float(numofclassifiedkraken)) * 100
percentageofcorrectinclassifiedclark = (numcorrectseqclark/float(numofclassifiedclark)) * 100

percentageoffalseinclassifiedkraken = (numoffalseclassifiedkraken/float(numofclassifiedkraken)) * 100
percentageoffalseinclassifiedclark = (numoffalseclassifiedclark/float(numofclassifiedclark)) * 100

differencepseudokoreenskraken = abs(numofpseudokoreensout - numofpseudokoreenskraken)
differencepseudokoreensclark = abs(numofpseudokoreensclark - numofpseudokoreensout)

errorpseudokoreenskraken = differencepseudokoreenskraken/(float(numofpseudokoreensout))
errorpseudokoreensclark = differencepseudokoreensclark/(float(numofpseudokoreensout))

differenceklebipneumkraken = abs(numofklebipneumkraken - numofklebipneumout)
differenceklebipneumclark = abs(numofklebipneumclark - numofklebipneumout)

errorklebipneumkraken= differenceklebipneumkraken/(float(numofklebipneumout))
errorklebipneumclark = differenceklebipneumclark/(float(numofklebipneumout))

differencepantoagglokraken = abs(numofpantoagglokraken - numofpantoaggloout)
differencepantoaggloclark = abs(numofpantoaggloclark - numofpantoaggloout)

errorpantoagglokraken = differencepantoagglokraken/(float(numofpantoaggloout))
errorpantoaggloclark = differencepantoaggloclark/(float(numofpantoaggloclark))


f = open(outputfile + ".txt", "w")
f.write("KRAKEN:"+ "\n")
f.write("-------------------------------------------" + "\n")

f.write("NUMBER OF SEQUENCES: " + str(numseqkraken) + "\n")
f.write("NUMBER OF CLASSIFIED SEQUENCES: " + str(numofclassifiedkraken) + "\n")
f.write("PERCENTAGE OF CLASSIFIED SEQUENCES: " + str(percentageofclassifiedkraken) + "\n")
f.write("NUMBER OF UNCLASSIFIED SEQUENCES: " + str(numofunclassifiedkraken) + "\n")
f.write("PERCENTAGE OF UNCLASSIFIED SEQUENCES: " + str(percentageofunclassifiedkraken) + "\n")
f.write("NUMBER OF CORRECT CLASSIFIED SEQUENCES: " + str(numcorrectseqkraken) + "\n")
f.write("PERCENTAGE OF CORRECT CLASSIFIED SEQUENCES: " + str(percentageofsameseqkraken) + "\n")
f.write("NUMBER OF FALSE CLASSIFIED SEQUENCES: " + str(numoffalseclassifiedkraken) + "\n")
f.write("PERCENTAGE OF FALSE CLASSIFIED SEQUENCES: " + str(percentageoffalseclassifiedkraken) + "\n")
f.write("PERCENTAGE OF CORRECT CLASSIFIED SEQUENCES IN SET OF CLASSIFIED SEQUENCES: " + str(percentageofcorrectinclassifiedkraken) + "\n")
f.write("PERCENTAGE OF FALSE CLASSIFIED SEQUENCES IN SET OF CLASSIFIED SEQUENCES: " + str(percentageoffalseinclassifiedkraken) + "\n")
f.write("NUMBER OF PSEUDOMONAS KOREENSIS SEQUENCES: " + str(numofpseudokoreenskraken) + "\n")
f.write("PSEUDOMONAS KOREENSIS ERROR: " + str(errorpseudokoreenskraken) + "\n")
f.write("NUMBER OF KLEBSIELLA PNEUMONIAE SEQUENCES: " + str(numofklebipneumkraken) + "\n")
f.write("KLEBSIELLA PNEUMONIAE ERROR: " + str(errorklebipneumkraken) + "\n")
f.write("NUMBER OF PANTONEA AGGLOMERANS SEQUENCES: " + str(numofpantoagglokraken) + "\n")
f.write("PANTONEA AGGLOMERANS ERROR: " + str(errorpantoagglokraken) + "\n")

f.write("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" + "\n")
f.write("CLARK:"+ "\n")
f.write("-------------------------------------------"+ "\n")
f.write("NUMBER OF SEQUENCES: " + str(numseqclark) + "\n")
f.write("NUMBER OF CLASSIFIED SEQUENCES: " + str(numofclassifiedclark) + "\n")
f.write("PERCENTAGE OF CLASSIFIED SEQUENCES: " + str(percentageofclassifiedclark) + "\n")
f.write("NUMBER OF UNCLASSIFIED SEQUENCES: " + str(numofunclassifiedclark) + "\n")
f.write("PERCENTAGE OF UNCLASSIFIED SEQUENCES: " + str(percentageofunclassifiedclark) + "\n")
f.write("NUMBER OF CORRECT CLASSIFIED SEQUENCES: " + str(numcorrectseqclark) + "\n")
f.write("PERCENTAGE OF CORRECT CLASSIFIED SEQUENCES: " + str(percentageofsameseqclark) + "%\n")
f.write("NUMBER OF FALSE CLASSIFIED SEQUENCES: " + str(numoffalseclassifiedclark) + "\n")
f.write("PERCENTAGE OF FALSE CLASSIFIED SEQUENCES: " + str(percentageoffalseclassifiedclark) + "\n")
f.write("PERCENTAGE OF CORRECT CLASSIFIED SEQUENCES IN SET OF CLASSIFIED SEQUENCES: " + str(percentageofcorrectinclassifiedclark) + "\n")
f.write("PERCENTAGE OF FALSE CLASSIFIED SEQUENCES IN SET OF CLASSIFIED SEQUENCES: " + str(percentageoffalseinclassifiedclark) + "\n")
f.write("NUMBER OF PSEUDOMONAS KOREENSIS SEQUENCES: " + str(numofpseudokoreensclark) + "\n")
f.write("PSEUDOMONAS KOREENSIS ERROR: " + str(errorpseudokoreensclark) + "\n")
f.write("NUMBER OF KLEBSIELLA PNEUMONIAE SEQUENCES: " + str(numofklebipneumclark) + "\n")
f.write("KLEBSIELLA PNEUMONIAE ERROR: " + str(errorklebipneumclark) + "\n")
f.write("NUMBER OF PANTONEA AGGLOMERANS SEQUENCES: " + str(numofpantoaggloclark) + "\n")
f.write("PANTONEA AGGLOMERANS ERROR: " + str(errorpantoaggloclark) + "\n")

f.close()





