#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

fq_folder = sys.argv[1]
path = os.getcwd()
temp_folder = os.path.join(path, 'temp_folder1')
filtered_folder = os.path.join(path, 'filtered_folder1')
spoa_folder = os.path.join(path, 'spoa_folder1')

os.mkdir(temp_folder)
os.mkdir(filtered_folder)
os.mkdir(spoa_folder)

min_len = 0
max_len = 0

for filename in os.listdir(fq_folder):
	sys.stderr.write('\n\nPROCESSING FILENAME : %s' % filename)
	fq_filename = os.path.join(fq_folder, filename)
	if filename.lower().endswith(".fastq"):
		fname, fext = os.path.splitext(filename)
		tmp_filename = os.path.join(temp_folder, fname + '_tmp.fastq')
		flt_filename = os.path.join(filtered_folder, fname + '_filtered.fastq')
		consensus_filename = os.path.join(spoa_folder, fname + '.consensus')
		msa_filename = os.path.join(spoa_folder, fname + '.msa')
		matrix_filename = os.path.join(spoa_folder, fname + '.matrix')
        
        filterPath = os.path.join(os.getcwd(), 'samscripts/src/fastqfilter.py')
        string = 'python %s' % filterPath
        cmd = string + ' length_distribution %s' % fq_filename
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)
        mydict = {}        
        for line in output.splitlines():
            x=line.split()
            if len(x) < 5: continue
            mydict[x[0]] = x[4]
        
        key = max(mydict, key=int)
        value = mydict[key]
        min_len = max_len = int(value)
        min_len = min_len - 5
        max_len = max_len + 5      

        # Filter short sequences
        cmd = 'python /home/sanja/Desktop/zavrsniRad/za_skosier2/scripts/samscripts/src/fastqfilter.py minlen %d %s %s' % (min_len, fq_filename, tmp_filename)
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)
          
        # Filter long sequences
        cmd = 'python /home/sanja/Desktop/zavrsniRad/za_skosier2/scripts/samscripts/src/fastqfilter.py maxlen %d %s %s' % (max_len, tmp_filename, flt_filename)
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)

		# SPOA - consensus
        cmd = '/home/sanja/FER/doradenaSPOA/spoa/build/bin/spoa -r 0 -l 1 %s > %s' % (flt_filename, consensus_filename)
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)

        # SPOA - MSA
        cmd = '/home/sanja/FER/doradenaSPOA/spoa/build/bin/spoa -m 0 -n -1 -g -1 -r 1 -l 1 %s > %s' % (flt_filename, msa_filename)
        sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
        (status, output) = commands.getstatusoutput(cmd)
	
  

