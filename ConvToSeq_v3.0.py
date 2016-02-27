#!/usr/bin/python

# Python2.7 ConvToSeq_v3.0.py - downloads fasta-files from NCBI using txt-file with GenBank Protein GI and Affymetrix ID as input-file.
# Author: Aksenova Marina
#Last modification: 23.02.2016, 20:08

import urllib.request
import os

inputName='bioDBnet_db2db_160222071622_1532789246.txt'
directory = 'ConvToSeq'
if not os.path.exists(directory):
    os.makedirs(directory)
inputfile=open(inputName, 'r')
firstList=inputfile.readlines()
giList=[]
idList=[]
for line in firstList:
	n=line.find('at')
	idList.append(line[:n+2])
	giList.append(line[n+3:-2].split('; '))

print ('Downloading *.fasta files with proteins/peptides for:')

fastafile=open(directory+'/base.fasta', 'w')
for i in range(len(giList)):
	for j in range (len(giList[i])):
		if giList[i][j]=='':
			pass
		else:
			url='http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=fasta&sort=&val='+str(giList[i][j])+'&'
			newfasta=urllib.request.urlopen(url).read().decode('utf-8')
			fstr = newfasta.strip().split('\n')
			
			print('>aff|' + idList[i] + '|' + fstr[0][1:], file=fastafile)
			print('\n'.join(fstr[1:] + ['']), file=fastafile)
			
			print(giList[i][j], end = '') 
fastafile.close()

print ('Finished')



# To get input-file: go to https://biodbnet-abcc.ncifcrf.gov/db/db2db.php, then choose input ("Affy ID"), output ("GenBank Protein GI")
# Then put your list of Affymetrix ID's to "ID list" and push "Submit"
# Then in new window with results download comlpete results ("Complete results can be downloaded here") and in this result-file delete the 1st string
