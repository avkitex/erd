#!/usr/bin/python

# python3 getde.py 

import re 
filename='masout1'
fileH = open(filename + '.txt.xls')
betweenDots=True

header = fileH.readline()
assert header[0] == '#'
headerArr = header[1:].strip().split('\t')
curEntry= {}
#Filename	Hit Number	Protein hit ID	Gene Symbol	Gene name	Query	Observed	Mr(expt)	Mr(calc)	ppm	Miss	Score	Expect	Peptide	RED
#F001405	Hit1	IPI00220327	KRT1	Keratin, type II cytoskeletal 1
outHandle = open(filename + '.fasta', 'w')
for line in fileH:
	lineArr = line.strip().split('\t')
	for tag in lineArr:
		tag = tag.strip()
	entry = dict(zip(headerArr, lineArr))
	if len(entry['Hit Number']):
		curEntry = entry
	else:
		if 'RED' not in entry:
			entry['RED'] = '-'
		print('>mascot|' + filename + '|' + curEntry['Filename'] + '|' + curEntry['Hit Number'] + '|' + curEntry['Protein hit ID'] + '|' + curEntry['Gene Symbol'] + '|' + entry['Query'] + '|' + entry['RED'], file=outHandle)
		peptide=entry['Peptide'] 
		if betweenDots:
			peptide = peptide[2:-2]
		else:
			if peptide[0] == '-':
				peptide = peptide[2:]
			else:
				peptide = peptide[0] + peptide[2:]
			if peptide[-1] == '-':
				peptide = peptide[:-2]
			else:
				peptide = peptide[:-2] + peptide[-1]
		print(peptide, file=outHandle)
file=outHandle.close()
fileH.close()