#!/usr/bin/python

# sudo apt-get install python3-bs4 python3-urllib
# python3 getde.py 

import re
import urllib.request
from bs4 import BeautifulSoup

minTreshold = 0.35
maxTreshold = 2.5
AFHandle = open('list.acc', 'r')
outHandle = open('output_' + str(minTreshold) + '_' + str(maxTreshold) + '.csv', 'w')



allRegulatedGenes = []

def visible(element):
	if element.parent.name in ['style', 'script', '[document]', 'head', 'title']:
		return False
	elif re.match('<!--.*-->', str(element)):
		return False
	return True

class dataset():
	id = ''
	genes = {}
	def __init__(self):
		self.id = ''
		self.genes = {}
	
class patientAnn():
	id = ''
	normal = None
	cancer = None
	upGenes = []
	downGenes = []
	age = 0
	stage = ''
	def __init__(self):
		self.id = ''
		self.normal = dataset()
		self.cancer = dataset()
		self.age = 0
		self.upGenes = []
		self.downGenes = []
		self.stage = ''
		

annotation = {}
#GSM494629	control	69	normal	97
#GSM494569	lung cancer	69	1A	97
for line in AFHandle:
	lspl = line.strip().split('\t')
	pId = int(lspl[4].strip())
	if pId not in annotation:
		annotation[pId] = patientAnn()
		annotation[pId].id = pId
		annotation[pId].age = lspl[2]
		annotation[pId].upGenes = []
		annotation[pId].downGenes = []
	if lspl[1] == 'control':
		if annotation[pId].normal.id != '':
			print('Double normal')
		annotation[pId].normal.id = lspl[0]
		annotation[pId].normal.genes = {}
	else:
		if annotation[pId].cancer.id != '':
			print('Double cancer')
		annotation[pId].cancer.id = lspl[0]
		annotation[pId].cancer.genes = {}
		annotation[pId].stage = lspl[3]
AFHandle.close()
#for patient in annotation.values():
sortedPatients = sorted(annotation.keys())
for pId in sortedPatients:
	patient = annotation[pId]
	print(patient.id, patient.normal.id, patient.cancer.id)

	html=urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=' + patient.normal.id).read()
	soup = BeautifulSoup(html)
	texts = soup.findAll(text=True)
	visible_texts = filter(visible, texts)
	for file in visible_texts:
		for line in file.split('\n'):
			spLine = line.strip()
			if len(spLine):
				spLine = spLine.split()
				if len(spLine) == 2 and spLine[0][0] != '#':
					gene = spLine[0].strip()
					expression = float(spLine[1].strip())
					if gene in patient.normal.genes:
						print('Double gene value1')
					patient.normal.genes[gene] = expression
					#print(gene, expression)

	html=urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=' + patient.cancer.id).read()
	soup = BeautifulSoup(html)
	texts = soup.findAll(text=True)
	visible_texts = filter(visible, texts)
	for file in visible_texts:
		for line in file.split('\n'):
			spLine = line.strip()
			if len(spLine):
				spLine = spLine.split()
				if len(spLine) == 2 and spLine[0][0] != '#':
					gene = spLine[0].strip()
					expression = float(spLine[1].strip())
					if gene in patient.cancer.genes:
						print('Double gene value2')
					patient.cancer.genes[gene] = expression
	#handle = open(patient.id + '_pat.dif', 'w')
	for gene in patient.normal.genes:
		if gene not in patient.cancer.genes:
			print('No ' + gene + ' in cancer')
			continue
		diff = patient.cancer.genes[gene] / patient.normal.genes[gene]
		if diff > maxTreshold:
			patient.upGenes.append(gene)
			if gene not in allRegulatedGenes:
				allRegulatedGenes.append(gene)
		elif diff < minTreshold:
			patient.downGenes.append(gene)
			if gene not in allRegulatedGenes:
				allRegulatedGenes.append(gene)
		#print(gene, patient.normal.genes[gene], patient.cancer.genes[gene], patient.cancer.genes[gene]/ patient.normal.genes[gene], sep = '\t', file=handle)
	#handle.close()
print('gene', end = '', file=outHandle)
for pId in sortedPatients:
	print('\t', pId, sep = '', end = '', file=outHandle)
print('\tcountU\tcountD', file=outHandle)
for pId in sortedPatients:
	print('\t' + annotation[pId].age, end = '', file=outHandle)
print('\t\t', file=outHandle)
for pId in sortedPatients:
	print('\t' + annotation[pId].stage, end = '', file=outHandle)
print('\t\t', file=outHandle)
allRegulatedGenes.sort()
for gene in allRegulatedGenes:
	countU = 0
	countD = 0
	print(gene, end = '', file=outHandle)
	for pId in sortedPatients:
		patient = annotation[pId]
		print('\t', end = '', file=outHandle)
		if gene in patient.downGenes:
			countD += 1
			print('d', end = '', file=outHandle)
		elif gene in patient.upGenes:
			countU += 1
			print('u', end = '', file=outHandle)
		else:
			print('-', end = '', file=outHandle)
	print('\t', countU, '\t', countD, sep = '', file=outHandle)
	
outHandle.close()