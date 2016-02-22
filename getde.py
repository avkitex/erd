#!/usr/bin/python

import re
import urllib.request
from bs4 import BeautifulSoup

minTreshold = 0.6
maxTreshold = 1.5
allRegulatedGenes = []
AFHandle = open('list.acc', 'r')
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
	def __init__(self):
		self.id = ''
		self.normal = dataset()
		self.cancer = dataset()
		self.age = 0
		self.upGenes = []
		self.downGenes = []
		

annotation = {}
#GSM494629	control	69	normal	97
#GSM494569	lung cancer	69	1A	97
for line in AFHandle:
	lspl = line.strip().split('\t')
	if lspl[4] not in annotation:
		annotation[lspl[4]] = patientAnn()
		annotation[lspl[4]].id = lspl[4]
		annotation[lspl[4]].age = lspl[2]
		annotation[lspl[4]].upGenes = []
		annotation[lspl[4]].downGenes = []
	if lspl[1] == 'control':
		if annotation[lspl[4]].normal.id != '':
			print('Double normal')
		annotation[lspl[4]].normal.id = lspl[0]
		annotation[lspl[4]].normal.genes = {}
	else:
		if annotation[lspl[4]].cancer.id != '':
			print('Double cancer')
		annotation[lspl[4]].cancer.id = lspl[0]
		annotation[lspl[4]].cancer.genes = {}
AFHandle.close()
outHandle = open('output.csv', 'w')
#for patient in annotation.values():
for pId in annotation:
	patient = annotation[pId]
	print(patient.id, patient.normal.id, patient.cancer.id)

	html1=urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=' + patient.normal.id).read()
	soup = BeautifulSoup(html1)
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

	html2=urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=' + patient.cancer.id).read()
	soup = BeautifulSoup(html2)
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
		diff = patient.cancer.genes[gene]/ patient.normal.genes[gene]
		if diff > maxTreshold:
			patient.upGenes.append(gene)
			allRegulatedGenes.append(gene)
		elif diff < minTreshold:
			patient.downGenes.append(gene)
			allRegulatedGenes.append(gene)
		#print(gene, patient.normal.genes[gene], patient.cancer.genes[gene], patient.cancer.genes[gene]/ patient.normal.genes[gene], sep = '\t', file=handle)
#handle.close()
print('gene', end = '', file=outHandle)
for pId in annotation:
	print('\t' + pId, end = '', file=outHandle)
print(file=outHandle)
for gene in allRegulatedGenes:
	print(gene, end = '', file=outHandle)
	for pId in annotation:
		patient = annotation[pId]
		print('\t', end = '', file=outHandle)
		if gene in patient.downGenes:
			print('d', end = '', file=outHandle)
		elif gene in patient.upGenes:
			print('u', end = '', file=outHandle)
		else:
			print('-', end = '', file=outHandle)
	print(file=outHandle)
	
	
	
	
	
	
outHandle.close()