from __future__ import print_function

import argparse, sys
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import SparseBitVect

from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo import draw

parser = argparse.ArgumentParser(prog='Clusterize.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('input', type=str, help='Input peptides.')
parser.add_argument('-o', '--output', default='out_tree.newick', type=str, help='Output newick tree.')
parser.add_argument('-t', '--output_ascii', default='out_tree_ascii.txt', type=str, help='Output newick tree.')
parser.add_argument('-a', '--annotation', type=str, help='Annotation')
args = parser.parse_args()

def get_similarity_from_bit_vectors(bit_vectors):
	simil = []
	count = 0
	n = 1
	for mol1 in range(len(bit_vectors)):
		simil.append([1-x for x in DataStructs.BulkTanimotoSimilarity(bit_vectors[mol1], bit_vectors[:mol1+1])])
		if count >= 100:
			#print(n * 100)
			n += 1
			count = 0
		count += 1
	return simil


def construct_distance_matrix(names, bit_vectors):
	return _DistanceMatrix(names, get_similarity_from_bit_vectors(bit_vectors))


def distance_matrix_to_tree(distance_matrix, method = "upgma"):
	constructor = DistanceTreeConstructor()
	if method == "upgma":
		return constructor.upgma(distance_matrix)
	else:
		return constructor.nj(distance_matrix)

def output_asci_tree(file_name, tree):
	with open(file_name, 'w') as handle:
		Phylo.draw_ascii(tree, file=handle)

def output_tree(file_name, tree):
	Phylo.write(tree, file_name, 'newick')

def read_patients_input(input_file):
	handle = open(input_file)
	header = handle.next().strip().split('\t')[1:]
	if 'Amount' in header:
		header.remove('Amount')
	#print(header)
	result = dict()
	for patient in header:
		result[patient] = dict()
	for line in handle:
		items = line.strip().split()
		peptide = items[0]
		for index, patient in enumerate(header):
			result[patient][peptide] = int(items[index + 1])

	handle.close()
	return result

def patient_params_to_vectors(patients_info_dict):
	result = dict()
	for patient in patients_info_dict:
		vector = SparseBitVect(len(patients_info_dict[patient]))
		index = 0
		for peptide_name in patients_info_dict[patient]:
			vector[index] = patients_info_dict[patient][peptide_name] > 0
			index += 1
		result[patient] = vector
	return result

class AnnotationBadFormat(Exception):

	def __init__(self, message):
		self.message = message

	def __str__(self):
		return 'AnnotationBadFormat ' + self.message


class DoubleAnotation(Exception):

	def __init__(self, message):
		self.message = message

	def __str__(self):
		return 'DoubleAnotation ' + self.message

def get_annotation(file_name):
	try:
		with open(file_name) as handle:
			header = handle.next().strip().split('\t')
			if 'ID' not in header:
				raise AnnotationBadFormat('ID not found')
			annotation = dict()
			for line in handle:
				entry = dict(zip(header, [x.strip() for x in line.split('\t')]))
				if len(entry['ID'].strip()) > 0:
					if entry['ID'] in annotation:
						#print(entry['ID'], annotation[entry['ID']])
						raise DoubleAnotation(entry['ID'])
					annotation[entry['ID']] = entry
			return annotation
	except Exception as e:
		print('There was error with parsing annotation ' + str(e))
	return None

def compose_annotation_id(sample):
	return sample['Diagnosis'] + ':' + sample['Diagnosis_add'] + ':' + sample['Folder'] + ':' + sample['FIO'] + ':' + sample['Gender'] + ':' + sample['Age'] + ':[' + sample['ID'] + ']'

patients_info = read_patients_input(args.input)
samples_peptides_vectors_dict = patient_params_to_vectors(patients_info)
# for patient in samples_peptides_vectors_dict:
# 	print(samples_peptides_vectors_dict[patient].ToBitString())

annotation = get_annotation(args.annotation)
# for item in annotation:
# 	print(compose_annotation_id(annotation[item]))

sample_names = []
sample_peptide_vectors = []
for sample in samples_peptides_vectors_dict:
	if annotation:
		if sample in annotation:
			sample_names.append(compose_annotation_id(annotation[sample]))
		else:
			print('There is no sample id ' + sample + ' in annotation')
			sample_names.append(sample)
	else:
		sample_names.append(sample)
	sample_peptide_vectors.append(samples_peptides_vectors_dict[sample])


# print(len(sample_names), len(sample_peptide_vectors))
distance_matrix = construct_distance_matrix(sample_names, sample_peptide_vectors)
tree = distance_matrix_to_tree(distance_matrix)
output_tree(args.output, tree)
output_asci_tree(args.output_ascii, tree)
