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

def draw_tree(tree):
	Phylo.draw_ascii(tree)
	# Phylo.draw_graphviz(tree)
	# pylab.show()


def read_patients_input(input_file):
	handle = open(input_file)
	header = handle.next().split('\t')[1:-1]
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

patients_info = read_patients_input(args.input)
patients_vectors = patient_params_to_vectors(patients_info)
# for patient in patients_vectors:
# 	print(patients_vectors[patient].ToBitString())
distance_matrix = construct_distance_matrix(patients_vectors.keys(), patients_vectors.values())
draw_tree(distance_matrix_to_tree(distance_matrix))
