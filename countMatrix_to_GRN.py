#!/usr/bin/env python
import networkx as nx
import argparse
import json
import sys

###################################################################################################################################################################################

# Read command-line arguments
parser = argparse.ArgumentParser(description = 'Convert a motif count matrix to a Gene Regulatory Network')
parser.add_argument('matrix', type = str, help = 'Motif count matrix file')
parser.add_argument('exprs', type = str, help = 'Gene expression data file')
parser.add_argument('tfanno', type = str, help = 'Transcription Factor annotation file')
parser.add_argument('out', type = str, help = 'Output file')

args = parser.parse_args()

###################################################################################################################################################################################

# Read TF annotation data
tf_annotations = dict()
motif_families = set()

with open(args.tfanno, 'r') as tfanno:
	for num,line in enumerate(tfanno):
		if num == 0:
			continue # Skip header line

		gene_id, motif_id = line.strip().split('\t')

		tf_annotations[gene_id] = motif_id
		motif_families.add(motif_id)

sys.stdout.write('Read TF annotation data for %d genes from %d motif families\n' % (len(tf_annotations), len(motif_families)))

###################################################################################################################################################################################

# Read gene expression data
# For each TF motif family, find the highest expressed member
# This will be used as the source node for that motif
gene_expression = dict()
motif_family_ref = dict()

with open(args.exprs, 'r') as exprs:
	for num,line in enumerate(exprs):
		if num == 0:
			continue # Skip header line

		gene_id, exprs_value = line.strip().split('\t')
		exprs_value = float(exprs_value)

		gene_expression[gene_id] = exprs_value

		# Check if gene is a transcription factor
		if gene_id in tf_annotations:
			motif_id = tf_annotations[gene_id]

			# Check if motif already found and add to motif_family_ref
			if motif_id not in motif_family_ref:
				motif_family_ref[motif_id] = ['', float('-Inf')]

			if exprs_value > motif_family_ref[motif_id][1]:
				motif_family_ref[motif_id] = [gene_id, exprs_value]

sys.stdout.write('Read gene expression data for %d genes\n' % len(gene_expression))

###################################################################################################################################################################################

# Read count matrix and create GRN
grn = nx.DiGraph()

with open(args.matrix, 'r') as matrix:
	for num,row in enumerate(matrix):
		columns = row.strip().split('\t')

		if num == 0:
			motifs = columns
		else:
			target_node = columns[0]

			# Check if gene is expressed and get expression value
			if target_node in gene_expression:
				exprs_value = gene_expression[target_node]
			else:
				continue # Skip non-expressed genes

			for i in range(1, len(columns)):
				motif_id = motifs[i]

				if motif_id in motif_family_ref:
					source_node = motif_family_ref[motif_id][0]
				else:
					continue

				count = int(columns[i])

				if count == 0:
					continue # Skip edges with 0 motifs

				# Check if source and target nodes are already present in network
				# Add them if they not
				if not grn.has_node(source_node):
					source_node_exprs = gene_expression[source_node]
					grn.add_node(source_node, expression = source_node_exprs)

				if not grn.has_node(target_node):
					grn.add_node(target_node, expression = exprs_value)

				# Create edge
				grn.add_edge(source_node, target_node, count = count)

sys.stdout.write('Created GRN with %d nodes and %d edges\n' % (grn.number_of_nodes(), grn.number_of_edges()))

###################################################################################################################################################################################

# Write network as a Cytoscape JSON file
cytoscape = nx.cytoscape_data(grn)

# Ensure output file has the correct file extension
outfile = args.out

if not outfile.endswith('.cyjs'):
	outfile = outfile + '.cyjs'

out = open(outfile, 'w')
json.dump(cytoscape, out)
out.close()

###################################################################################################################################################################################