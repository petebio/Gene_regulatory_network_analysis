#!/usr/bin/env python
import pybedtools as pb
import networkx as nx
import argparse
import json
import sys
import os

###################################################################################################################################################################################

# Read command line arguments
parser = argparse.ArgumentParser(description = 'Build a Gene Regulatory Network')
parser.add_argument('bed', type = str, help = 'Annotated BED file of peaks (with associated gene ID in 4th column) to create the GRN')
parser.add_argument('dir', type = str, help = 'Directory of BED files with the genomic coordinates for each TF motif to inlclude in the GRN')
parser.add_argument('exprs', type = str, help = 'Gene expression data file with gene ID (1st column) and gene expression value (e.g. FPKM; 2nd column)')
parser.add_argument('annot', type = str, help = 'Transcription Factor (TF) annotation file with gene ID and name of the motif it can bind to')
parser.add_argument('out', type = str, help = 'Output file')
parser.add_argument('-f', type = str, required = False, help = 'An optional BED file of DNaseI/ATAC-Seq footprints to filter motifs against')
parser.add_argument('-a', action = 'store_true', required = False, help = 'Include all genes in the network - not just transcription factor genes. Default: False')
parser.add_argument('-m', type = float, default = float(0), help = 'Minimum gene expression value for gene to include in the network. Default = 0')

args = parser.parse_args()

###################################################################################################################################################################################

# Read bed file of peaks from which to build the network
peaks = pb.BedTool(args.bed).sort()
sys.stdout.write('Read %d peaks from %s\n' % (peaks.count(), args.bed))

###################################################################################################################################################################################

# Read motif bed files
motif_positions = dict()

for motif_bed_file in os.listdir(args.dir):
	if not motif_bed_file.endswith('.bed'):
		continue # Skip file if it has the wrong file extension

	motif_id = motif_bed_file.replace('.bed', '')
	motif_positions[motif_id] = pb.BedTool('%s/%s' % (args.dir, motif_bed_file))

# Check if a footprint bed file is provided
# Only keep motifs that occur in footprints
if args.f:
	footprints = pb.BedTool(args.f).sort()

	for motif_id in motif_positions:
		motif_positions[motif_id] = motif_positions[motif_id].intersect(footprints, wa = True, u = True)

sys.stdout.write('Read motif positions for %d motifs\n' % len(motif_positions))

###################################################################################################################################################################################

# Read gene expression data
gene_expression = dict()

with open(args.exprs, 'r') as exprs:
	for num,line in enumerate(exprs):
		if num == 0:
			continue # Skip first line

		gene_id, value = line.strip().split('\t')
		value = float(value)

		if value >= args.m:
			gene_expression[gene_id] = value

# Read TF annotation data
# For each motif, find the TF gene with the highest expression. Use this as the reference source node for that motif
motif_associated_genes = []
motif_ref_gene = dict()

with open(args.annot, 'r') as annot:
	for num,line in enumerate(annot):
		if num == 0:
			continue # Skip first line

		gene_id, motif_id = line.strip().split('\t')

		if motif_id not in motif_positions:
			continue # Skip if motif is not included in the motif bed file directory

		# Check gene is expressed
		if gene_id in gene_expression:
			value = gene_expression[gene_id]
			motif_associated_genes.append(gene_id)

			# Check if motif has already been found
			# If not, add it to the motif_ref_gene dictionary with an empty gene ID and -Inf expression value
			if motif_id not in motif_ref_gene:
				motif_ref_gene[motif_id] = ['', float('-inf')]

			# Check if expression of gene is higher than the current motif reference gene. If yes, replace it
			if value > motif_ref_gene[motif_id][1]:
				motif_ref_gene[motif_id] = [gene_id, value]

sys.stdout.write('Read gene expression and annotation data for %d TF genes and %d TF families\n' % (len(motif_associated_genes), len(motif_ref_gene)))

###################################################################################################################################################################################

# Build the GRN
grn = nx.DiGraph()

for motif_id in motif_positions:
	if motif_id not in motif_ref_gene:
		continue # Skip if motif has no annotation data

	source_node, source_node_exprs = motif_ref_gene[motif_id]

	# Align motif to peaks and count the number of motifs associated with each gene
	peaks_with_motif = peaks.intersect(motif_positions[motif_id], wo = True)
	gene_motif_count = dict()

	for hit in peaks_with_motif:
		gene_id = hit[3]

		# Unless the full GRN has been requested, skip gene if not a TF
		if not args.a and gene_id not in motif_associated_genes:
			continue

		# Skip if gene is not expressed
		if gene_id not in gene_expression:
			continue

		# Check if gene already count and add if not
		if gene_id not in gene_motif_count:
			gene_motif_count[gene_id] = 0

		gene_motif_count[gene_id] += 1

	# Create nodes/edges
	for target_node in gene_motif_count:
		# Check if source and target nodes are already present in GRN. Add them if not
		if not grn.has_node(source_node):
			grn.add_node(source_node, expression = source_node_exprs)

		if not grn.has_node(target_node):
			target_node_exprs = gene_expression[target_node]
			grn.add_node(target_node, expression = target_node_exprs)

		grn.add_edge(source_node, target_node, count = gene_motif_count[target_node], source_motif = motif_id)

sys.stdout.write('Built network with %d nodes and %d edges\n' % (grn.number_of_nodes(), grn.number_of_edges()))

###################################################################################################################################################################################

# Write network as a Cytoscape JSON file
cytoscape = nx.cytoscape_data(grn)

# Ensure output file has correct file extension
outfile = args.out

if not outfile.endswith('.cyjs'):
	outfile = outfile + '.cyjs'

out = open(outfile, 'w')
json.dump(cytoscape, out)
out.close()

###################################################################################################################################################################################