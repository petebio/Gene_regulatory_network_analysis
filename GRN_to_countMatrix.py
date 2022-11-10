#!/usr/bin/env python
import networkx as nx
import argparse
import json
import sys

###################################################################################################################################################################################

# Read command line arguments
parser = argparse.ArgumentParser(description = 'Convert a GRN to a motif count matrix')
parser.add_argument('cyjs', type = str, help = 'Cytoscape JSON file of the GRN')
parser.add_argument('out', type = str, help = 'Output file')

args = parser.parse_args()

###################################################################################################################################################################################

# Read network
grn_json_file = open(args.cyjs, 'r')
grn_json = json.load(grn_json_file)
grn_json_file.close()

# Convert JSON to networkx
grn = nx.cytoscape_graph(grn_json)
sys.stdout.write('Read network with %d nodes and %d edges\n' % (grn.number_of_nodes(), grn.number_of_edges()))

###################################################################################################################################################################################

# Create count matrix
motifs = set()
gene_motif_count = dict()

for source, target, metadata in grn.edges(data = True):
	motif_id = metadata['source_motif']
	count = int(metadata['count'])

	motifs.add(motif_id)

	# Check if target gene is already in count table
	if target not in gene_motif_count:
		gene_motif_count[target] = dict()

	gene_motif_count[target][motif_id] = count

###################################################################################################################################################################################

# Write output file
out = open(args.out, 'w')

# Write file header
motifs = sorted(motifs)
out.write('Motif\t%s\n' % '\t'.join(motifs))

for gene_id in sorted(gene_motif_count):
	out.write(gene_id)

	for motif_id in motifs:
		if motif_id in gene_motif_count[gene_id]:
			count = gene_motif_count[gene_id][motif_id]
		else:
			count = 0

		out.write('\t%d' % count)

	out.write('\n')

out.close()

###################################################################################################################################################################################