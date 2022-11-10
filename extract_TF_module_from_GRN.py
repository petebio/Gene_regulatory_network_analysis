#!/usr/bin/env python
import networkx as nx
import argparse
import json
import sys

###################################################################################################################################################################################

# Read command-line arguments
parser = argparse.ArgumentParser(description = 'Extract a transcription factor module from a gene regulatory network')
parser.add_argument('cyjs', type = str, help = 'Cytoscape JSON file of the GRN')
parser.add_argument('motif', type = str, help = 'Name of the TF motif to extract the module for')
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

# Extract targets from GRN
module = []

for source_node, target_node, metadata in grn.edges(data = True):
	if metadata['source_motif'] == args.motif:
		module.append(target_node)

# Count number of targets
# If none are found, return an error and exit
if len(module) == 0:
	sys.stderr.write('Warning: No genes found in TF module. Check spelling of motif ID and try again. Exiting\n')
	sys.exit(1)
else:
	sys.stdout.write('Found %d genes\n' % len(module))

###################################################################################################################################################################################

# Write module genes to file
out = open(args.out, 'w')

for gene_id in sorted(module):
	out.write(gene_id + '\n')

out.close()

###################################################################################################################################################################################