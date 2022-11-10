#!/usr/bin/env python
import pybedtools
import argparse
import numpy
import os

###############################################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('bed', type = str, help = 'BED file to use for motif search')
parser.add_argument('dir', type = str, help = 'Directory of motif files to use with Homer')
parser.add_argument('genome', type = str, help = 'Genome version to use with Homer')
parser.add_argument('out', type = str, help = 'Output directory')
parser.add_argument('-d', '--dist', dest = 'd', type = int, default = 2, help = 'Maximum distance between duplicate motif pairs. Default = 2')

args = parser.parse_args()

###############################################################################################################################################################

# Convert relative path of bed file and input directory to an absolute path
bed = os.path.abspath(args.bed)
motif_dir = os.path.abspath(args.dir)

# If the output directory does not already exists, create it, then go there

if not os.path.exists(args.out):
	os.mkdir(args.out)

os.chdir(args.out)

###############################################################################################################################################################

for motif_file in os.listdir(motif_dir):
	if not motif_file.endswith('.motif'):
		continue

	motif = motif_file.replace('.motif', '')

	# Check if BED file already exists
	if os.path.exists('%s.bed' % motif):
		continue

	# Search for motif in BED file using Homer
	homer = 'annotatePeaks.pl %s %s -noann -m %s/%s -mbed %s_raw.bed > Homer_output.tsv' % (bed, args.genome, motif_dir, motif_file, motif)
	os.system(homer)

	# Homer has trouble dealing with palindromic motifs which result in each motif being found twice.
	# These need to be removed.

	motif_bed = pybedtools.BedTool('%s_raw.bed' % motif).sort()
	prev_line = list()

	out = open('%s.bed' % motif, 'w')

	for site in motif_bed:
		chrom = site[0]
		start = int(site[1])
		end = int(site[2])
		strand = site[5]

		# Get the center position of the motif
		pos = int(numpy.ceil(numpy.mean([start, end])))

		is_duplicate = False

		# A motif is considered a duplicate if:
		# The motif occurs on the same chromosome but opposite strand than the motif found on the previous line, and
		# the distance between the center positions of these motifs is < 2

		if len(prev_line) != 0:
			prev_chrom, prev_pos, prev_strand = prev_line

			if chrom == prev_chrom and strand != prev_strand:
				dist = abs(pos - prev_pos)

				if dist <= args.d:
					is_duplicate = True

		if not is_duplicate:
			out.write('%s\t%d\t%d\n' % (chrom, start, end))

		prev_line = [chrom, pos, strand]

	# Clean up
	os.system('rm %s_raw.bed' % motif)
	os.system('rm Homer_output.tsv')

###############################################################################################################################################################
