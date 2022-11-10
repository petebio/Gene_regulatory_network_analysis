#!/usr/bin/env python
from numpy import random
import pybedtools as pb
import argparse
import sys
import os

#############################################################################################################################################################################################

parser = argparse.ArgumentParser(description = 'Annotate a BED file to an associated gene using CHiC data')
parser.add_argument('bed', type = str, help = 'BED file to annotate')
parser.add_argument('chic', type = str, help = 'BED file of CHiC interactions. Target gene should be 6th column')
parser.add_argument('out', type = str, help = 'Output file')
parser.add_argument('-g', '--genome', dest = 'g', type = str, default = 'hg38', help = 'Genome version to use with Homer (for sites with no CHiC annotation). Default = hg38')
parser.add_argument('-d', '--dist', dest = 'd', type = int, default = 200000, help = 'Maximum distance between peak and target gene (for sites with no CHiC annotation). Default = 200000')

args = parser.parse_args()

#############################################################################################################################################################################################

# Read CHiC annotations
chic_annotated_bed = list()

with open(args.chic, 'r') as chic:
	for num,line in enumerate(chic):
		if num == 0:
			continue

		fields = line.strip().split('\t')

		chrom = fields[0]
		start = fields[1]
		end = fields[2]
		geneID = fields[4]

		site = '\t'.join([chrom, start, end, geneID])
		chic_annotated_bed.append(site)

chic_annotated_bed = pb.BedTool(chic_annotated_bed).sort()
sys.stderr.write('Read %d sites annotated by CHiC\n' % chic_annotated_bed.count())

#############################################################################################################################################################################################

# Read bed file and separate sites based on:
# - Sites with a CHiC annotation
# - Sites with no ChIC annotation

bed = pb.BedTool(args.bed).sort()
bed_with_CHiC = bed.intersect(chic_annotated_bed, wa = True, u = True)
bed_without_CHiC = bed.intersect(chic_annotated_bed, v = True)

perc_with_CHiC = bed_with_CHiC.count() / bed.count() * 100
perc_without_CHiC = bed_without_CHiC.count() / bed.count() * 100

sys.stderr.write('Read %d sites from %s\n' % (bed.count(), args.bed))
sys.stderr.write('\t- %d (%.2f%%) of these could be annoted by CHiC\n' % (bed_with_CHiC.count(), perc_with_CHiC))
sys.stderr.write('\t- %d (%.2f%%) of these could not be annoted by CHiC\n' % (bed_without_CHiC.count(), perc_without_CHiC))

#############################################################################################################################################################################################

# Sites annotated by CHiC
# Extract CHiC annotated gene from CHiC annotation

bed_with_CHiC = bed_with_CHiC.intersect(chic_annotated_bed, wo = True)
annotated_peaks = dict()

for site in bed_with_CHiC:
	chrom = site[0]
	start = site[1]
	end = site[2]
	geneID = site[6]

	peak = '%s\t%s\t%s' % (chrom, start, end)

	annotated_peaks[peak] = geneID

#############################################################################################################################################################################################

# Sites that cannot be annotated by CHiC
# Write sites as a bed file and use Homer's annotatePeaks.pl to find closest gene

sys.stderr.write('\nAnnotating remaining peaks with Homer\n')

bed_tmp_file = 'tmp.%d' % int(random.randint(1e6))
homer_tmp_output = 'tmp.%d' % int(random.randint(1e6))

bed_without_CHiC.saveas(bed_tmp_file)

homer = 'annotatePeaks.pl %s %s -noann > %s' % (bed_tmp_file, args.g, homer_tmp_output)
os.system(homer)

annotated_count = 0

with open(homer_tmp_output, 'r') as homer:
	for num,line in enumerate(homer):
		if num == 0:
			continue

		fields = line.strip().split('\t')

		chrom = fields[1]
		start = int(fields[2])
		end = int(fields[3])

		try:
			dist_to_tss = abs(int(fields[9]))
		except:
			dist_to_tss = args.d + 1

		if dist_to_tss <= args.d:

			if len(fields) == 19:
				geneID = fields[15]
			else:
				continue

			peak = '%s\t%d\t%d' % (chrom, start - 1, end) # Homer automatically adds 1 to every start position - remove this
			annotated_peaks[peak] = geneID
			annotated_count = annotated_count + 1

sys.stderr.write('Successfully annotated %d sites to closest gene within %d bp\n' % (annotated_count, args.d))

os.remove(bed_tmp_file)
os.remove(homer_tmp_output)

#############################################################################################################################################################################################

# Convert annotated peak table back to a BED file and write to output

annotated_bed = []

for site in annotated_peaks:
	annotated_bed.append('%s\t%s' % (site, annotated_peaks[site]))

annotated_bed = pb.BedTool(annotated_bed).sort()
annotated_bed.saveas(args.out)

sys.stderr.write('Wrote %d annotated peaks to %s\n' % (annotated_bed.count(), args.out))

#############################################################################################################################################################################################
