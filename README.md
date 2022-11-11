# Gene regulatory network analysis

Software for the creation and analysis of Gene Regulatory Networks (GRNs)

Requirements
----------
Python packages:
* Numpy
* Pybedtools
* Networkx
* Json

Software:
* Homer
* Cytoscape

Basic Usage
----------
<p><b>build_gene_regulatory_network.py</b> - The main script for building a GRN
##### Required files
1. A BED file of regulatory elements, such as DNaseI or ATAC-Seq sites that are to be included in the GRN. Sites must be annotated to their associated gene, which should be included as the 4th column in the BED file.
2. A directory of BED files with the genomic coordinates for each TF motif to be included in the GRN. This can be created using the <b>findMotifs.py</b> script included with this software
3. Gene expression data file (tab-delimited) with gene ID (1st column) and gene expression value (e.g. FPKM; 2nd column)
4. Transcription Factor (TF) annotation file with gene ID and name of the motif it can bind to. An annotation file (TF_family_gene_annotation.tsv) is provided with this software.

##### Optional Arguments
<b>-a</b>By default this will build a GRN which only includes transcription factor genes. If you would like to build the full network, which includes all genes, use the <b>-a</b> option

</p>

Required file formats
----------

### Gene expression data
A tab-delimited file with two columns: 1st - gene ID, 2nd - Gene expression value (e.g. FPKM)


