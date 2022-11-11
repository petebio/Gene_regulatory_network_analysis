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
<p><b>build_gene_regulatiry_network.py</b> - The main script for building a GRN. This script requires a BED file of regulatory elements, such as DNaseI or ATAC-Seq sites that are to be included in the GRN</p>

Required file formats
----------

### Gene expression data
A tab-delimited file with two columns: 1st - gene ID, 2nd - Gene expression value (e.g. FPKM)


