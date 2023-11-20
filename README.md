[![DOI](https://sandbox.zenodo.org/badge/564330528.svg)](https://sandbox.zenodo.org/doi/10.5072/zenodo.268)

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
### <b>build_gene_regulatory_network.py</b> - The main script for building a GRN

python build_gene_regulatory_network.py \<BED file of regulatory elements\> \<Motif position directory\> \<Gene expression file\> \<TF annoation file\> \<Output file\>

##### Required files
1. A BED file of regulatory elements, such as DNaseI or ATAC-Seq sites that are to be included in the GRN. Sites must be annotated to their associated gene, which should be included as the 4th column in the BED file.
2. A directory of BED files with the genomic coordinates for each TF motif to be included in the GRN. This can be created using the <b>findMotifs.py</b> script included with this software
3. Gene expression data file (tab-delimited) with gene ID (1st column) and gene expression value (e.g. FPKM; 2nd column)
4. Transcription Factor (TF) annotation file with gene ID and name of the motif it can bind to. An annotation file (TF_family_gene_annotation.tsv) is provided with this software.

##### Optional Arguments
<b>-a</b> By default this will build a GRN which only includes transcription factor genes. If you would like to build the full network, which includes all genes, use the <b>-a</b> option
<br>
<b>-f</b> An optional BED file of DNaseI/ATAC-Seq footprints to filter the motif positions against
<br>
<b>-m</b> Minimium gene expression value for genes to include in the network. Default = 0 (include everything)

##### Output file
The output file is a Cytoscape JSON file (.cyjs) which can be opened and manipulated in Cytoscape

### <b>findMotifs.py</b> - Find the genomic positions for a set of transcription factor binding motifs

python findMotifs.py \<BED file to use for motif search\> \<Directory of motif PWMs to search for\> \<Genome version (e.g. hg38, mm10)\> \<Output directory\>

##### Required files
1. BED file of sites to use for the motif search
2. A directory of Homer compatible motif probabilty weight matrices files. A set of PWMs used in Assi et al (2019) are provided in motif_PWM_database
3. Genome version to use with Homer (e.g. hg38, mm10 etc.)

##### Output
A directory of BED files containing the aligned motif positions


### <b>extract_TF_module_from_GRN.py</b> - Extract a list of nodes connected to a given TF binding motif

python extract_TF_module_from_GRN.py \<GRN in cyjs format\> \<TF motif name\> \<output file\>

##### Required files
1. A gene regulatory network in cyjs format
2. The name of the TF motif to extract the module for
3. Output file

### <b>GRN_to_countMatrix.py</b> - Convert a GRN in cyjs format to a motif count matrix that is more readable to tools such as R and excel

python GRN_to_countMatrix.py \<Cytoscape cyjs file\> \<Output file (.tsv)\>

##### Required files
1. GRN in cytoscape json (cyjs) format
2. Output file

The output of this script is a n x m matrix where columns represent TF motifs and rows gene targets. The entry in each row/column represents the number of TF motifs associated with the gene target

### <b>countMatrix_to_GRN.py</b> - Convert a motif count matrix to a GRN cytoscape json (cyjs) file

python countMatrix_to_GRN.py \<GRN count matrix file\> \<Gene expression file\> \<TF motif annotation file\> \<Output file\>

##### Required files
1. GRN count matrix file
2. Gene expression data file (tab-delimited) with gene ID (1st column) and gene expression value (e.g. FPKM; 2nd column)
3. Transcription Factor (TF) annotation file with gene ID and name of the motif it can bind to. An annotation file (TF_family_gene_annotation.tsv) is provided with this software.

Annotating DNaseI/ATAC sites to their rightful gene with HiC
----------
<p>In order to construct an accurate gene regulatory network we need to be able to associate a cis-regulatory element with the correct gene. Here, we will use promoter-capture HiC from AML patient cells and from healthy CD34+ cells. A set of annotated DNaseI sites are provided for FLT3-ITD, t(8;21), CEBPAx2 and healthy cells are provided in HiC_annoation_data along with a python script <b>annotateBed_with_CHiC.py</b> to help with annotation.</p>

<p>This script will read a bed3 file and attempt to map them to their proper target gene using the annotation file provided. If a peak cannot be annotated with HiC, they will then be assigned to their closest gene using the annotatePeaks.pl function in Homer</p>

### Usage

python annotateBed_with_CHiC.py \<BED file to annotate\> \<CHiC annotation file\> \<Output file\>

##### Required files
1. BED file to annotate
2. HiC annotation file. Files are provided in HiC_annotation_data/data

#### Optional arguments
<b>-g</b> - Genome version to use with Homer. Default = hg38
<br>
<b>-d</b> - Maximum distance for closest gene. Genes with distances greater than this value will not be annotated. Default = 200000
