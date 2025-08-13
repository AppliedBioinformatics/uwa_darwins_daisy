# Phylogenomic Insights Into the Diversification and Population Structure of _Scalasia_ Within the Galápagos Islands.
## Abstract


## Materials and Methods.

### Bioinformatic processing of _Scalesia_ samples.
The genome for _Scalesia atractyloides_ was used as the reference dataset for this study, obtained from an 
existing DOI: Dryad: 10.5061/dryad.8gtht76rh <ref>. The FASTA file was indexed using Bowtie2 (v2.4.5)
<ref> and used as the reference to align each of the thirty-four _scalasia_ genomes detailed below.

Illumina paired end reads were obtained from <complete> for thirty-four individuals within the Galápagos Islands, 
collectively representing eighteen different _scalasia_ species <supplementary table with meta-data>. Initial read 
quality metrics were obtained using FastQC (v0.11.9) <ref>. Paired samples were then processed with Trimmomatic using
the parameters:

`ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:99 TOPHRED33`

Post-trimmed quality was reassessed with FastQC and asessed using MultiQC (v1.12) <ref>. Trimmed FASTA files were 
aligned to the reference genome using Bowtie2 (v2.4.5) with the --sensitive flag and insert size parameters -I 0 -X 1000.
sorted BAM files for the thirty-four samples were generated using the Samtools (v0.5.0) with `view` and `sort` and 
`index` commands with default parameters. Initial alignment metrics were generated using Samtools `stats` and visualised
using a custom Python script <maybe more info>.

### Generation of Presence/Abcence Variation Matrix.
The GFF annotation file for the reference dataset was downloaded using the same DOI above and used in tandem with the
SGS Gene Loss (v0.1) <ref> package to call gene loss of each sample compared to the reference dataset for each of the 
thirty-four _scalasia_ BAM files. SGS Gene Loss was run with default parameters to obtain presence/abcence calls across
43093 features from the GFF file.

A custom Python script <more info> was used to parse the SGS Gene Loss output excov files to determine whether each 
sample contained adequate read depth/coverage to prevent the calling of false abcence of features. # Table here? #.
All thirty-four samples were considered suitable for further analysis.

### Functional Annotation of _Scalesia atractyloides_.



