# indelfindR
finding possible combinations of indel alleles from inhomogenous sanger sequencing results

Genotyping of CRISPR-cas9-edited embryonic cells can be done through PCR and DNA sequencing. DNA sequencing data may turn out inhomogeneous, if the alleles arise late, or cells from multple zygotes were mixed together.
This is a short R script to anaylze the sequencing ab1 file, call multiple peaks at each position, and screen the possible indel combinations that may result in those results.

Currently the script detects up to 6 different alleles, outputing whether each allele may contain insertions or deletions less than 10bp. Point mutations, and the specific insertions, can then be determined by eye.
