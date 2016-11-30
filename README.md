# indelfindR
finding possible combinations of indel alleles from inhomogenous sanger sequencing results

Genotyping of CRISPR-cas9-edited embryonic cells can be done through PCR and DNA sequencing. DNA sequencing data may turn out inhomogeneous, if the alleles arise late, or cells from multple zygotes were mixed together.
This is a short R script to anaylze the sequencing ab1 file, call multiple peaks at each position, and screen the possible indel combinations that may result in those results.

Currently the script detects up to 6 different alleles (a much quicker 4x version is ran first, around 1-3 min. if no suitable results are found, the 6x version will take 1-2 hrs), outputing whether each allele may contain insertions or deletions less than 10bp. Point mutations, and the specific insertions, can then be determined by eye.

Different verisions of the script:
dnax4or6.r is the original script where the user has to define the signal threshold (in my own experience 0.1 works well)
dnax4or6autocall.r considers the ~40 nt signals around the area of interest, and picks the largest drop-off for each of A,C,G, and T as the threshold, this is the recommended script to start with
dnax4or6autocall-inputpath.r can be executed on mac or linux command line directly, where the user will be prompted to enter input path and filename (strangely some mac users with foreign languages have problems editing the previous scripts to their custom path and file names, this is the workaround)
