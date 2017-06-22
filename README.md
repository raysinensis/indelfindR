# indelfindR
finding possible combinations of indel alleles from inhomogenous sanger sequencing results
# can now be ran as shiny app (indelFindR.shiny folder)

Genotyping of CRISPR-cas9-edited embryonic cells can be done through PCR and DNA sequencing. DNA sequencing data may turn out inhomogeneous, if the alleles arise late, or cells from multple zygotes were mixed together.
This is a short R script to anaylze the sequencing ab1 file, call multiple peaks at each position, and screen the possible indel combinations that may result in those results.

Currently the script searches with 2 methods: 1) indelfindR-search script, looks at whether indels of up to 20nt can fit into a 50nt section of those chromatographs. This method is quick but only lists all possible candidates. 2) Alternatively, the dnax4or6 scripts detect up to 6 different alleles (a much quicker 4x version is ran first, around 1-3 min. if no suitable results are found, the 6x version will take 1-2 hrs), outputing whether each allele may contain insertions or deletions less than 10bp. Results are only listed when the combinations match sequencing results exactly. Point mutations, and the specific insertions, can then be determined by eye. (coming later)

other variants can be found in the variants folder:
dnax4or6.r is the original script where the user has to define the signal threshold (in my own experience 0.1 works well)
dnax4or6autocall.r considers the ~40 nt signals around the area of interest, and picks the largest drop-off for each of A,C,G, and T as the threshold, this is the recommended script to start with
dnax4or6autocall-inputpath.r can be executed on mac or linux command line directly, where the user will be prompted to enter input path and filename (strangely some mac users with foreign languages have problems editing the previous scripts to their custom path and file names, this is the workaround)

As an alternative approach for simpler cases, BackgroundRead.r reads out different allelles as seperate files, where when no secondary base read is over the automatically determined threshold, the primary read is given instead. Suitable for quick blast afterwards
