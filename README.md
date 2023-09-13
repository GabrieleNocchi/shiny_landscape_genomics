# shiny_landscape_genomics
My first attempt at building a shiny app :)


Does a quick structure assessment with snmf to help choosing number of latent factors.

Then it download the 19 Bioclim variables for your dataset (0.5 res) and  run LFMM, allowing the user to save final manhattan plot and a table with significant SNPs results at the chosen FDR cutoff.

Unzip plink executable before running!
Needs Plink PED/MAP file for genotypes and a geographic coordinates tab separated file with 3 columns and the following header:  ID  LON  LAT
