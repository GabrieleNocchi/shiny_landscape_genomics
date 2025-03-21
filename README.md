My first attempt at building a shiny app :)

First it does a quick structure assessment with snmf to help choosing number of latent factors for LFMM.

Then it downloads the 19 Bioclim variables for your dataset (0.5 res), it calculates correlation between variables, it plots the correlation and removes highly correlated variables using a user defined correlation threshold.

Then it runs LFMM using the number of latent factors (K) chosen by the user.

Finally it allows to save final LFMM manhattan plot and as well as a table with significant SNPs results at the chosen FDR cutoff.

Unzip plink executable before running! Needs Plink PED/MAP file for genotypes and a geographic coordinates tab separated file with 3 columns and the following header: ID LON LAT

Still work in progress ...

Have fun!
