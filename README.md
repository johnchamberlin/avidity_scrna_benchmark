# Avidity paired-end scRNA-seq benchmark

### A comparison of 150bp x 150bp 3' scRNA-seq sequencing configuration using Illumina and Element sequencers 

This repository contains analysis code used by Chamberlin et al 2024: https://doi.org/10.1101/2024.07.10.602909

This analysis entails a comparison of five libraries sequenced with 150bp x 150bp configuration on
Illumina and Element short-read sequencers. Single and paired-end alignment are performed with STARsolo,
and polyadenylation site analysis was conducted with scraps (https://github.com/rnabioco/scraps) and polyApipe (https://github.com/MonashBioinformaticsPlatform/polyApipe). 

Site-level polyadenylation site quantifications from scraps and polyApipe are hosted in this repository (data/).

Raw sequencing data are available on GEO at GSE143363 (Library 1) and GSE232559 (Libraries 2,3,4,5)

Our STAR and cutadapt parameters are available in the 'scripts' folder
