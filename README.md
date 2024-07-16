# Avidity paired-end scRNA-seq benchmark

### A comparison of 150bp x 150bp 3' scRNA-seq sequencing configuration using Illumina and Element sequencers 

This repository contains analysis code used by Chamberlin et al 2024: https://doi.org/10.1101/2024.07.10.602909

This analysis entails a comparison of five libraries sequenced with 150bp x 150bp configuration on
Illumina and Element short-read sequencers. Single and paired-end alignment are performed with STARsolo,
and polyadenylation site analysis was conducted with scraps and polyApipe. 

Site-level polyadenylation site quantifications from scraps and polyApipe are hosted in this repository (data/).

Raw Illumina data are available on GEO (see preprint); Element data are in the process of submission.
Please contact the authors if data are not yet available through GEO, direct access can be provided. 

Scraps installation is described at https://github.com/rnabioco/scraps
PolyApipe is available at https://github.com/MonashBioinformaticsPlatform/polyApipe

Our STAR and cutadapt parameters are available in the 'scripts' folder
