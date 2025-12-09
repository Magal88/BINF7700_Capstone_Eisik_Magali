#!/bin/bash

# Script to download the normalized beta values matrix from GEO

wget -O GSE55763_normalized_betas.txt.gz \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55763&format=file&file=GSE55763%5Fnormalized%5Fbetas%2Etxt%2Egz"

