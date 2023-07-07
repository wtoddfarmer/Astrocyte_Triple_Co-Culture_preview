#! /bin/bash

docker run -v /Volumes/two/DATA/Astro_Endo_Neuron_RNAseq/:/EAN -w /EAN  xinglab/rmats \
--task both \
--b1 splicing/A.txt \
--b2 splicing/EA.txt \
--gtf mm10_withEGFP.gtf \
-t paired \
--libType fr-firststrand \
--readLength 101 \
--nthread 8 \
--od A_AE_out \
--tmp A_AE_temp 
