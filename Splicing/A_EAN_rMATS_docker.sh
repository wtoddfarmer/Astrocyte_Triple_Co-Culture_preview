#! /bin/bash

docker run -v /Volumes/two/DATA/Astro_Endo_Neuron_RNAseq/:/EAN -w /EAN  xinglab/rmats \
--task both \
--b1 splicing/A.txt \
--b2 splicing/EAN.txt \
--gtf mm10_withEGFP.gtf \
-t paired \
--libType fr-firststrand \
--readLength 101 \
--nthread 20 \
--od A_EAN_out \
--tmp A_EAN_temp 
