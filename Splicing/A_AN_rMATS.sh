#! /bin/bash

docker run -v /Volumes/two/DATA/Astro_Endo_Neuron_RNAseq/:/EAN -w /EAN  xinglab/rmats \
--task both \
--b1 splicing/A.txt \
--b2 splicing/AN.txt \
--gtf mm10_withEGFP.gtf \
-t paired \
--libType fr-firststrand \
--readLength 101 \
--nthread 10 \
--od A_AN_out \
--tmp A_AN_temp 
