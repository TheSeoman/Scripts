#!/bin/sh
cd /data/media/Masterarbeit/data/SNPs
indices=$(<snp_indices.txt)
sed -i.bak -e $indices withexpr_sorted > S2_sorted