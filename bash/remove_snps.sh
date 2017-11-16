#!/bin/sh
cd /data/media/Masterarbeit/data/SNPs
indices=$(<snp.indices.txt)
sed -i.bak -e $indices withexpr_sorted > S2_sorted