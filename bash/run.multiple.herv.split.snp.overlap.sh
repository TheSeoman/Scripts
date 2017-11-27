#!/bin/sh
for i in {30..50}; do
	qsub -q long_fed25 -hard -l job_mem=4G -V -o /home/icb/julian.schmidt/qsub/out.herv.split.snp.overlap-1-$i.txt -e /home/icb/julian.schmidt/qsub/err.herv.split.snp.overlap-1-$i.txt /home/icb/julian.schmidt/Scripts/qsub/run.herv.split.snp.overlap.sh $i 697
done