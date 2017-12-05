#!/bin/sh
for i in {1..7}; do
	qsub -q long_fed25 -hard -l job_mem=4G -V -o /home/icb/julian.schmidt/qsub/out.herv.split.snp.overlap-2-$i.txt -e /home/icb/julian.schmidt/qsub/err.herv.split.snp.overlap-2-$i.txt /home/icb/julian.schmidt/Scripts/qsub/run.herv.split.snp.overlap.sh $(( (i-1)*100 )) $(( i*100-1 ))
done