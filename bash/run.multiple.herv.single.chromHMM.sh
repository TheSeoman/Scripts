#!/bin/sh
for i in {1..27}; do
	qsub -q long_fed25 -hard -l job_mem=4G -V -o /home/icb/julian.schmidt/qsub/out.herv.single.chromHMM-2-$i.txt -e /home/icb/julian.schmidt/qsub/err.herv.single.chromHMM-2-$i.txt /home/icb/julian.schmidt/Scripts/qsub/run.herv.single.chromHMM.sh $i
done