#!/bin/sh
for i in {1..27}; do
	qsub -q long_fed25 -hard -l job_mem=4G -V -o /home/icb/julian.schmidt/qsub/out.single.chromHMM-1-$i.txt -e /home/icb/julian.schmidt/qsub/err.single.chromHMM-1-$i.txt /home/icb/julian.schmidt/Scripts/qsub/run.single.chromHMM.sh $i
done