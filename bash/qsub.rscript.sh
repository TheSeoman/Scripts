#!/bin/sh
script_name=$1
mem_gb=$2
num_cores=$3
from=$4
to=$5
r_script_dir="/home/icb/julian.schmidt/Scripts/R/"
qsub_out_dir="/home/icb/julian.schmidt/qsub/"
qsub_script_dir="/home/icb/julian.schmidt/Scripts/qsub/"
qsub_calls_file="/home/icb/julian.schmidt/qsub.calls.txt"

qsub_script_file=$qsub_script_dir'run.'$script_name'.sh'
r_script_file=$r_script_dir$script_name'.R'

if [ ! -f $qsub_script_file ]; then
	printf 'Generating '$qsub_script_file'\n'
	printf '#!/bin/sh\nRscript '$r_script_file' $1\n' > $qsub_script_file
fi

qsub_out_file=$qsub_out_dir'out.'$script_name'-1.txt'
qsub_err_file=$qsub_out_dir'err.'$script_name'-1.txt'
i=1
while [ -f $qsub_out_file ]; do
	i=$[i+1]
	qsub_out_file=$qsub_out_dir'out.'$script_name'-'$i'.txt'
	qsub_err_file=$qsub_out_dir'err.'$script_name'-'$i'.txt'
done


cores=""
if [ $num_cores > 1 ]; then
	cores="-pe smp $num_cores"
fi



if [ ! -z $from ] && [ ! -z $to  ]; then
	qsub_out_file=$qsub_out_dir'out.'$script_name'-1:1.txt'
	qsub_err_file=$qsub_out_dir'err.'$script_name'-1:1.txt'
	i=1
	while [ -f $qsub_out_file ]; do
		i=$[i+1]
		qsub_out_file=$qsub_out_dir'out.'$script_name'-'$i'#1.txt'
		qsub_err_file=$qsub_out_dir'err.'$script_name'-'$i'#1.txt'
	done
	for j in `seq $from $to`; do
		qsub_out_file=$qsub_out_dir'out.'$script_name'-'$i'#'$j'.txt'
		qsub_err_file=$qsub_out_dir'err.'$script_name'-'$i'#'$j'.txt'
		qsub -q long_fed25 -hard -l job_mem="$mem_gb"G $cores -V -o $qsub_out_file -e $qsub_err_file $qsub_script_file $j
	done
else	
	qsub_out_file=$qsub_out_dir'out.'$script_name'-1.txt'
	qsub_err_file=$qsub_out_dir'err.'$script_name'-1.txt'
	i=1
	while [ -f $qsub_out_file ]; do
		i=$[i+1]
		qsub_out_file=$qsub_out_dir'out.'$script_name'-'$i'.txt'
		qsub_err_file=$qsub_out_dir'err.'$script_name'-'$i'.txt'
	done
	qsub -q long_fed25 -hard -l job_mem="$mem_gb"G $cores -V -o $qsub_out_file -e $qsub_err_file $qsub_script_file 
fi

DATE=`date '+%Y-%m-%d %H:%M:%S'`
printf "$DATE qsub -q long_fed25 -hard -l job_mem="$mem_gb"G $cores -V -o $qsub_out_file -e $qsub_err_file $qsub_script_file\n" >> $qsub_calls_file