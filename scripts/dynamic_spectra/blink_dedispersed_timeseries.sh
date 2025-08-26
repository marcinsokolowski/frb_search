#!/bin/bash

# INFO : use - to skip parameters (use default values)

# on SETONIX also :
# salloc --mem 64g --time 02:00:00 --nodes=1 # request interactive session on CPU node 
# module use /software/projects/pawsey1154/msok/setonix/2025.08/modules/zen3/gcc/14.2.0
# module load frb-search/main # load relevant module 

dynaspec_fits=0377_0896.fits
if [[ -n "$1" && "$1" != "-" ]]; then
   dynaspec_fits=$1
fi

dts_out_fits=0377_0896_series.fits
if [[ -n "$2" && "$2" != "-" ]]; then
   dts_out_fits=$2
fi

dm_min=50
if [[ -n "$3" && "$3" != "-" ]]; then
   dm_min=$3
fi

dm_max=65 
if [[ -n "$4" && "$4" != "-" ]]; then
   dm_max=$4
fi

dm_step=0.1
if [[ -n "$5" && "$5" != "-" ]]; then
   dm_step=$5
fi

threshold_in_sigmas=7
if [[ -n "$6" && "$6" != "-" ]]; then
   threshold_in_sigmas=$6
fi

first_coarse_channel=109
if [[ -n "$7" && "$7" != "-" ]]; then
   first_coarse_channel=$7
fi

freqres=0.120
if [[ -n "$8" && "$8" != "-" ]]; then
   freqres=$8
fi

obsid=1192477696
if [[ -n "$9" && "$9" != "-" ]]; then
   obsid=$9
fi

start_freq=`echo "$first_coarse_channel $freqres" | awk '{start_freq=$1*1.28-1.28/2.00+$2/2.00;print start_freq;}'`

echo "Start frequency = $start_freq [MHz]"

echo "dynaspec_search ${dynaspec_fits} ${dts_out_fits} -o ${obsid} -l ${dm_min} -m ${dm_max} -s ${dm_step} -n ${threshold_in_sigmas} -S 1 -A -T -C ${freqres} -I 1 -B ${start_freq} -E ${first_coarse_channel} -P 0"
dynaspec_search ${dynaspec_fits} ${dts_out_fits} -o ${obsid} -l ${dm_min} -m ${dm_max} -s ${dm_step} -n ${threshold_in_sigmas} -S 1 -A -T -C ${freqres} -I 1 -B ${start_freq} -E ${first_coarse_channel} -P 0
