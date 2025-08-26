#!/bin/bash

# INFO : use - to skip parameters (use default values)

# on SETONIX also :
# salloc --mem 64g --time 02:00:00 --nodes=1 # request interactive session on CPU node 
# module use /software/projects/pawsey1154/msok/setonix/2025.08/modules/zen3/gcc/14.2.0
# module load frb-search/main # load relevant module 

pixel_x=377
if [[ -n "$1" && "$1" != "-" ]]; then
   pixel_x=$1
fi

pixel_y=896
if [[ -n "$2" && "$2" != "-" ]]; then
   pixel_y=$2
fi

n_times=500
if [[ -n "$3" && "$3" != "-" ]]; then
   n_times=$3
fi

first_coarse_channel=109
if [[ -n "$4" && "$4" != "-" ]]; then
   first_coarse_channel=$4
fi

n_fine_ch=10
if [[ -n "$5" && "$5" != "-" ]]; then
   n_fine_ch=$5
fi

timeres=0.02 # time resolution on seconds :
if [[ -n "$6" && "$6" != "-" ]]; then
   timeres=$6
fi

outdir=dynamic_spectrum
if [[ -n "$7" && "$7" != "-" ]]; then
   outdir="$7"
fi

obsid=1192477696
if [[ -n "$8" && "$8" != "-" ]]; then
   obsid=$8
fi

start_second=1508442485
if [[ -n "$9" && "$9" != "-" ]]; then
   start_second=$9
fi

echo "create_dynaspec -p \"(${pixel_x},${pixel_y})\" -o ${obsid} -S ${start_second} -f start_time_%d_int_%02d_coarse_%03d_fine_ch%02d_image_real.fits -v 10 -N ${n_fine_ch} -X ${timeres} -I 1 -C ${first_coarse_channel} -T ./ -t ${n_times} -d ${outdir} -P "
create_dynaspec -p "(${pixel_x},${pixel_y})" -o ${obsid} -S ${start_second} -f start_time_%d_int_%02d_coarse_%03d_fine_ch%02d_image_real.fits -v 10 -N ${n_fine_ch} -X ${timeres} -I 1 -C ${first_coarse_channel} -T ./ -t ${n_times} -d ${outdir} -P 

