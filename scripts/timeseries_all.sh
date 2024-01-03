#!/bin/bash

# ls 05??_????.fits  | awk '{print substr($1,1,4);}' | sort -u
prefix="????_????"
if [[ -n "$1" && "$1" != "-" ]]; then
   prefix=$1
fi

rm_series_fits=0
if [[ -n "$2" && "$2" != "-" ]]; then
   rm_series_fits=$2
fi

force=0
if [[ -n "$3" && "$3" != "-" ]]; then
   force=$3
fi

snr_threshold=3 # was 5 , but can be filtered out later (if all saved >3 sigma)
if [[ -n "$4" && "$4" != "-" ]]; then
   snr_threshold=$4
fi

options=""
if [[ -n "$5" && "$5" != "-" ]]; then
   options="$5"
fi

outdir="./"
use_subdir=1
if [[ -n "$6" && "$6" != "-" ]]; then
   use_subdir=$6
   mkdir -p ${outdir}
fi

obsid=1217495184
if [[ -n "$7" && "$7" != "-" ]]; then
   obsid=$7
fi

dm_min=100
if [[ -n "$8" && "$8" != "-" ]]; then
   dm_min=$8
fi

dm_max=3000
if [[ -n "$9" && "$9" != "-" ]]; then
   dm_max=$9
fi

dm_step=1 # should be calculated with a script : python ./dispersion_dm_step.py 100.00 230.00 .925926
if [[ -n "${10}" && "${10}" != "-" ]]; then
   dm_step=${10}
fi

datadir=`pwd`

echo "################################################"
echo "PARAMETERS:"
echo "################################################"
echo "obsid      = $obsid"
echo "DM range   = $dm_min - $dm_max"
echo "DM step    = $dm_step"
echo "use_subdir = $use_subdir"
echo "datadir    = $datadir"
echo "################################################"


for fits in `ls ${prefix}.fits`
do
  outdir=${fits%%.fits}
  echo "Processing FITS file $fits -> outdir = ${outdir}"
  
  if [[ $use_subdir -gt 0 ]]; then
     mkdir -p ${outdir}
     cd ${outdir}
     ln -s ${datadir}/${fits}
  fi

  # 5 sigma :
  out_fits=${fits%%.fits}_series.fits
  out_stat=${fits%%.fits}_series.stat

  if [[ ! -s ${out_fits} || $force -gt 0 ]]; then   # ! -s ${out_fits}
     echo "which dynaspec_search"
     which dynaspec_search
     echo "dynaspec_search ${fits} ${out_fits}  -o ${obsid} -l $dm_min -m $dm_max -s ${dm_step} -n ${snr_threshold} ${options}"
     dynaspec_search ${fits} ${out_fits}  -o ${obsid} -l $dm_min -m $dm_max -s ${dm_step} -n ${snr_threshold} ${options}
  
     echo "calcstat ${out_fits} > ${out_stat}"
     calcstat ${out_fits} > ${out_stat}
  else   
     echo "WARNING : file ${out_stat} exists -> not re-creating ! use option force to force re-creation !"
  fi
  
  if [[ $rm_series_fits -gt 0 ]]; then
     echo "rm -f *_series.fits *_count.fits"
     rm -f *_series.fits *_count.fits
  else
     echo "WARNING : removal of files *_series.fits *_count.fits is not required -> they may use a lot of space ..."
  fi
  
  if [[ use_subdir -gt 0 ]]; then
     cd -
  fi
done
