#!/bin/bash

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
#SBATCH --account=pawsey0809
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64gb
#SBATCH --output=./process_frb_search.o%j
#SBATCH --error=./process_frb_search.e%j
#SBATCH --export=NONE

echo "module use  /software/projects/director2183/msok/setonix/modules/"
module use  /software/projects/director2183/msok/setonix/modules/

echo "module load msfitslib/devel"
module load msfitslib/devel

echo "Checks:"
echo "which dynaspec_search"
which dynaspec_search
dynaspec_search

dataset=20230601_100213_100ms_ch294
if [[ -n "$1" && "$1" != "-" ]]; then
   dataset="$1"
fi
dtm=`echo $dataset | cut -b 1-15`
ux=`date2date -ut2ux=${dtm} | awk '{print $3;}'`
gps=$(($ux-315964783))

create_dynaspec=1
if [[ -n "$2" && "$2" != "-" ]]; then
   create_dynaspec=$2
fi

freq_ch=294
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_ch=$3
fi

# 229.23900462962962962963 for ch=294
# center of the first fine channel:
freq_mhz=`echo $freq_ch | awk '{freq_mhz=$1*(400.00/512.00)-(400.00/512.00)*(32.00/27.00)*0.5 + (((400.00/512.00)*(32.00/27.00))/32.00)*0.5;printf("%.20f\n",freq_mhz);}'`


echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "dataset = $dataset ( dtm = $dtm , ux = $ux , gps = $gps )"
echo "create_dynaspec = $create_dynaspec"
echo "freq_ch = $freq_ch -> freq = $freq_mhz [MHz]"
echo "##########################################"


if [[ ! -d ${dataset} ]]; then
   echo "ERROR : $dataset not found in location:"
   pwd
   exit -1
fi

cd ${dataset}
pwd

if [[ -d fits_images ]]; then
   cd fits_images
   
   # create fits_list files for every frequency channel imaged:
   for ch in `ls -d ?????`; 
   do    
      cd $ch;    
      ls dirty_image*_real.fits > fits_list;    
      cd ..; 
   done
   
   # -X - time resolution in seconds
   # -A - frequency resolution in MHz 
   # -I - telescope type 2 is EDA2/AAVS2 -> calculation of start frequency based on this :
   # -C ${freq_ch} : channel ${freq_ch} of EDA2/AAVS2 :
   if [[ $create_dynaspec -gt 0 ]]; then
      echo "create_dynaspec -w \"(50,50)-(150,150)\" -o ${gps} -f dirty_image_%dT%d_real.fits -v 10 -F -N 32 -X 0.100 -I 2 -C ${freq_ch} "
      create_dynaspec -w "(50,50)-(150,150)" -o ${gps} -f dirty_image_%dT%d_real.fits -v 10 -F -N 32 -X 0.100 -I 2 -C ${freq_ch} 
   else
      echo "WARNING : create_dynaspec it not required"
   fi
      
   cd dynamic_spectra/   
   # 229.239004629 
   # 229.23900462962962962963
   # -B is start frequency = ch*(400/512)-(400/512)*(32/27)/2 - (400/512)*(32/27)/(32)/2 - center of the first fine channel
   # for example for ch=294 , freq = 294*(400.00/512.00)-(400.00/512.00)*(32.00/27.00)/2 +  (400/512)*(32/27)/(32)/2 = 229.23900462962962962963
   #                                                                                                                   229.23900462962962962963
   # was -B 229.23900462962962962963 -E 294
   echo "timeseries_all.sh \"????_????\" 0 0 5 \"-S 1 -A -T -C 0.02893518518518518518 -I 2 -B ${freq_mhz} -E ${freq_ch} -P 0\" 1 ${gps} 0 900 150 > all.out 2>&1"
   timeseries_all.sh "????_????" 0 0 5 "-S 1 -A -T -C 0.02893518518518518518 -I 2 -B ${freq_mhz} -E ${freq_ch} -P 0" 1 ${gps} 0 900 150 > all.out 2>&1 
   cd ..
else
   echo "ERROR : fits_images directory does not exist"
fi
