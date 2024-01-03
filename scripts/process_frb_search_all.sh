#!/bin/bash

template="2*"
if [[ -n "$1" && "$1" != "-" ]]; then
   template="$1"
fi

datadir=`pwd`
if [[ -n "$2" && "$2" != "-" ]]; then
   datadir="$2"
fi

echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "template = $template"
echo "datadir  = $datadir"
echo "##########################################"


if [[ ! -d ${datadir} ]]; then
   echo "ERROR : $datadir not found in location:"
   pwd
   exit -1
fi

cd ${datadir}
pwd

for dataset in `ls ${template}`
do
   echo "sbatch process_frb_search.sh ${dataset}"
   sbatch process_frb_search.sh ${dataset}
   echo "submitted at :"
   date
   echo
   echo
   sleep 1    
done

