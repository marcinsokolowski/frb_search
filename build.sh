#!/bin/bash -e

build_dir=build
cmake_options=""

build_opt="Release" # or  "Debug"
if [[ -n "$1" && "$1" != "-" ]]; then
   build_opt="$1"
fi


# default options decided based on build type can be over-written by the 2nd parameter:
if [[ -n "$2" && "$2" != "-" ]]; then
   cmake_options=$2
fi

version=""
if [[ -n "$3" && "$3" != "-" ]]; then
   version="$3"
   build_dir=${build_dir}_${version}
fi

# First, you need to source the bash library
module load bash-utils
echo "source ${BASH_UTILS_DIR}/build_utils.sh"
source "${BASH_UTILS_DIR}/build_utils.sh"


PROGRAM_NAME=frb-seaarch
PROGRAM_VERSION=main

if [[ -n "$4" && "$4" != "-" ]]; then
   PROGRAM_VERSION="$4"
fi


echo "############################################"
echo "PARAMETERS (build.sh scripts) :"
echo "############################################"
echo "PROGRAM_NAME = $PROGRAM_NAME"
echo "PROGRAM_VERSION = $PROGRAM_VERSION"
echo "############################################"


 
# the following function sets up the installation path according to the
# cluster the script is running on and the first argument given. The argument
# can be:
# - "group": install the software in the group wide directory
# - "user": install the software only for the current user
# - "test": install the software in the current working directory 
process_build_script_input user


# load all the modules required for the program to compile and run.
# the following command also adds those module names in the modulefile
# that this script will generate.
echo "Loading required modules ..."
echo "Loading modules for PAWSEY_CLUSTER = $PAWSEY_CLUSTER"
module reset
module_load cfitsio/4.4.0  msfitslib/master-ittkjmq  fftw/3.3.10  pal/0.9.8-n3thcaw  libnova/0.15.0-iwh6cpn 
   
# cmake is only required at build time, so we use the normal module load
module load cmake/3.30.5
# build your software..
echo "Building the software.."

[ -d ${build_dir} ] || mkdir ${build_dir}
cd ${build_dir}
cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-O3 ${cmake_options}
make -j 12 VERBOSE=1
# make test
# Install the software
# make install

# test:

echo "Create the modulefile in $MODULEFILE_DIR (or $INSTALL_DIR)"
export ADDITIONAL_MODULEFILE_COMMANDS="prepend_path('BLINK_IMAGER_PATH', root_dir )"
create_modulefile

echo "Done."


