#!/bin/bash
#SBATCH --job-name=nobecostrict
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

module load vim/7.4 
module load java/1.8.0_121
module load emacs/25.1 
module load cmake/3.15.1
module load python/3.6 
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

export CC=`which gcc`
export CXX=`which c++`

export LD_LIBRARY_PATH=/global/home/users/cbrook/beagle-lib/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/global/home/users/cbrook/beagle-lib/lib/pkgconfig:$PKG_CONFIG_PATH
export BEAST_EXTRA_LIBS=/global/home/users/cbrook/beagle-lib/lib/:$BEAST_EXTRA_LIBS

/global/home/users/cbrook/beast/bin/beast -beagle_GPU -beagle_CPU -beagle_SSE -threads 12  /global/home/users/cbrook/RdRpBEASTxml/NobecoStrict/AllNobecoStrict.xml
