#!/bin/bash
echo " "
#export INTEL_LICENSE_FILE=28518@skynet.chem.ttu.edu
#source  /opt/intel/Compiler/11.1/064/bin/ifortvars.sh intel64
#source  /opt/intel/Compiler/11.1/064/bin/iccvars.sh intel64
#export PATH=/lustre/work/apps/mvapich-test/bin:$PATH
#export LD_LIBRARY_PATH=/lustre/work/apps/mvapich-test/lib:/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
export mpidir=/cm/shared/apps/openmpi/intel/64/1.3.3
cd /home/carolp/work/CPC-D-12-00495-REV/NWC/nwchem-6.1.1-src
export NWCHEM_TOP=`pwd`
length=`pwd | wc -m`
if [ $length -gt 64 ]
then
  echo "The directory name chosen for NWCHEM_TOP is longer"
  echo "than the maximum allowed value of 65 characters"
  echo "Please shorten the path"
  echo "Exiting"
  exit 1
fi
cd src
export TCGRSH=/usr/bin/ssh
export NWCHEM_TARGET=LINUX64
export USE_MPI=y
export USE_MPIF=y
export USE_MPIF4=y
export CFLAGS="-pc80 -axSSE4.2 -axSSE4.2 -O3 -ip -g"
export FFLAGS="-pc80 -axSSE4.2 -axSSE4.2 -O3 -ip -g"
export F90FLAGS="-pc80 -axSSE4.2 -axSSE4.2 O3 ip -g"
export MPI_LOC=$mpidir
export MPI_LIB=$MPI_LOC/lib64
export LIBMPI="-L$MPI_LIB -L/usr/lib64 -lmpi -lmpi_f77 -lmpi_f90"
export MPI_INCLUDE=$MPI_LOC/include
export ARMCI_DEFAULT_SHMMAX=256
export MA_USE_ARMCI_MEM=1
export LARGE_FILES="true"
export NWCHEM_MODULES="all"
export F77=ifort
export F90=ifort
export FC=ifort
export FL=ifort
export CC=icc
export CXX=icpc
make realclean 
find . -name include_stamp -exec rm {} \; -print
#make clean
make nwchem_config FC=ifort CC=icc > make.log 
make CC=icc FC=ifort CXX=icpc >> make.log
echo "Done Making Nwchem"
export NWCHEM_MODULES="qm geninterface"
make nwchem_config NWCHEM_MODULES=venus
make nwchem_config
make CC=icc FC=ifort >> make.log
make stubs.o >> make.log
cd ../../
cd venus-nwchem
make clean > make-venus.log
make FC=ifort CC=icc CXX=icpc >> make-venus.log




