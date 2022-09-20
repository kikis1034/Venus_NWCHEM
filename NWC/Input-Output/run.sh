#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -o $JOB_NAME.o$JOB_ID
#$ -q all.q
#$ -notify
#$ -pe mpi 12

unset SGE_ROOT
ulimit -l unlimited
export NWCHEM_TOP=/lustre/work/skohale/CPC-D-12-00495a/NWC-V/nwchem-6.1.1
export NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries/
export NWCHEM_TARGET=LINUX64
export USE_MPI=y
export USE_MPIF=y
export MPI_LOC=/lustre/work/apps/mvapich2-1.6/
export MPI_LIB=$MPI_LOC/lib
export LIBMPI="-L $MPI_LIB -lmpich"
export MPI_INCLUDE=$MPI_LOC/include
export ARMCI_NETWORK=OPENIB
export IB_HOME=/usr
export IB_INCLUDE=$IB_HOME/include
export IB_LIB=$IB_HOME/lib64
export IB_LIB_NAME="-libverbs -libumad -lrdmacm -lpthread -lrt"
export LARGE_FILES
export NWCHEM_MODULES=all


FILENAME=Input
CWD=`pwd`
ulimit -l unlimited

cp $TMPDIR/machines machinefile.$JOB_ID
export NWSCRATCH=/state/partition1
mkdir $NWSCRATCH/$JOB_ID
mkdir $WORK/$JOB_ID
cp $FILENAME.dt $NWSCRATCH/$JOB_ID
cp $FILENAME.nw $NWSCRATCH/$JOB_ID
cp $SGE_O_WORKDIR/machinefile.$JOB_ID  $NWSCRATCH/$JOB_ID
#cp fort.50 $SCRATCH/$JOB_ID
cd $NWSCRATCH/$JOB_ID
echo "/lustre/work/apps/mvapich2-1.6/bin/mpdboot -n `expr $NSLOTS / 12` -f $SGE_O_WORKDIR/machinefile.$JOB_ID" > script1
chmod 755 script1
./script1
run="/lustre/work/apps/mvapich2-1.6/bin/mpirun  -machinefile $SGE_CWD_PATH/machinefile.$JOB_ID -np $NSLOTS /lustre/work/skohale/CPC-D-12-00495a/NWC-V/nwchem-6.1.1/venus-nwchem/ven_nw.e < $FILENAME.dt > Output.dat"

echo run=$run
echo "ulimit -l unlimited " > script
echo $run >>  script
chmod 755 script
date >> $FILENAME.timing
./script
date >> $FILENAME.timing
mpdallexit

cp  Output.dat $WORK/$JOB_ID/Output.dat.$JOB_ID
cp  fort.* $WORK/$JOB_ID/.
cd $NWSCRATCH
rm -rf $JOB_ID
cd $WORK/$JOB_ID
cp Output.dat.$JOB_ID   $CWD/.
cp fort.* $CWD/.
