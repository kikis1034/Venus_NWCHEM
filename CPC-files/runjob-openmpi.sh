#!/bin/bash
#PBS -N VENUS-ompi
#PBS -q short_eth
#PBS -l nodes=2:ppn=4
#PBS -l walltime=23:59:00
#PBS -j oe

module add openmpi/gcc/64/1.3.3
module add shared
module add torque/2.3.7
module add maui/3.2.6p21
module add intel/compiler/64

cd /home/carolp/work/CPC-D-12-00495-REV/NWC/Input-Output
mpirun -np 8 -machinefile $PBS_NODEFILE ./ven_nw.e < Input.dt > Output.dat

