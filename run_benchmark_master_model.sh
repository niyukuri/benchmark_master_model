#!/bin/bash -l
#
#PBS -N benchmarkMM
#PBS -P CBBI0949
#PBS -q smp
#PBS -l select=1:ncpus=24:mpiprocs=24:nodetype=haswell_reg
#PBS -l walltime=24:00:00
#PBS -m be
#PBS -M niyukuri@aims.ac.za



ulimit -s unlimited

module add /apps/chpc/scripts/modules/bio/app/simpact-cyan/0.21.0
module add /apps/chpc/scripts/modules/bio/app/Seq-Gen/1.3.4
module add /apps/chpc/scripts/modules/bio/app/FastTree/2.1.10
module add chpc/R/3.4.4-gcc7.2.0


cd /mnt/lustre/users/dniyukuri/benchmark_master_model


Rscript wrapper.benchmark.master.model.R
