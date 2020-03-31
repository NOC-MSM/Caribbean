#!/bin/bash
#SBATCH -N 11
#SBATCH --job-name=O1_AGRIF
#SBATCH --exclude=node[42]
echo node42 excluded
#SBATCH --time=3:00:00

scontrol wait_job $SLURM_JOB_ID

module purge
export MODULEPATH=/home/acc/MyMods:$MODULEPATH
module use /home/acc/MyMods
module load nemo-PrgEnv/4.0
module list

#
  export OMP_NUM_THREADS=1
  export OCORES=156 #96
  export XCORES=2 # 13 # 8 # 2
#
  export EXE_DIR=$PWD
#
# end of set up
###############################################################
#
# change to the working directory
#
  date
  squeue

  cd $EXE_DIR
       echo time `which mpirun` --report-bindings -x MALLOC_MMAP_MAX_=-1 \
                                -x MALLOC_TRIM_THRESHOLD_=33554432       \
                                -np $XCORES ./xios_server.exe            \
                              : -np $OCORES --map-by node  --mca mpi_paffinity_alone 1 ./nemo
#
            time `which mpirun` --report-bindings -x MALLOC_MMAP_MAX_=-1 \
                                -x MALLOC_TRIM_THRESHOLD_=33554432       \
                                -np $XCORES ./xios_server.exe            \
                              : -np $OCORES --map-by node  --mca mpi_paffinity_alone 1 ./nemo
#

