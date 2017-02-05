#!/usr/bin/env bash                                                           
#SBATCH --job-name=jacobiResult                                    
#SBATCH --output=jacobiResult.out                                     
#SBATCH --error=jacobiResult.err                                       
#SBATCH --time=04:00:00                                                      
#SBATCH --exclusive                                                         
#SBATCH --nodes=1 --ntasks-per-node=16                                         
module load slurm/14.11.11 compiler/gcc/5.1.0 intel/mkl/64/11.2/2016.0.0 mpi/openmpi/gcc/1.10.0-tm slurm/14.11.11
make

rm -f data/jacobiResult.dat
touch data/jacobiResult.dat

for nbProcs in 1 2 4 8 16
do
    line="${nbProcs} "
    t=`mpiexec -n ${nbProcs} jacobi_test.exe`
    line="${line} ${t}"
    echo "${line}" >> data/jacobiResult.dat
done