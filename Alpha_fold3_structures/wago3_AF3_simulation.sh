#!/bin/bash
#-----------------------------------------------------------------
# Example SLURM job script to run GROMACS with MPI on MOGON.
#
# This script requests 128 cores on two node. The job
# will have access to all the memory in the nodes.  
#-----------------------------------------------------------------

#SBATCH -J Wago3_simulation          # Job name
#SBATCH -o Wago3_simulation.out   # Specify stdout output file (%j expands to jobId)
#SBATCH -e Wago3_simulation.err   # Specify stderr output file (%j expands to jobId)
#SBATCH -p parallel              # Partition/Queue name
#SBATCH -C skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N 8                     # Total number of nodes requested (64 cores/node)
#SBATCH -t 96:00:00              # Run time (hh:mm:ss) - 0.5 hours

#SBATCH -A nhr-dyndisphase             # Specify allocation to charge against

# Load all necessary modules if needed (these are examples)
# Loading modules in the script ensures a consistent environment.
module load bio/GROMACS/2018.1-intel-2018.02 #you can select a specific version, too

# Launch the executable

mkdir equil
cd equil || exit

gmx grompp -r ../em.gro -c ../em.gro  -p ../topol.top  -f ../equilibration-W.mdp  -o equilibration-W.tpr -maxwarn 3

srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-W 

gmx grompp -r equilibration-W.gro -c equilibration-W.gro  -p ../topol.top  -f ../equilibration-nvt.mdp  -o equilibration-nvt.tpr -maxwarn 3

srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-nvt

gmx grompp -r equilibration-nvt.gro -c equilibration-nvt.gro  -p ../topol.top  -f ../equilibration-npt.mdp  -o equilibration-npt.tpr -maxwarn 3

srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-npt

cd ../ 
mkdir dynamics
cd dynamics || exit

gmx grompp -r ../equil/equilibration-npt.gro -c ../equil/equilibration-npt.gro  -p ../topol.top  -f ../dynamics.mdp  -o dynamics.tpr -maxwarn 3

srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -cpi dynamics.cpt -v -deffnm dynamics

# Extend from 10-20 microseconds

#gmx convert-tpr -s dynamics.tpr -extend 1000000000 -o dynamics_10-20ms.tpr

#srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -cpi dynamics.cpt -deffnm dynamics_10-20ms -noappend