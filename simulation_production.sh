#!/bin/bash
#-----------------------------------------------------------------
# Example SLURM job script to run GROMACS with MPI on MOGON.
#
# This script requests 128 cores on two node. The job
# will have access to all the memory in the nodes.  
#-----------------------------------------------------------------

#SBATCH -J wago3_simulation       # Job name
#SBATCH -o wago3_simulation.out   # Specify stdout output file (%j expands to jobId)
#SBATCH -e wago3_simulation.err   # Specify stderr output file (%j expands to jobId)
#SBATCH -p parallel                       # Partition/Queue name
#SBATCH -C skylake                        # select either 'broadwell' or 'skylake'
#SBATCH -N 8                              # Total number of nodes requested (64 cores/node)
#SBATCH -t 96:00:00                       # Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A nhr-dyndisphase                # Specify allocation to charge against

# Load the gromacs module
module load bio/GROMACS/2018.1-intel-2018.02

#Input paramer
input_pdb="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/em.gro"
equilibration_W_mdp="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/equilibration_W.mdp"
equilibration_NVT_mdp="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/equilibration_NVT.mdp"
equilibration_NPT_mdp="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/equilibration_NPT.mdp"
dynamics_mdp="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/dynamics.mdp"
topology="/lustre/miifs01/project/nhr-dyndisphase/kugaurav/wago_simulations/wago3_simulations/wago3_without_RNA/starting_simulation/topol.top"

#Equilibrate the system
mkdir equil
cd equil || exit

#Equilibration NVT with protein restrained
gmx grompp -f "${equilibration_W_mdp}" -p "{$topology}" -c "${input_pdb}" -o equilibration_W.tpr -maxwarn 3 -r "${input_pdb}"
srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-W

#Equilibration NVT without protein restrained
gmx grompp -f "${equilibration_NVT_mdp}" -p "{$topology}" -c "equilibration_W.gro" -o equilibration_NVT.tpr -maxwarn 3 -r "equilibration_W.gro"
srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-NVT

#Equilibration NPT without protein restrained
gmx grompp -f "${equilibration_NPT_mdp}" -p "{$topology}" -c "equilibration_NVT.gro" -o equilibration_NPT.tpr -maxwarn 3 -r "equilibration_NVT.gro"
srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm equilibration-NPT

cd ../

#Run dynamics
mkdir dynamics
cd dynamics || exit

gmx grompp -f "${dynamics_mdp}" -p "{$topology}" -c "../equil/equilibration_NPT.gro" -o dynamics.tpr -maxwarn 3 -r "../equil/equilibration_NPT.gro"
srun -n 256 -c 2 gmx_mpi mdrun -ntomp 2 -v -deffnm dynamics

