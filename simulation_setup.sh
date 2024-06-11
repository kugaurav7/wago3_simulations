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
#SBATCH -N 1                              # Total number of nodes requested (64 cores/node)
#SBATCH -t 10:00:00                       # Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A nhr-dyndisphase                # Specify allocation to charge against

# Load the gromacs module
module load bio/GROMACS/2018.1-intel-2018.02

#Input paramer
input_pdb="mut16_rbr_ffr_cg.pdb"
output_prefix="output"
salt="0.15"
ffdir="martini_v300"
minmdp="minimization.mdp"

#Convert pdb to gro
gmx pdb2gmx -f "${input_pdb}" -o "${output_prefix}.gro" -water spce <<EOF
1
EOF

#Put protein in the box
gmx editconf -f "${output_prefix}.gro" -o "${output_prefix}_box.gro" -c -d 10.0

#Edit the topology 

#Solvate the protein
gmx solvate -cp "${output_prefix}_box.gro" -cs spc216.gro -o "${output_prefix}_solv.gro" -p topol.top

#Add ions
gmx grompp -f ions.mdp -c "${output_prefix}_solv.gro" -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o "${output_prefix}_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral -conc "${salt}" <<EOF
SOL 
EOF

#Energy Minimization
gmx grompp -f "$minmdp" -p topol.top -c "${output_prefix}_solv_ions.gro" -o em.tpr -pp all_PRO.top -maxwarn 3 -r "${output_prefix}_solv_ions.gro"
gmx mdrun -v -deffnm em


