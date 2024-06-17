#!/bin/bash
#-----------------------------------------------------------------
# Example SLURM job script to run GROMACS with MPI on MOGON.
#
# This script requests 128 cores on two node. The job
# will have access to all the memory in the nodes.  
#-----------------------------------------------------------------

#SBATCH -J wago3_AF3_setup       # Job name
#SBATCH -o wago3_AF3_setup.out   # Specify stdout output file (%j expands to jobId)
#SBATCH -e wago3_AF3_setup.err   # Specify stderr output file (%j expands to jobId)
#SBATCH -p parallel                       # Partition/Queue name
#SBATCH -C skylake                        # select either 'broadwell' or 'skylake'
#SBATCH -N 1                              # Total number of nodes requested (64 cores/node)
#SBATCH -t 02:00:00                       # Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A nhr-dyndisphase                # Specify allocation to charge against

# Load the gromacs module
module load bio/GROMACS/2018.1-intel-2018.02

#parameters
output_prefix="output"
ionmdp="../parameter_files/ions.mdp"
salt="0.15"
minmdp="../parameter_files/em.mdp"

for i in {1..5}
do
    mkdir simulation_"$i"
    cd simulation_"$i" || exit

    #Copy the force field directory
    cp -r ../parameter_files/amber99sb-star-ildn-q-OL3-tip4pd.ff ./

    #Convert pdb to gro
    gmx pdb2gmx -f "../parameter_files/wago3_AF3_${i}.pdb" -o "${output_prefix}.gro" -water tip3p -ignh <<EOF
1
EOF

    #Put protein in the box
    gmx editconf -f "${output_prefix}.gro" -o "${output_prefix}_box.gro" -c -d 2.0

    #Solvate the protein
    gmx solvate -cp "${output_prefix}_box.gro" -cs spc216.gro -o "${output_prefix}_solv.gro" -p topol.top

    #Add ions
    gmx grompp -f "$ionmdp" -c "${output_prefix}_solv.gro" -p topol.top -o ions.tpr
    gmx genion -s ions.tpr -o "${output_prefix}_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral -conc "${salt}" <<EOF
SOL 
EOF

#Energy Minimization
gmx grompp -f "$minmdp" -p topol.top -c "${output_prefix}_solv_ions.gro" -o em.tpr -pp all_PRO.top -maxwarn 3 -r "${output_prefix}_solv_ions.gro"
gmx mdrun -v -deffnm em

cd ..

done