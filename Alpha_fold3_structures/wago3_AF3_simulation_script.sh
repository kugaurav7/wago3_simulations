#!/bin/bash

#A for loop to run simulation in each folders
for i in {1..5}
do
    cd simulation_"$i" || exit

    #Run the simulation
    sbatch wago3_AF3_simulation.sh

    cd ..
done