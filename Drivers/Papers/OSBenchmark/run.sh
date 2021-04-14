#!/bin/bash

# Simulation parameters
# - Set the number of omp threads. E.g. 'opt=-omp 12' will use 12 threads per simulation.
# - Use 'opt=-writeOutput' to make the code output additional data that can be used for visualisation and analysis.
opt='-omp 12'

# Run all simulations for the OS paper:

# - Silo with small orifice, M1 particles
./OSSilo -useSmallOrifice -useM1 $opt

# - Silo with large orifice, M1 particles
./OSSilo -useM1 $opt

# - Silo with small orifice, M2 particles
./OSSilo -useSmallOrifice $opt

# - Silo with large orifice, M2 particles
./OSSilo $opt

# - Penetration with 25K particles
./OSPenetration -size 25K $opt

# - Penetration with 50K particles
./OSPenetration -size 50K $opt

# - Penetration with 100K particles
./OSPenetration -size 100K $opt

# - Rotating drum with M1/M2 particles
./OSDrum $opt

# - Rotating drum with softer particles
./OSDrum -soft $opt
