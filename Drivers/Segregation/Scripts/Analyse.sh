#!/bin/bash
#Step 1 create (build) and copy in the centerofmass code
echo Creating centerofmass of mass code
(cd ../../MercuryCG && make  centerofmass && make fstatistics)
mkdir -p bin
cp ../../MercuryCG/centerofmass bin/
cp ../../MercuryCG/fstatistics bin/

# Step 2 create stat and com files
mkdir -p stat
#set -x
cd stat
#for num in 1 5 7 9; do
#for num in 2 3; do
for num in 1; do
#for num in 1 4 6 8 10; do
	restartfile=$(ls ../../Segregation.$num.restart)
  ../bin/centerofmass ../../Segregation.$num 0.0003001

	../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 90 -w 0.0003 -z -0.002 0.016 -StressTypeForFixedParticles 1
	sleep 10s
	 ../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 90 -w 0.0003 -z -0.002 0.016 -StressTypeForFixedParticles 1 \
	 -rmin 0.0003001 -o Segregation.$num.large.stat
	sleep 10s
	../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 90 -w 0.0003 -z -0.002 0.016 -StressTypeForFixedParticles 1 \
	 -rmax 0.0003001 -o Segregation.$num.small.stat
	sleep 10s
done