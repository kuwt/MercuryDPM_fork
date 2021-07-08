#!/bin/bash
#Step 1 create (build) and copy in the centerofmass code
echo
echo
echo Creating centerofmass of mass code and fstatistics...
echo
echo

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

echo
echo
echo Starting centreofmass calculation ...
echo
echo
sleep 10s

for num in 1; do
#for num in 1 4 6 8 10; do
  ../bin/centerofmass ../../Segregation.$num 1.00001
  mv ../../Segregation.$num.com .

echo
echo
echo Starting global stats ...
echo
echo
sleep 10s

	../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 1500 -w 1.0 -z 0 100 -StressTypeForFixedParticles 1 \
	 -o Segregation.$num.stat
echo
echo
echo Starting large particle stats ...
echo
echo
	sleep 10s
	 ../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 1500 -w 1.0 -z 0 100 -StressTypeForFixedParticles 1 \
	 -rmin 1.00001 -o Segregation.$num.large.stat
echo
echo
echo Starting small particle stats ....
echo
echo
	sleep 10s
	../bin/fstatistics ../../Segregation.$num -stattype Z \
	 -n 200 -tmin 1500 -w 1.0 -z 0 100 -StressTypeForFixedParticles 1 \
	 -rmax 1.00001 -o Segregation.$num.small.stat

done