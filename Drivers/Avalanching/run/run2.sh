#!/bin/bash

make -C .. ContinuousTilt
make -C .. SteppedTilt
mkdir -p Test2
cp ../ContinuousTilt.exe ../SteppedTilt.exe ../ContinuousTilt.cpp ../SteppedTilt.cpp Test2
cd Test2

for mu in 0 0.1 0.1 0.2; do
for mur in 0 0.02; do

c=Mu${mu}Mur${mur}
nohup nice -19 \
./ContinuousTilt.exe -name ContinuousTilt$c -mu $mu -murolling $mur &

done
done
