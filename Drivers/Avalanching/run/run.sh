#!/bin/bash

make -C .. ContinuousTilt
make -C .. SteppedTilt
mkdir -p Test
cp ../ContinuousTilt.exe ../SteppedTilt.exe ../ContinuousTilt.cpp ../SteppedTilt.cpp Test
cd Test

for mu in 0.1; do
for mur in 0 0.02; do

c=Mu${mu}Mur${mur}
nohup nice -19 \
./ContinuousTilt.exe -name ContinuousTilt$c -mu $mu -murolling $mur &
#nohup nice -19 \
#./SteppedTilt.exe -name SteppedTilt$c -mu $mu -murolling $mur &

done
done
