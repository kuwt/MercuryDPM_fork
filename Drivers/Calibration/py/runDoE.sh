#!/bin/bash
# run a parameter study based on materialA
rm -f runDoE2.sh

for Exe in ./CalibrationHeap ./CalibrationDrum ./CalibrationShearCell; do
rc0=0.5
mus0=0.5
mur0=0.1
bo0=0.0
Exe='nohup '$Exe' -speciesType MaterialA -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1500 -psd cumulative volume diameter 2e-4 0 2e-4 0.5 3e-4 0.5 3e-4 1 -collisionTime 0.000067997 -torsionFriction 0 -normalStress 1000 800 600 400 200 -output '
echo $Exe

for rc in 0.3 0.5 0.8; do
  ExeParam=$Exe' -restitutionCoefficient '$rc' -slidingFriction '$mus0' -rollingFriction '$mur0' -bondNumber '$bo0' -param '_${rc}_${mus0}_${mur0}_${bo0}' &'
  echo $ExeParam >> runDoE2.sh
done
echo 'sleep 10' >> runDoE2.sh

for mus in 0 1; do
  ExeParam=$Exe' -restitutionCoefficient '$rc0' -slidingFriction '$mus' -rollingFriction '$mur0' -bondNumber '$bo0' -param '_${rc0}_${mus}_${mur0}_${bo0}' &'
  echo $ExeParam >> runDoE2.sh
done
echo 'sleep 10' >> runDoE2.sh

for mur in 0 0.5; do
  ExeParam=$Exe' -restitutionCoefficient '$rc0' -slidingFriction '$mus0' -rollingFriction '$mur' -bondNumber '$bo0' -param '_${rc0}_${mus0}_${mur}_${bo0}' &'
  echo $ExeParam >> runDoE2.sh
done
echo 'sleep 10' >> runDoE2.sh

for bo in 1 5 25; do
  ExeParam=$Exe' -restitutionCoefficient '$rc0' -slidingFriction '$mus0' -rollingFriction '$mur0' -bondNumber '$bo' -param '_${rc0}_${mus0}_${mur0}_${bo}' &'
  echo $ExeParam >> runDoE2.sh
done
echo 'sleep 10' >> runDoE2.sh

done

scp runDoE2.sh $cloud:Code/Calibration/build/Drivers/Calibration
