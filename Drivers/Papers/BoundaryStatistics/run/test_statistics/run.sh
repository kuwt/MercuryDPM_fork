direc=$(pwd)
cd ..
sc/quick_run static2d &&
cd $direc
make statXY -C ../.. &&
for type in 0 1 2 3
do
mkdir $type
cd $type
../../../statXY.exe ../../static2d/static2d 0.01 -StressTypeForFixedParticles $type -options_data 1 -n 200 -w 0.15 -superexact
cd ..
done
