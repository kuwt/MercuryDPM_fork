direc=$(pwd)
cd ..
cd $direc
make -C ../../../FlowRulePaper/ statistics_while_running
for type in 0 1 2 3
do
mkdir $type
cd $type
../../../../FlowRulePaper/statistics_while_running.exe ../H10A20L1M0.5B0.5 0.003 -StressTypeForFixedParticles $type -options_data 1 -nz 200 -w 0.15 -superexact
cd ..
done
