#rm -f H* 
t=.01
for f in H10A12L0M0.5B0.5 H10A20L2M0.5B0.5
do
for StressType in 1 #0 2 3
do
for w in 0.2
do 
../../../FlowRulePaper/statistics_while_running.exe ../stats2/restart/$f $t -w $w -n 500 \
-superexact -StressTypeForFixedParticles $StressType -o ${f}W${w}Stress${StressType}.stat &
sleep 2s
cp $f.restart ${f}W${w}Stress${StressType}.restart
cp $f.ene ${f}W${w}Stress${StressType}.ene
done
done
done
