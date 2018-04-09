#rm -f H* 
t=10
for f in H20A14L0M0.5B0.5 H20A26L1M0.5B0.5
do
for StressType in 0 2 3
do
for w in 0.5
do 
../../../FlowRulePaper/statistics_while_running.exe restart/$f $t -w $w -n 1000 \
-superexact -StressTypeForFixedParticles $StressType -o ${f}W${w}Stress${StressType}.stat &
sleep 2s
cp $f.restart ${f}W${w}Stress${StressType}.restart
cp $f.ene ${f}W${w}Stress${StressType}.ene
done
done
done
