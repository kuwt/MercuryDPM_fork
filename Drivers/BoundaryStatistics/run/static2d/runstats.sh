cd /storage/usr/people/weinhartt/DRIVERS/BoundaryStatistics/run
#rm -f static2d/static2d.* 
#sc/quick_run static2d  
#make  -C .. statXZ
#make  -C ../../FlowRulePaper statistics_while_running
cd static2d
w=0.125
cmd="../../../FlowRulePaper/statistics_while_running.exe static2d .02 -w "$w" -n 200  -StressTypeForFixedParticles 2 -o static2d.z.stat -superexact"
#$cmd
cmd="../../statXZ.exe static2d .02 -StressTypeForFixedParticles 1 -w "$w" -n 200 -o static2d.tra.stat -superexact"
$cmd
cmd="../../statXZ.exe static2d .02 -StressTypeForFixedParticles 3 -w "$w" -n 200 -o static2d.inf.stat -superexact"
#$cmd
