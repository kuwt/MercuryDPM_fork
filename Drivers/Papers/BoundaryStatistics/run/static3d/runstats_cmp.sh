cd /storage/usr/people/weinhartt/DRIVERS/BoundaryStatistics/run
rm -f static3d/static3d.* 
sc/quick_run static3d  
cd static3d
cmd="../../../FlowRulePaper/statistics_while_running.exe static3d .01 -w .2 -superexact"
$cmd -StressTypeForFixedParticles 3 -o static3d.inf.stat
$cmd
