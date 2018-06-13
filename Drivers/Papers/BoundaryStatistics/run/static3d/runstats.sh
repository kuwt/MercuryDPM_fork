cd /storage/usr/people/weinhartt/DRIVERS/BoundaryStatistics/run
#rm -f static3d/static3d.* 
sc/quick_run static3d  
cd static3d
cmd="../../../FlowRulePaper/statistics_while_running.exe static3d .5 -w .2 -n 400 -superexact "
$cmd -StressTypeForFixedParticles 1 -o static3d.1.stat &
$cmd -StressTypeForFixedParticles 0 -o static3d.0.stat &
#$cmd
