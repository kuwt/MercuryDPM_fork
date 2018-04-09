#this is the simulation from the BoundaryPaper
dir=H10A26L2M0.5B0.5NEW
rm -rf $dir
mkdir $dir
make flowRule_studyHeightAngle -C ../../FlowRulePaper
cp ../../FlowRulePaper/flowRule_studyHeightAngle.exe $dir
cd $dir
~/clusterscriptexecute ./flowRule_studyHeightAngle.exe 5 10 26 -tmax 3000
#./flowRule_studyHeightAngle.exe 5 10 26 -tmax 1e-4 -oldValues
