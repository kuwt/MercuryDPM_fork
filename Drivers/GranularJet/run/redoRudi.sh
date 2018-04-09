#!/bin/bash
#rm -rf redoRudi
mkdir -p redoRudi
cd redoRudi
make -C ../.. GranularJet && cp ../../GranularJet.exe .


#redo Rudi's small simulation
~/bin/clusterscriptexecute ./GranularJet.exe -tmax 6 -number_of_saves 200

#redo Rudi's large simulation
~/bin/clusterscriptexecute ./GranularJet.exe -chuteAngle 25.4 -chuteLength 0.6 -chuteWidth 0.25 -funHf 0.43 -funD 0.017 -funnz 50 -tmax 6 -number_of_saves 200

#new large simulation
~/bin/clusterscriptexecute ./GranularJet.exe -max_failed 10 -chuteLength 0.3 -chuteWidth 0.25 -funHf 0.43 -funD 0.017 -funnz 50 -tmax 20 -number_of_saves 400

#new large simulation
~/bin/clusterscriptexecute ./GranularJet.exe -max_failed 10 -chuteLength 0.5 -chuteWidth 0.25 -funHf 0.43 -funD 0.02 -funnz 75 -tmax 20 -number_of_saves 400

#next time: add -options_data 2
