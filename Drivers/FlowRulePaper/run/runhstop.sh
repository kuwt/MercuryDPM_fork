#!/bin/bash
#Note: this is how you run hstop curves (Aug-17 2012)
exe=./hstop_StudyHeightHmaxAngle.exe
for S in 48 49 50 51; do
 mkdir S$S
 cd S$S
 cp ../hstop_StudyHeightHmaxAngle/$exe .
 #nohup nice -19 
 echo $S 4 60 24 > arg
 echo '
#/bin/bash
set -x
cd '$(pwd)' 
for i in 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 
do
	sleep 2s
	a1=$(awk "{ print $1}" arg)
	echo '$exe' $a1 -test
	'$exe' $a1 
done
' > exe
 chmod +x exe
 ~/bin/clusterscriptexecute ./exe
 cd ..
done
