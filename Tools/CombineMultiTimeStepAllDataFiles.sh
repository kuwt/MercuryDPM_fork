#!/bin/bash
#This script will merge time steps together assume 16 cores. Use ./AntCombine filename
#
#
#First step complier the files across process to create one file per time step

#Work out the directory location of where the script is
dir="$(dirname $0)"

(cd $dir && make CombineParallelDataFiles)

for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
	rm $1.data$i
	echo
	echo
	echo "Doing processor $i"
	echo "------------------"
	sleep 1s
	for (( timeStep =0; timeStep<=1007; timeStep++ ))
		do
		echo "Doing timestep $timeStep"
        	cat $1.data$i.$timeStep >> $1.data$i
	
	done
	
done

echo Combined process files to one data to rule them all
$dir/CombineParallelDataFiles $1
