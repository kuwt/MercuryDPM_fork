#/bin/bash
#set -x
make -j 16 -C ~/MercuryLab/serial/Tools/
#ids=('1' '14' '19' '2' '22' '3' '4' '6' '82' '11' '17' '19b' '21' '2b' '30' '5' '8')

## construct data and vtu files of all finished simulations

# loop through all v9 simulations
for id in '1' '14' '19' '2' '22' '3' '4' '6' '82' '11' '17' '19b' '21' '2b' '30' '5' '8' '21b' '22b'; do
	cd v9run${id}
	#only continue if finished
	tmax=`tail -n1 GCG${id}.ene0|awk '{print $1}'`
	finished=`echo "$tmax > 50" | bc`
	if [ $finished = '1' ]; then 
		echo $id finished $tmax $finished
		# if some vtu files are missing, reconstruct the data and vtu files
		if [ ! -f GCG${id}Particle_10000.vtu ]; then 
			rm GCG${id}.data
			~/MercuryLab/serial/Tools/CombineParallelDataFiles GCG${id}; 
			~/MercuryLab/serial/Tools/data2pvd GCG${id}.data GCG${id}Particle
		fi
	else
		echo !!! $id unfinished $tmax $finished
	fi
	cd ..
done
