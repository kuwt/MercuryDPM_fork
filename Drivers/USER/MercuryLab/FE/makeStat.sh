#/bin/bash
#set -x
make -j 16 -C ~/MercuryLab/Build/Drivers/MercuryCG/ MercuryCG

for id in '1' '14' '19' '2' '22' '3' '4' '6' '82' '11' '17' '19b' '21' '2b' '30' '5' '8'; do
	cd v9run$id
	
	#only continue if finished
	tmax=`tail -n1 GCG${id}.ene0|awk '{print $1}'`
	finished=`echo "$tmax > 50" | bc`

	if [ -f GCG${id}.TXZ1.stat ]; then
		tmaxStat=`tail -n1 GCG${id}.TXZ1.stat|awk '{print $1}'`
		finishedStat=`echo "$tmax > 50" | bc`
	else
		finishedStat='0'
	fi

	echo id $id tmax $tmax finished $finished finishedStat $finishedStat
	

	if [ $finished = '1' ] && [ $finishedStat != '1' ]; then 
		#if [ ! -f GCG${id}.restart ]; then 
			cp GCG${id}.restart0 GCG${id}.restart
			sed -ie '{s/numberOfProcessors 16 numberOfDomains 16/numberOfProcessors 1 numberOfDomains 1/g}' GCG${id}.restart
			sed -ie '{s/numberOfProcessors 15 numberOfDomains 15/numberOfProcessors 1 numberOfDomains 1/g}' GCG${id}.restart
			sed -ie '{s/numberOfProcessors 8 numberOfDomains 8/numberOfProcessors 1 numberOfDomains 1/g}' GCG${id}.restart
			sed -ie '{s/numberOfProcessors 7 numberOfDomains 7/numberOfProcessors 1 numberOfDomains 1/g}' GCG${id}.restart
		#fi

		cmd=~/'MercuryLab/Build/Drivers/MercuryCG/MercuryCG GCG'$id' -w 0.005 -timesmooth -function Heaviside -tmax '$tmax
		echo creating stats for id $id
		$cmd -o GCG${id}.T.stat    -coordinates O
		$cmd -o GCG${id}.T1.stat   -coordinates O -species 1
		$cmd -o GCG${id}.TX.stat   -coordinates X  -nx 100
		$cmd -o GCG${id}.TX1.stat  -coordinates X  -nx 100 -species 1
		$cmd -o GCG${id}.TXZ.stat  -coordinates XZ -nx 100 -nz 37 -z -0.036 0.036  
		$cmd -o GCG${id}.TXZ1.stat -coordinates XZ -nx 100 -nz 37 -z -0.036 0.036 -species 1
	fi

	cd ..
done

#rsync -avz *v9run*/GCG*.T*.stat $home:MercuryLab/Presentations/v9stat/
