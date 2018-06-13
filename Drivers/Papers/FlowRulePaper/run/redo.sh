#this is the simulation from the hstop data, 
set -x 
S=$1
dir=$(pwd)/hstopS${S}
exe=hstop_StudyHeightHmaxAngle
rm -rf $dir
mkdir $dir
make $exe -C .. || exit
cp ../$exe.exe $dir 
cd $dir 
echo $S 4 25 24 > $dir/arg
#echo $S 30 20 24 > $dir/arg

echo '
#/bin/bash
set -x
cd '$dir' 
for i in 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 
do
	sleep 2s
	a1=$(awk "{ print $1}" arg)
	a2=$(awk "{ print $2}" arg)
	a3=$(awk "{ print $3}" arg)
	a4=$(awk "{ print $4}" arg)
	echo $argum	
	'$(pwd)/$exe.exe' $a1 $a2 $a3 $a4
done
' > exe
chmod +x ./exe
sleep 3s

~/clusterscriptexecute ./exe


