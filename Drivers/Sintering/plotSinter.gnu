#evaluation; create data and fstat files
#file=Sinter2H8;
#for i in 0010 0020 0030 0040 0050 0060 0070 0080 0090 0100 0000; do 
#  ./restartToDataFstat $file.restart.$i; mv $file.data.0000 $file.data.$i;  mv $file.fstat.0000 $file.fstat.$i;
#done

set term png 
set termoption enhanced
set output 'SinterForce.png'
set title 'depth D=0.05, varying time step'
set ylabel 'f^{ind} [\mu N]'
set xlabel '\delta [\mu m]'
set key left top 
min(x,y) = (x < y) ? x : y
max(x,y) = (x > y) ? x : y
p \
'SinterH8.fstat.0000' u ($7*1e6):($9*1e6) t 't=0.0 ms' w p pt 7 ps 1, \
'SinterH8.fstat.0010' u ($7*1e6):($9*1e6) t 't=0.1 ms' w d, \
'SinterH8.fstat.0020' u ($7*1e6):($9*1e6) t 't=0.2 ms' w d, \
'SinterH8.fstat.0030' u ($7*1e6):($9*1e6) t 't=0.3 ms' w d, \
'SinterH8.fstat.0040' u ($7*1e6):($9*1e6) t 't=0.4 ms' w d, \
'SinterH8.fstat.0100' u ($7*1e6):($9*1e6) t 't=0.9 ms' w d, \
max(max(0,1346*x),5*1346*(x-0.11)) notitle lc rgb "black", \
max(-5*1346*x,5*1346*(x-0.11)) notitle lc rgb "black"


set output 'SinterCOM.png'
set title 'depth D=0.05, varying time step'
set xlabel 't [ms]'
set ylabel 'COM_z [\mu m]'
set key left top 
p [0:0.5] 'SinterH10.ene' u ($1*1e3):($8*1e6),  'Sinter2H10.ene' u ($1*1e3):($8*1e6)
