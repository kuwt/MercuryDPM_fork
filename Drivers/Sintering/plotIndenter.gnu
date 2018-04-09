#D=0.1; for N in 0 25 50 75 100 125 150; do ./Indenter $N $D; done; N=100; for D in 0.2 0.4 0.8 1.6; do ./Indenter $N $D; done

set term png
set output 'sphericalIndenterVaryN.png'
set title 'depth D=0.1, varying time step'
set ylabel 'f^{ind} [\mu N]'
set xlabel '\delta [\mu m]'
set key left top 
p [-0:0.12][-30:80] \
'sphericalIndenterN000D0.1.ene' u (-$10*1e6):($9*1e6) t 'N000' w l, \
'sphericalIndenterN025D0.1.ene' u (-$10*1e6):($9*1e6) t 'N025' w l, \
'sphericalIndenterN050D0.1.ene' u (-$10*1e6):($9*1e6) t 'N050' w l, \
'sphericalIndenterN075D0.1.ene' u (-$10*1e6):($9*1e6) t 'N075' w l, \
'sphericalIndenterN100D0.1.ene' u (-$10*1e6):($9*1e6) t 'N100' w l, \
'sphericalIndenterN125D0.1.ene' u (-$10*1e6):($9*1e6) t 'N125' w l, \
'sphericalIndenterN150D0.1.ene' u (-$10*1e6):($9*1e6) t 'N150' w l, \
1346*x*(x>0) notitle lc rgb "gray", \
-5*1346*x*(x>0) notitle lc rgb "gray", \
5*1346*(x-0.11) notitle lc rgb "gray"

set output 'sphericalIndenterVaryD.png'
set key left top 
set title 'time step N=100'
set ylabel 'f^{ind}/(k1*Depth)'
set xlabel '\delta/Depth'
#p 'sphericalIndenter.ene' u (-$10*1e6):($9*1e6) w lp
p [-0:1][-.100:1] \
'sphericalIndenterN100D0.1.ene' u (-$10*1e6/0.1):($9*1e6/0.1/1346) t 'D0.1' w l, \
'sphericalIndenterN100D0.2.ene' u (-$10*1e6/0.2):($9*1e6/0.2/1346) t 'D0.2' w l, \
'sphericalIndenterN100D0.4.ene' u (-$10*1e6/0.4):($9*1e6/0.4/1346) t 'D0.4' w l, \
'sphericalIndenterN100D0.8.ene' u (-$10*1e6/0.8):($9*1e6/0.8/1346) t 'D0.8' w l, \
'sphericalIndenterN100D1.6.ene' u (-$10*1e6/1.6):($9*1e6/1.6/1346) t 'D1.6' w l, \
x*(x>0) notitle lc rgb "gray", \
-5*x*(x>0) notitle lc rgb "gray", \
5*(x-1) notitle lc rgb "gray"

set output 'sphericalIndenterN0VaryD.png'
set key left top 
set title 'time step N=0'
set ylabel 'f^{ind}/(k1*Depth)'
set xlabel '\delta/Depth'
#p 'sphericalIndenter.ene' u (-$10*1e6):($9*1e6) w lp
p [-0:1][-.00001:.0001] \
'sphericalIndenterN000D0.1.ene' u (-$10*1e6/0.1):($9*1e6/0.1/1346) t 'D0.1' w l, \
'sphericalIndenterN000D0.2.ene' u (-$10*1e6/0.2):($9*1e6/0.2/1346) t 'D0.2' w l, \
'sphericalIndenterN000D0.4.ene' u (-$10*1e6/0.4):($9*1e6/0.4/1346) t 'D0.4' w l, \
'sphericalIndenterN000D0.8.ene' u (-$10*1e6/0.8):($9*1e6/0.8/1346) t 'D0.8' w l, \
'sphericalIndenterN000D1.6.ene' u (-$10*1e6/1.6):($9*1e6/1.6/1346) t 'D1.6' w l, \
x*(x>0) notitle lc rgb "gray", \
-5*x*(x>0) notitle lc rgb "gray", \
5*(x-1) notitle lc rgb "gray"

set output 'sintering.png'
set key left top 
set title 'Sintering process'
set ylabel 'COM_z [\mu m]'
set xlabel 'time [ms]'
p 'Sinter.ene' u ($1*1e3):($8*1e6)
