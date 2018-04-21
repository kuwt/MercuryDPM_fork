#!/usr/bin/env gnuplot
#set term png size 1440,1024;
#set out "~/public_html/CT2D.png";
set term pdfcairo size 24cm,18cm;
set out "~/public_html/CT2D.pdf";
# set multiplot layout 4,2;
set multiplot layout 4,1;
# set multiplot layout 2,2;
do for [probname in "ballotini sand grit frictionless"] {
    set title probname;
    set key top right opaque box;
    set xrange [0:1200];
    set xtics 50;
    set grid;
#    set yrange [0:1];
#    set logscale y;
    set xlabel "time";
    set ylabel "energy";
    enefn = "CT2D-".probname."/".probname.".ene";
#    plot enefn u 1:3 w l title "GPE",\
#         enefn u 1:4 w l title "TKE",\
#         enefn u 1:14 w l title "ThKE";
    plot enefn u 1:3 w l title "GPE" lc rgb "red" ,\
         enefn u 1:14 w l title "ThKE" lc rgb "blue";
#    plot enefn u 8:11 w l; 
    

#    set yrange [1e-3:1e3];
#    set logscale y;
#    set ylabel "speed";
#    plot "CT2D-".probname."/".probname.".ene" u 1:($11/$10) w l title "speed";
#    unset logscale;
}
unset multiplot;
