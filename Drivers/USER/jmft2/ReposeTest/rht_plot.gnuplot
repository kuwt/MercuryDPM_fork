#!/usr/bin/env gnuplot

istart = 0;
iend = 0;
ijump = 5;
isample = 60;

unset key;

problemname = "ReposeHeapTest-run2";
problemname = "ReposeHeapTest-binbinlarge-run0";
problemname = "ReposeHeapTest-elasticparticles-run3";
indir = "~/MercuryDPM/MercuryBuild/Drivers/ReposeTest/";
enefn = indir.problemname.".ene";
haffn = indir.problemname.".haf";
datafn = indir.problemname.".data.";

outname = "rht-hra";
outname = problemname;
outdir = "~/MercuryDPM/MercuryBuild/Drivers/ReposeTest/";
outdir = "~/public_html/ReposeTest-results/";

set term png size 1024,768;
set output outdir.outname.".png";
set multiplot layout 3,2;

set grid;

set yrange [0:10];
plot haffn u 1:($2 * 180 / pi) w l;

set yrange [*:*];
set logscale y;
set ytics format "10^{%T}";
plot haffn u 1:(abs($3) * 180 / pi) w l; 
unset logscale;

unset title;
set key bottom right;
set yrange [-.1:.1];
set ytics format "%.2f";
plot enefn u 1:3 w l lc rgb "green" title "PE";

set yrange [*:*];
set logscale y;
set ytics format "10^{%T}";
plot enefn u 1:($4-$14) w l lc rgb "black" title "bulk KE",\
     enefn u 1:14 w l lc rgb "red" title "thermal KE",\
     enefn u 1:4 w l lc rgb "blue" title "total KE";

#set title "ela";
#plot enefn u 1:6 w l lc rgb "orange" title "EPE";

set yrange [*:*];
set ytics format "10^{%T}";
unset title;
set key top left;
plot enefn u 1:($4/$14 - 1) w l lc rgb "black" title "(bulk ke)/(thermal ke)";
unset logscale;
set ytics format "%.2f";

#set key bottom right;
#set title "total energy";
#plot enefn u 1:($3+$4+$5+$6) w l lc rgb "black" title "total energy";
#unset key;

set yrange [*:*];
set key bottom right;
plot enefn u 1:13 w l lc rgb "black" title "pz";
unset grid;

unset multiplot;

#set term png size 1280,512;
#set output outdir.outname."-slopes.png";
#set title datafn.isample;
#set xrange [-2.5:2.5];
#set yrange [0:.5];
#set size ratio -1;
#set mxtics;
#set mytics;
#set grid mxtics mytics;
#print datafn.isample;
#plot datafn.isample u 1:3 w d;
#unset grid;
#
#set term gif animate delay 0 size 800,600;
#set output outdir.outname.".gif";
#set view equal xyz;
#
#do for [i=istart:iend:ijump] {
#    set title datafn.i
#    set xrange [-5:5];
#    set yrange [-4:4];
#    # splot datafn.i u 1:2:3;
#    set size ratio -1;
#    plot datafn.i u 1:3:7 w circles;
#    #velscale = 5;
#    #plot datafn.i u 1:3:($4/velscale):($6/velscale) w vectors;
#}
#
