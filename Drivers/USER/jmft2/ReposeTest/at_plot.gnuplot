#!/usr/bin/env gnuplot

anim = 0;
istart = 3000;
iend = 5000;
ijump = 8;
isample = 60;

unset key;

#problemname = "AvalancheTest"
problemname = "AvalancheTest-larger"
#indir = "~/MercuryDPM/MercuryBuild/Drivers/ReposeTest/AvalancheTest-preliminary/";
indir = "~/MercuryDPM/MercuryBuild/Drivers/ReposeTest/";
enefn = indir.problemname.".ene";
haffn = indir.problemname.".haf";
datafn = indir.problemname.".data.";

outname = problemname;
outdir = "~/MercuryDPM/MercuryBuild/Drivers/ReposeTest/";
outdir = "~/public_html/ReposeTest-results/";

set term png size 1024,768;
set output outdir.outname.".png";
set multiplot layout 2,2;

set grid;
unset title;
set key bottom right;

set yrange [*:*];
set ylabel "N";
set y2range [0:*];
set y2tics;
set y2label "mass";
plot enefn u 1:2 w l axes x1y1 lc rgb "blue" title "N",\
     enefn u 1:10 w l axes x2y2 lc rgb "black" title "mass";
unset y2tics;
unset y2label;

set yrange [*:*];
set ylabel "EPE";
plot enefn u 1:6 w l lc rgb "green" title "EPE";

set yrange [0:*];
#set logscale y;
#set ytics format "10^{%T}";
set ylabel "KE";
plot enefn u 1:($4-$14) w l lc rgb "black" title "bulk KE",\
     enefn u 1:14 w l lc rgb "red" title "thermal KE",\
     enefn u 1:4 w l lc rgb "blue" title "total KE";
#unset logscale;

set key top left;
set key box;
set yrange [0:25]
set ylabel "heap angle /deg";
set ytics format "%.1f";
plot haffn u 1:3 w l lc rgb "black" title "heap angle";
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

if (anim) {
    set term gif animate delay 0 size 800,600;
    set output outdir.outname.".gif";
#    set lmargin 5;
#    set rmargin 5;
#    set tmargin 0;
#    set bmargin 0;
    unset key;

    #set xrange [-1.2:+1.2];
    #set yrange [-1.2:+1.2];
    #set zrange [-.1:1.1];
    #set view equal xyz;
    #set view 80,45;
    #do for [i=istart:iend:ijump] {
    # set title datafn.i
    # splot datafn.i u 1:2:3:7 w dots;
    #}

    do for [i=istart:iend:ijump] {
        set title datafn.i;
        # set multiplot layout 2,1;

        set xrange [-1.2:1.2];
        set yrange [-.1:1.1];
        set size ratio -1;
        plot datafn.i u 1:3:7 w circles;

        # set xrange [-1.2:1.2];
        # set yrange [-1.2:1.2];
        # set cbrange [0:1];
        # set size ratio -1;
        # plot datafn.i u 1:2:3 w p palette;

        # unset multiplot;
    }
}
