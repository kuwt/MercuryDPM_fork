#!/usr/bin/gnuplot

expname="Exp9";
indir="/home/bj268/build/".expname."/";
set term pngcairo size 800,600
set xrange [0:9]
set yrange [-4:4];
set cbrange [0:1]
set palette model RGB defined (0 'red', 1 'blue')
set size ratio -1;
unset key;
unset colorbox;

set output "~/public_html/binbin".expname."-top600.png";

plot indir.'Exp.data.600' u 1:3:($14+1)/5:14 every ::2 w p pal pt 7 ps var tit "600"
