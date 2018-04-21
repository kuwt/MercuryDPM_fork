#!/usr/bin/env gnuplot
simdir = "~/MercuryDPM/MercuryBuild/Drivers";
drivername = "ReposeTest";
seriesname = "muICal3D-ballotini-single";
probnames = "ballotini-18 ballotini-20 ballotini-22 ballotini-24 ballotini-25 ballotini-26 ballotini-28 ballotini-30 ballotini-32";
#seriesname = "muICal3D-ballotini";
#probnames = "ballotini-20-26";
seriesname = "muICal3D-grit";
probnames = "grit-20 grit-22 grit-24 grit-25 grit-28 grit-30 grit-32 grit-34 grit-35 grit-36 grit-38 grit-40";
i = 12500;

do for [probname in probnames] {
    parsfn = sprintf("%s/%s/%s/%s.pars", simdir, drivername, seriesname, probname);
    enefn = sprintf("%s/%s/%s/%s.ene", simdir, drivername, seriesname, probname);
    datafn = sprintf("%s/%s/%s/%s.data.%d", simdir, drivername, seriesname, probname, i);
    muIfn = sprintf("%s/%s/%s/%s.muI", simdir, drivername, seriesname, probname);

    xmin = -0.56;
    xmax = +0.56;
    ymin = -0.21;
    ymax = +0.21;
    zmin = -0.01;
    zmax = +0.65;
    thetamin = 20;
    thetamax = 26;

    vsc = 2;
    vl = 0.01;
    magn(x,y) = sqrt(x**2 + y**2);
    magn(x,y,z) = sqrt(x**2 + y**2 + z**2);

# 2D plots (side view), colour by velocity components

    if (0) {
        set term pngcairo size 1400,1400;
        set out sprintf("~/public_html/%s/%s/%s-data.png", drivername, seriesname, probname);
        set title sprintf("%s : %s : %s : %d", drivername, seriesname, probname, i);
        set multiplot layout 3,1;

        set cbrange [0:vsc];
        set size ratio -1;
        load "../../Visualisation/colours-rgb.palette";
        set xlabel "x - downstream";
        set ylabel "z - perpendicular";
        set xrange [xmin:xmax];
        set yrange [zmin:zmax];
# plot datafn u 1:3:7 w circles lw 1 lc rgb "gray" title "",\
#     datafn u 1:3:($4*vl/magn($4,$5,$6)):($6*vl/magn($4,$5,$6)):(magn($4,0,$6)) w vectors palette title "";
        plot datafn u 1:3:($4*vl/magn($4,$5,$6)):($6*vl/magn($4,$5,$6)):(magn($4,0,$6)) w vectors palette title "";

        set xlabel "y - crossstream";
        set ylabel "z - perpendicular";
        set xrange [ymin:ymax];
        set yrange [zmin:zmax];
# plot datafn u 2:3:7 w circles lw 1 lc rgb "gray" title "",\
#      datafn u 2:3:($5*vl/magn($4,$5,$6)):($6*vl/magn($4,$5,$6)):(magn(0,$5,$6)) w vectors palette title "";
        plot datafn u 2:3:($5*vl/magn($4,$5,$6)):($6*vl/magn($4,$5,$6)):(magn($4,$5,$6)) w vectors palette title "";

        set xlabel "x - downstream"; 
        set ylabel "y - crossstream";
        set xrange [xmin:xmax];
        set yrange [ymin:ymax];
# plot datafn u 1:2:7 w circles lw 1 lc rgb "gray" title "",\
#      datafn u 1:2:($4*vl/magn($4,$5,$6)):($5*vl/magn($4,$5,$6)):(magn($4,$5,0)) w vectors palette title 
        plot datafn u 1:2:($4*vl/magn($4,$5,$6)):($5*vl/magn($4,$5,$6)):(magn($4,$5,0)) w vectors palette title "";

        set size ratio 0;
        unset multiplot

# colour by radius
# set cbrange [0.003:.007];
# splot datafn u 1:2:3:7 w p palette;

# 3D plot, colour by speed
        set out sprintf("~/public_html/%s/%s/%s-data3d.png", drivername, seriesname, probname);
        set view equal xyz;
        set view 50,15,1;
        set xlabel "x - downstream";
        set ylabel "y - crossstream";
        set zlabel "z - perpendicular";
        set xrange [-.25:.25];
        set yrange [-.1:.1];
        set zrange [0:.3];
        set cbrange [0.0:vsc];
        set colorbox horizontal user origin .1,.05 size .8,.04;

        splot datafn u 1:2:3:(magn($4,$5,$6)) w p palette;

        unset view;
        unset xrange;
        unset yrange; 
        unset zrange;
        unset cbrange;
        unset colorbox;
    }

    set term pngcairo size 1280,700;
    set out sprintf("~/public_html/%s/%s/%s-ene.png", drivername, seriesname, probname);
    set title sprintf("%s : %s : %s", drivername, seriesname, probname);
    set multiplot layout 2,1;
    set xlabel "time";
    set xrange [0:200];
    set xtics 20;
    set ylabel "energy / mass";
    set yrange [0:2.0];
#plot \
#     enefn u 1:($3/$10) w l lw 2 lc 1 title "GPE",\
#     enefn u 1:($4/$10) w l lw 2 lc 2 title "KE",\
#     enefn u 1:($14/$10) w l ls 0 lc 2 lw 2 title "ThKE",\
#     enefn u 1:($5/$10) w l lw 2 lc 3 title "RKE",\
#     enefn u 1:($6/$10) w l lw 2 lc 4 title "EPE";
    plot \
        enefn u 1:($4/$10) w l lw 2 lc 2 title "KE",\
        enefn u 1:($14/$10) w l ls 0 lc 2 lw 1 title "ThKE",\
        enefn u 1:($5/$10) w l lw 2 lc 3 title "RKE",\
        enefn u 1:($6/$10) w l lw 2 lc 4 title "EPE";

    set ylabel "xmom";
    set yrange [0:0.1];
    plot enefn u 1:11 w l lw 2 title "xmom";
    unset multiplot;

    if (0) {
        set out sprintf("~/public_html/%s/%s/%s-muI.png", drivername, seriesname, probname);
        set title sprintf("%s : %s : %s", drivername, seriesname, probname);
        set multiplot layout 3,1
            set xlabel "theta";
        set xrange [thetamin:thetamax];

        set ylabel "xmom";
        set yrange [0:0.05];
        set ytics auto;
        set ytics format "%.2e";
        plot muIfn u 2:4 w l lw 2;

        set ylabel "ke";
        set yrange [0:*];
        plot muIfn u 2:5 w l lw 2;

        set ylabel "heat content";
        set yrange [0:1];
        set ytics format "%.2f";
        plot muIfn u 2:($6/$5) w l lw 2;
    }
}
