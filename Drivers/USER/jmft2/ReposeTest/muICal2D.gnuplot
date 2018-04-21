#!/usr/bin/env gnuplot
simdir = "~/MercuryDPM/MercuryBuild/Drivers";
drivername = "ReposeTest";
seriesname = "muICal2D-ballotini";
speciesname = "ballotini";
# probnames = "ballotini-35 ballotini-36 ballotini-37 ballotini-40 ballotini-41 ballotini-42 ballotini-43 ballotini-44 ballotini-45"

# Plotting results for mu(I) calibration tests, using muICal2D.
# do for [probname in probnames] {
do for [i = 30:45] {
    probname = sprintf("%s-%d", speciesname, i);
    parsfn = sprintf("%s/%s/%s/%s.pars", simdir, drivername, seriesname, probname);
    enefn = sprintf("%s/%s/%s/%s.ene", simdir, drivername, seriesname, probname);
    muIfn = sprintf("%s/%s/%s/%s.muI", simdir, drivername, seriesname, probname);

    print enefn;
     # set term pdfcairo size 18cm,12cm;
    set term pngcairo size 1280,700;
     
    set out sprintf("~/public_html/%s/%s/%s-ene.png", drivername, seriesname, probname);
    set title sprintf("%s : %s : %s", drivername, seriesname, probname);
    set grid;

    set xlabel "time";
    set xrange [0:1500];
    set xtics 100;

    set multiplot layout 3,1;

    if (0) {
        set ylabel "energy / mass";
        set yrange [0:20];
        plot \
            enefn u 1:($4/$10) w l lw 2 lc 2 title "KE",\
            enefn u 1:($14/$10) w l ls 0 lc 2 lw 1 title "ThKE",\
            enefn u 1:($5/$10) w l lw 2 lc 3 title "RKE",\
            enefn u 1:($6/$10) w l lw 2 lc 4 title "EPE";
    }

    set ylabel "xmom";
    set yrange [0:50];
    plot enefn u 1:11 w l lw 2 title "xmom";

    set ylabel "depth";
    set yrange [0:6];
    plot enefn u 1:(2*$8) w l lw 2 title "depth";

    set ylabel "heat content";
    set yrange [0:1];
    plot enefn u 1:($14/$4) w l lw 2 title "ThKE / KE";

    unset multiplot;

}
