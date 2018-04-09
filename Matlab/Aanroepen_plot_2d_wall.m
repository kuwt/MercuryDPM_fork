clear all
close all
clc

figure
hold on;
s1='numFiniteWalls 3 normal -0.34202014 0 -0.93969262 position -0.04341525 normal -0.98480775 0 -0.17364818 position -0.025522814 normal 0.87215428 0 0.48923094 position 0.022603225';
s2='numFiniteWalls 3 normal 0.64278761 0 -0.76604444 position -0.015392436 normal 0.98480775 0 0.17364818 position 0.030522814 normal -0.97368212 0 0.22791034 position -0.02624273';
s3='numFiniteWalls 3 normal -1 0 -0 position -0.005 normal 0 0 1 position 0 normal 0.81811048 0 -0.57506107 position 0';




% s1='numFiniteWalls 3 normal -0.34202014 0 -0.93969262 position -0.037600226 normal -0.98480775 0 -0.17364818 position -0.024448241 normal 0.68477196 0 0.72875741 position 0.029160007';
% s2='numFiniteWalls 3 normal 0.64278761 0 -0.76604444 position -0.010651985 normal 0.98480775 0 0.17364818 position 0.029448241 normal -0.75737073 0 0.65298513 position 0.00082910992';
% s3='numFiniteWalls 3 normal -1 0 -0 position -0.1 normal 0 0 1 position 0 normal 0.8660254 0 -0.5 position 0';

plot_2d_wall(s1);
plot_2d_wall(s2);
plot_2d_wall(s3);
axis equal