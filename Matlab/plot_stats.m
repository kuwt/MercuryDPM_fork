function plot_stats
close all;

data=load_stats('../DRIVERS/Simple_MD/run/free_cooling/free_cooling.stat');

lin_engy_fig=figure;
plot(data(:,1),data(:,2));
rot_engy_fig=figure;
plot(data(:,1),data(:,3))
x_COM_fig=figure;
plot(data(:,1),data(:,4))
Y_COM_fig=figure;
plot(data(:,1),data(:,5))
sym=figure 
plot(data(:,4),data(:,5));

data=load_stats('../DRIVERS/Simple_MD/run/free_cooling/free_cooling_hgrid.stat');

figure(lin_engy_fig)
hold on;
plot(data(:,1),data(:,2),'r');
figure(rot_engy_fig);
hold on;
plot(data(:,1),data(:,3),'r');
figure(x_COM_fig);
hold on;
plot(data(:,1),data(:,4),'r')
figure(Y_COM_fig);
hold on;
plot(data(:,1),data(:,5),'r')

figure(sym)
hold on
plot(data(:,4),data(:,5),'r');

