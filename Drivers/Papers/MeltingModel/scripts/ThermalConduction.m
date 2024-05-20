clear all

%% parameters
R = 50e-6;
rho = 1050;
L=56.4e3;
cS = 1200;
cL = 1.2*cS;
Tm = 451.15;
dT = 20;
m = rho*4/3*pi*R^3;
dir = '/Users/weinhartt/Code/melting/cmake-build-debug/Drivers/Papers/MeltingModel/';
name='ThermalConductionSelfTest';
k = 0.12;
overlap = 0.1*R;
s = pi*R*overlap;
dist = 2*R-overlap;

%% read sim data
data = importdata([dir name '.out']);
%data = data(1:savecount:end,:);
sim.t = data(:,1);
sim.T0 = data(:,2);
sim.T1 = data(:,4);
T0=sim.T0(1);
T1=sim.T1(1);
tMax = max(sim.t);

%% plot temperature
figure(1)
%DTdt = -2*k*s/dist/m/cS*DT 
K = 2*k*s/dist/m/cS;
T0_ = 0.5*(T0+T1)+0.5*(T0-T1)*exp(-K*sim.t);
T1_ = 0.5*(T0+T1)-0.5*(T0-T1)*exp(-K*sim.t);
clf
set(gcf,'Position',[00 600 200 180]) 
plot(sim.t,sim.T0,'.','DisplayName','T_i')
hold on
plot(sim.t,sim.T1,'.','DisplayName','T_j')
h(1)=plot(sim.t,T0_,'Color',[.5 .5 .5],'DisplayName','analytic')
h(2)=plot(sim.t,T1_,'Color',[.5 .5 .5],'HandleVisibility','off')
uistack(h,'up',1);
hold off
axis tight
legend('show')
set(legend,'Location','SouthEast')
xlabel('t [s]')
ylabel('T [K]')
set(gca,'Box','off')
saveas(gcf,[name '.png'])
system(['/opt/homebrew/bin/convert ' name '.png -trim ' name '.png'])




