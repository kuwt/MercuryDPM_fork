clear all

%% parameters
R = 50e-6;
rho = 1050;
L = 56.4e3;
cS  =1200;
cL = 1.2*cS;
Tm = 451.15;
dT = 20;
H = 0.06072; %from code
m = rho*4/3*pi*R^3;
cMelt  = (cS+cL)/2+L/dT;
tMax = (dT*(cS+cL)+L)*m/H;
iota=0.9;
f=1/(1-(1-iota)^3);
dir = '/Users/weinhartt/Code/melting/cmake-build-debug/Drivers/Papers/MeltingModel/';
name='OneParticleHeatingSelfTest';

%% read sim data
data = importdata([dir name '.out']);
sim.t = data(:,1);
sim.T = data(:,2);
sim.X = data(:,3)/R;
sim.EMelt = data(:,4);
sim.EHeat = data(:,5);

%% plot temperature
figure(1) 
set(gcf,'Position',[0 600 600 180]) 
subplot(1,3,1)
T0 = Tm-dT;
tS = (Tm-dT/2-T0)/(H/m/cS) %time when particle starts to melt
tL = tS + dT/(H/m/cMelt) %time when particle melting ends 
t = linspace(0,tMax,1000);
T(t<tS)=T0+(H/m/cS)*t(t<tS);
T(t>tS)=Tm-dT/2+(H/m/cMelt)*(t(t>tS)-tS);
T(t>tL)=Tm+dT/2+(H/m/cL)*(t(t>tL)-tL);
plot(1e3*t,T,'DisplayName','analytic')
hold on
plot(1e3*sim.t,sim.T,'.','DisplayName','simulation')
plot(1e3*t,(Tm-dT/2)*ones(size(t)),'k:','HandleVisibility','off')
plot(1e3*t,(Tm+dT/2)*ones(size(t)),'k:','HandleVisibility','off')
axis tight
plot(1e3*tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(1e3*tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest')
xlabel('t [ms]')
ylabel('T_i [K]')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
yticks([Tm-dT/2 Tm+dT/2])
%yticklabels({'T_m-\Delta{T}/2','T_m+\Delta{T}/2'})
set(gca,'Box','off')
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);

%% plot melt thickness
subplot(1,3,2)
X(t>tS)=1-(1-H*(t(t>tS)/f-tS)/m/(cMelt*dT)).^(1/3);
X(t>tL)=iota;
plot(1e3*t,X,'DisplayName','analytic')
hold on
plot(1e3*sim.t,sim.X,'.','DisplayName','simulation')
axis tight
plot(1e3*tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(1e3*tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest')
xlabel('t [ms]')
ylabel('\chi_i/R_i')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
yticks([0 iota])
set(gca,'Box','off')
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);

%% plot energies
subplot(1,3,3)
EInput = H*t;
T0=Tm-dT;
EHeat =          cS*m*(min(T,Tm-dT/2)-T0) ...
     +(cL+cS)/2*m*(max(min(T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
     +           cL*m*(max(T,Tm+dT/2)-(Tm+dT/2)); 
Rs = R-X*R;
EMelt = L*(m-rho*4/3*pi*Rs.^3)/f;
EMax = H*tMax;
sim.EMelt = L*(m-rho*4/3*pi*(R-sim.X*R).^3);
sim.EHeat =           cS*m*(min(sim.T,Tm-dT/2)-T0) ...
    +(cL+cS)/2*m*(max(min(sim.T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
    +           cL*m*(max(sim.T,Tm+dT/2)-(Tm+dT/2)); 
%plot(t,EInput,'DisplayName','$E^{\rm in}$') %=\int H dt
plot(t,1e6*EMelt,'DisplayName','$E^{\rm m} (ana.)$') %+=L_m(m_0-(4/3)\pi\rho {R_0^s}^3)
hold on 
plot(1e3*t,1e6*EHeat,'DisplayName','$E^{\rm h} (ana.)$') %=\int c^p\dot{T} dt
plot(1e3*sim.t,1e6*sim.EMelt,'.','Color',[0 0.4470 0.7410],'DisplayName','$E^{\rm m} (sim.)$')
plot(1e3*sim.t,1e6*sim.EHeat,'.','Color',[0.8500 0.3250 0.0980],'DisplayName','$E^{\rm h} (sim.)$')
axis tight
%ylim([0 EMax])
plot(1e3*tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(1e3*tL*[1 1],ylim,'k:','HandleVisibility','off')
%plot(xlim,1e3*(H*tMax-L*m)*[1 1],'k:','HandleVisibility','off')
plot(xlim,1e6*L*m*[1 1],'k:','HandleVisibility','off')
plot(xlim,1e6*(H*tMax-L*m)*[1 1],'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest','Interpreter','latex')
xlabel('t [ms]')
ylabel('Energy [uJ]')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
%yticks(1e3*[0 H*tMax-L*m L*m])
%yticklabels({'0','E^H','E^M'})
set(gca,'Box','off')
axis tight
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);
saveas(gcf,[name '.png'])
system(['/opt/homebrew/bin/convert ' name '.png -trim ' name '.png'])



