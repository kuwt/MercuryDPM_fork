clear variables

%% parameters
R = 50e-6;
rho = 1050;
L=56.4e3;
cSolid  =1200;
cLiquid = 1.2*cSolid;
Tm = 451.15;
dT = 20;
tMax = 1.1; % from code
iota=0.1;
f = 1-iota^3;
cMelt  = (cSolid+cLiquid)/2+L/dT;
m = rho*4/3*pi*R^3;
H = (2*dT*m*(cLiquid+cSolid)/2 + L*m)/tMax;
dir = '/Users/weinhartt/Code/melting/cmake-build-debug/Drivers/Papers/MeltingModel/';
n=1;
name='OneParticleCoolingSelfTest';

%% read sim data
data = importdata([dir name '.out']);
data = data(1:n:end,:);
sim.t = data(:,1);
data(sim.t>tMax,:)=[];
sim.t = data(:,1);
sim.T = data(:,2);
sim.X = data(:,3)/R;
sim.EMelt = data(:,4);
sim.EHeat = data(:,5);

%% plot temperature
figure(1) 
set(gcf,'Position',[0 600 600 180]) 
subplot(1,3,1)
T0 = Tm+dT;
t = linspace(0,tMax,1000);
opt = odeset("RelTol",1e-15);
[~,T]=ode45(@(t,T) dTCooling(t,T,cSolid,cLiquid,dT,R,Tm,m,L), t, T0, opt);
T=T';
tL=interp1(T,t,Tm-dT/2);
tS=interp1(T,t,Tm+dT/2);
plot(t,T,'DisplayName','ode45')
hold on
plot(sim.t,sim.T,'.','DisplayName','simulation')
plot(t,(Tm-dT/2)*ones(size(t)),'k:','HandleVisibility','off')
plot(t,(Tm+dT/2)*ones(size(t)),'k:','HandleVisibility','off')
axis tight
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthEast')
xlabel('t [s]')
ylabel('T_i [K]')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
yticks([Tm-dT/2 Tm+dT/2])
%yticklabels({'T_m-\Delta{T}/2','T_m+\Delta{T}/2'})
set(gca,'Box','off')
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);

%% plot melt thickness
subplot(1,3,2)
c = (cSolid+cLiquid)/2+L/dT;
X=1-(1-(T-Tm+dT/2)/dT*f).^(1/3);
X(t<tS)=1-iota;
X(t>tL)=0;
plot(t,X,'DisplayName','analytic')
hold on
plot(sim.t,sim.X,'.','DisplayName','simulation')
axis tight
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthEast')
xlabel('t [s]')
ylabel('\chi_i/R_i')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
yticks([0 1-iota])
%yticklabels({'0','\iota'})
set(gca,'Box','off')
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);

%% plot energies
subplot(1,3,3)
EInput = H*t;
T0=Tm-dT;
EHeat =          cSolid*m*(min(T,Tm-dT/2)-T0) ...
     +(cLiquid+cSolid)/2*m*(max(min(T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
     +           cLiquid*m*(max(T,Tm+dT/2)-(Tm+dT/2)); 
Rs = R-X*R;
EMelt = L*(m-rho*4/3*pi*Rs.^3)/f;
EMax = H*tMax;
sim.EMelt = L*(m-rho*4/3*pi*(R-sim.X*R).^3);
sim.EHeat =           cSolid*m*(min(sim.T,Tm-dT/2)-T0) ...
    +(cLiquid+cSolid)/2*m*(max(min(sim.T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
    +           cLiquid*m*(max(sim.T,Tm+dT/2)-(Tm+dT/2)); 
%plot(t,EInput,'DisplayName','$E^{\rm in}$') %=\int H dt
plot(t,1e6*EMelt,'DisplayName','$E^{\rm m} (ana.)$') %+=L_m(m_0-(4/3)\pi\rho {R_0^s}^3)
hold on 
plot(t,1e6*EHeat,'DisplayName','$E^{\rm h} (ana.)$') %=\int c^p\dot{T} dt
plot(sim.t,1e6*sim.EMelt,'.','Color',[0 0.4470 0.7410],'DisplayName','$E^{\rm m} (sim.)$')
plot(sim.t,1e6*sim.EHeat,'.','Color',[0.8500 0.3250 0.0980],'DisplayName','$E^{\rm h} (sim.)$')
axis tight
%ylim([0 EMax])
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
%plot(xlim,1e3*(H*tMax-L*m)*[1 1],'k:','HandleVisibility','off')
plot(xlim,1e6*L*m*[1 1],'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthEast','Interpreter','latex')
xlabel('t [s]')
ylabel('Energy [uJ]')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
%yticks([0 L*m H*tMax])
%yticklabels({'0','L\cdot{}m_i','H_i\cdot{}t^{max}'})
set(gca,'Box','off')
pos = get(gca, 'Position'); set(gca, 'Position', [pos(1) pos(2)+0.06 pos(3) pos(4)-0.06]);

saveas(gcf,[name '.png'])
system(['/opt/homebrew/bin/convert ' name '.png -trim ' name '.png'])

