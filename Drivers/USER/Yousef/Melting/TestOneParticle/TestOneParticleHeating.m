clear all

%% parameters
R = 115e-6;
rho = 1050;
L=56.4e3;
cSolid  =1200;
cLiquid = 1.2*cSolid;
Tm = 451.15;
dT = 20;
tMax = 1;

cMelt  = (cSolid+cLiquid)/2+L/dT;
m = rho*4/3*pi*R^3;
H = (2*dT*m*(cLiquid+cSolid)/2 + L*m)/tMax;

%% read sim data
data = importdata('OneParticleHeatingSelfTest.out');
sim.t = data(:,1);
sim.T = data(:,2);
sim.X = data(:,3)/R;
sim.EMelt = data(:,4);
sim.EHeat = data(:,5);

%% plot temperature
figure(1) 
set(gcf,'Position',[600 600 1200 400]) 
subplot(1,3,1)
T0 = Tm-dT;
tS = (Tm-dT/2-T0)/(H/m/cSolid) %time when particle starts to melt
tL = tS + dT/(H/m/cMelt) %time when particle melting ends 
t = linspace(0,tMax,1000);
T(t<tS)=T0+(H/m/cSolid)*t(t<tS);
T(t>tS)=Tm-dT/2+(H/m/cMelt)*(t(t>tS)-tS);
T(t>tL)=Tm+dT/2+(H/m/cLiquid)*(t(t>tL)-tL);
plot(t,T,'DisplayName','analytic')
hold on
plot(sim.t,sim.T,'.','DisplayName','simulation')
plot(t,(Tm-dT/2)*ones(size(t)),'k:','HandleVisibility','off')
plot(t,(Tm+dT/2)*ones(size(t)),'k:','HandleVisibility','off')
axis tight
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest')
xlabel('t')
ylabel('T_0')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
% yticks([Tm-dT/2 Tm+dT/2])
% yticklabels({'T_m-\Delta T/2','T_m+\Delta T/2'})
set(gca,'Box','off')

%% plot melt thickness
subplot(1,3,2)
X(t>tS)=1-(1-H*(t(t>tS)-tS)/m/(cMelt*dT)).^(1/3);
X(t>tL)=1;
plot(t,X,'DisplayName','analytic')
hold on
plot(sim.t,sim.X,'.','DisplayName','simulation')
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest')
xlabel('t')
ylabel('\chi/R_0')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
% yticks([Tm-dT/2 Tm+dT/2])
% yticklabels({'T_m-\Delta T/2','T_m+\Delta T/2'})
set(gca,'Box','off')

%% plot energies
subplot(1,3,3)
EInput = H*t;
EHeat =          cSolid*m*(min(T,Tm-dT/2)-T0) ...
    +(cLiquid+cSolid)/2*m*(max(min(T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
    +           cLiquid*m*(max(T,Tm+dT/2)-(Tm+dT/2)); 
Rs = R-X*R;
EMelt = L*(m-rho*4/3*pi*Rs.^3);
EMax = H*tMax;
sim.EMelt = L*(m-rho*4/3*pi*(R-sim.X*R).^3);
sim.EHeat =           cSolid*m*(min(sim.T,Tm-dT/2)-T0) ...
    +(cLiquid+cSolid)/2*m*(max(min(sim.T,Tm+dT/2),Tm-dT/2)-(Tm-dT/2)) ...
    +           cLiquid*m*(max(sim.T,Tm+dT/2)-(Tm+dT/2)); 
%plot(t,EInput,'DisplayName','$E^{\rm in}$') %=\int H dt
plot(t,EMelt,'DisplayName','$E^{\rm m} (analytic)$') %+=L_m(m_0-(4/3)\pi\rho {R_0^s}^3)
hold on 
plot(t,EHeat+EMelt,'DisplayName','$E^{\rm m}+E^{\rm h} (analytic)$') %=\int c^p\dot{T} dt
plot(sim.t,sim.EMelt,'.','Color',[0 0.4470 0.7410],'DisplayName','$E^{\rm m} (simulation)$')
plot(sim.t,sim.EHeat+sim.EMelt,'.','Color',[0.8500 0.3250 0.0980],'DisplayName','$E^{\rm m}+E^{\rm h} (simulation)$')
ylim([0 EMax])
plot(tS*[1 1],ylim,'k:','HandleVisibility','off')
plot(tL*[1 1],ylim,'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest','Interpreter','latex')
xlabel('t')
ylabel('Energy')
% xticks([tS tL])
% xticklabels({'t_s','t_l'})
% yticks([0 L*m H*tMax])
% yticklabels({'0','L\cdot{}m','H\cdot{}t^{max}'})
set(gca,'Box','off')
saveas(gcf,'Heating.png')



