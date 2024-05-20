clear variables

%% parameters
dir = '/Users/weinhartt/Code/melting/cmake-build-debug/Drivers/Papers/MeltingModel/';
name='TestCoalescence';
R = 50e-6;
rho = 1050;
L=56.4e3;
cSolid  =1200;
cLiquid = 1.2*cSolid;
Tm = 451.15;
dT = 20;
tMax = 0.23;
cMelt  = (cSolid+cLiquid)/2+L/dT;
m = rho*4/3*pi*R^3;
savecount = 20;
%deltab = 0.55*0.1*R;
surfaceTension = 34.4e-3;
viscosity = 48;

%% read sim data
scaleFactors = 10.^[0 6 8];
scaleStrings={'1','10^6','10^8'};
%scaleFactors = 10.^[8];
%scaleStrings={'10^8'};
sim=struct;
for i =1:length(scaleFactors)
    data = importdata([dir name '_' num2str(scaleFactors(i),"%g") '.out']);
    data = data(1:savecount:end,:);
    %data = data(data(:,1)<tMax,:);
    sim(i).t = data(:,1);
    sim(i).moltenLayerThickness0 = data(:,2);
    sim(i).moltenLayerThickness1 = data(:,3);
    sim(i).delta = data(:,4);
    sim(i).force = -data(:,5);
    sim(i).neckRadius = min(2^(1/3)*R,sqrt(2.0*sim(i).delta*R));
end
%% solve
odefun = @(t,delta) 0.83/2*surfaceTension*(2*R-delta)/viscosity/R;
opt = odeset("RelTol",1e-13);
t = linspace(0,max(sim(1).t));
[~,delta]=ode45(odefun, t, 0, opt);
delta=delta';
delta = min(delta,1.8*R);
neckRadius = min(2^(1/3)*R,sqrt(2.0*delta*R));
%% plot
figure(1) 
set(gcf,'Position',[000 600 400 180]) 
%% plot overlap
subplot(1,2,1)
cla
%plot(sim.delta/R,sim.force,'-','DisplayName','\zeta=0')
plot(1e3*t,delta/R,'-','Linewidth',1,'DisplayName','ana.')
%plot(sim.t,sim.moltenLayerThickness0/R,'-','DisplayName','\zeta=0')
hold on
for i=1:length(scaleStrings)
    %scaleStrings{i}
    plot(1e3*sim(i).t,sim(i).delta/R,'.','DisplayName',['S=' scaleStrings{i}])
end
%plot(1e3*t,T,'DisplayName','analytic')
axis tight
% plot(deltab/R*[1 1],ylim,'k:','HandleVisibility','off')
% plot(xlim,(sim.moltenLayerThickness0+sim.moltenLayerThickness1)/R*[1 1],'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','SouthEast')
ylabel('\delta_{ij}/R_{ij}')
xlabel('t [ms]')
set(gca,'Box','off')

%% plot neck radius 
subplot(1,2,2)
% cla
plot(1e3*t,neckRadius/R,'-','Linewidth',1,'DisplayName','ana.')
%plot(sim.delta/R,sim.force,'-','DisplayName','\zeta=0')
%plot(sim.t,sim.moltenLayerThickness0/R,'-','DisplayName','\zeta=0')
hold on
for i=1:length(scaleStrings)
    %scaleStrings{i}
    %plot(1e3*sim(i).t,sim(i).delta/R,'-','DisplayName',['S=' scaleStrings{i}])
    plot(1e3*sim(i).t,sim(i).neckRadius/R,'.','DisplayName',['S=' scaleStrings{i}])
end

%plot(1e3*t,T,'DisplayName','analytic')
axis tight
% plot(deltab/R*[1 1],ylim,'k:','HandleVisibility','off')
%plot(xlim,(sim.moltenLayerThickness0+sim.moltenLayerThickness1)/R*[1 1],'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','SouthEast')
ylabel('a_{ij}/R_{ij}')
xlabel('t [ms]')
set(gca,'Box','off')
% %% plot force
% subplot(1,3,3)
% cla
% %plot(sim.delta/R,sim.force,'-','DisplayName','\zeta=0')
% plot(1e3*sim.t,1e6*sim.force,'-','DisplayName','\zeta=0')
% %plot(sim.t,sim.moltenLayerThickness0/R,'-','DisplayName','\zeta=0')
% hold on
% %plot(sim.t,sim.moltenLayerThickness1/R,'-','DisplayName','\zeta=0')
% %plot(1e3*t,T,'DisplayName','analytic')
% axis tight
% %ylim(1e6*[min(sim.force),-min(sim.force)])
% % plot(deltab/R*[1 1],ylim,'k:','HandleVisibility','off')
% % plot(xlim,[0 0],'k:','HandleVisibility','off')
% hold off
% legend('sim.')
% set(legend,'Location','NorthWest')
% ylabel('f_{ij} [\mu{N}]')
% xlabel('t [ms]')
% set(gca,'Box','off')
exportgraphics(gcf,[name '.png'],'Resolution',1000);
system(['/opt/homebrew/bin/convert ' name '.png -trim ' name '.png']);




