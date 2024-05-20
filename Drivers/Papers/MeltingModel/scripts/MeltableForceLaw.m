clear variables

%% parameters
dir = '/Users/weinhartt/Code/melting/cmake-build-debug/Drivers/Papers/MeltingModel/';
name='MeltableForceLawSelfTest';
R = 50e-6;
rho = 1050;
L=56.4e3;
cSolid  =1200;
cLiquid = 1.2*cSolid;
Tm = 451.15;
dT = 20;
tMax = 0.0120299;
cMelt  = (cSolid+cLiquid)/2+L/dT;
m = rho*4/3*pi*R^3;
H = (2*dT*m*(cLiquid+cSolid)/2 + L*m)/tMax;
deltab = 0.55*0.1*R;

%% read sim data
data = importdata([dir name '_0_1.out']);
%data = data(1:savecount:end,:);
sim0.t = data(:,1);
sim0.delta = data(:,3);
sim0.force = -data(:,4);

% scaleFactors = 10.^[4 6 8];
% scaleStrings={'10^4','10^6','10^8'};
scaleFactors = [1 1e8];
scaleStrings={'1','10^8'};
sim=struct;
for i =1:length(scaleFactors)
    filename = [dir name '_0.1_' num2str(scaleFactors(i),"%g") '.out']
    data = importdata(filename);
    %data = data(1:savecount:end,:);
    sim(i).t = data(:,1);
    sim(i).delta = data(:,3);
    sim(i).force = -data(:,4);
end

%% plot
figure(1) 
cla
set(gcf,'Position',[600 600 220 200]) 
plot(sim0.delta/R,sim0.force,'-','DisplayName','\zeta=0')
axis tight
%xlim([0,max(xlim)])
hold on
plot(sim(i).delta/R,sim(i).force,'-','DisplayName','\zeta=0.1, S=1')
for i=2:length(scaleFactors)
    n=30;
    plot(sim(i).delta(1:n:end)/R,sim(i).force(1:n:end),'.','DisplayName',['\zeta=0.1, S=' scaleStrings{i}])
end
%plot(1e3*t,T,'DisplayName','analytic')
plot(deltab/R*[1 1],ylim,'k:','HandleVisibility','off')
plot(xlim,[0 0],'k:','HandleVisibility','off')
hold off
legend('show')
set(legend,'Location','NorthWest')
xlabel('\delta_{ij}^s/R_{ij}^s')
ylabel('f [N]')
set(gca,'Box','off')
saveas(gcf,[name '.png'])
system(['/opt/homebrew/bin/convert ' name '.png -trim ' name '.png'])




