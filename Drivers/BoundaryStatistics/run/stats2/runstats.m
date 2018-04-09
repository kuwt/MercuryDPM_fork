close all
if ~exist('dataStaticL2','var')
  addpath('../../matlab')
  addpath('~/code/MD/matlab/thomas')
%   dataStaticL0=loadstatistics('H20A12L0M0.5B0.5W0.1*.stat');
%   dataStaticL2=loadstatistics('H10A20L2M0.5B0.5W0.1*.stat');
  dataSteadyL0=loadstatistics('../stats3/H10A22L0M0.5B0.5W0.2Stress*.stat');
  dataSteadyL2=loadstatistics('../stats3/H10A26L2M0.5B0.5W0.2Stress*.stat');
end
color=lines(7);
P=[0 0 560 .6*420];

figure(1); clf; hold on
set(gcf,'Position',P)
set(gcf,'FileName','L2Stress')

data=dataSteadyL2;
%i=1;
%plot(data{i}.z,data{i}.StressZZ,'k:','DisplayName','$\sigma_{zz}^0$')
i=2;
plot(data{i}.z,data{i}.StressZZ,'k--','DisplayName','$\sigma_{zz}~~$')
%i=3;
% plot(data{i}.z,data{i}.StressZZ,'k:','DisplayName','$\sigma_{zz}$')
i=4;
plot(data{i}.z,data{i}.StressZZ,'k','DisplayName','$\sigma_{zz}''$')
plot(data{i}.z,data{i}.Gravity(3)*(cumsum(data{i}.Density)-sum(data{i}.Density))*diff(data{i}.z(1:2)),'r:','LineWidth',2,'DisplayName','\rg')
xlabel('z')
legend('show')
%set(legend,'Box','off')
axis tight; 
plot(data{i}.Base*[1 1],ylim,'Color',[1 1 1]*.8,'LineWidth',.7)
plot(data{i}.Surface*[1 1],ylim,'Color',[1 1 1]*.8,'LineWidth',.7)
xlim([data{i}.Base data{i}.Surface]+0.1*data{i}.FlowHeight*[-1 1])

% figure(2); clf; hold on
% set(gcf,'Position',P)
% set(gcf,'FileName','L0Stress')
% 
% data=dataSteadyL0;
% i=2;
% plot(data{i}.z,data{i}.StressZZ,'k--','DisplayName','$\sigma_{zz}$')
% i=4;
% plot(data{i}.z,data{i}.StressZZ,'k','DisplayName','$\sigma_{zz}+\int_z^\infty t_z dz$')
% plot(data{i}.z,data{i}.Gravity(3)*(cumsum(data{i}.Density)-sum(data{i}.Density))*diff(data{i}.z(1:2)),'r:','LineWidth',2,'DisplayName','$\int_z^\infty g_z \rho dz$')
% xlabel('z')
% legend('show')
% axis tight; 
% plot(data{i}.Base*[1 1],ylim,'k:')
% plot(data{i}.Surface*[1 1],ylim,'k:')
% xlim([data{i}.Base data{i}.Surface]+0.2*data{i}.FlowHeight*[-1 1])
% 
% 
% figure(3); clf; hold on
% set(gcf,'Position',P)
% set(gcf,'FileName','L2Traction')
% 
% data=dataSteadyL2;
% i=4;
% nablaS=sqrt(...
%   (deriv(data{i}.StressXZ,data{i}.z)+data{i}.TractionX).^2+...
%   (deriv(data{i}.StressYZ,data{i}.z)+data{i}.TractionY).^2+...
%   (deriv(data{i}.StressZZ,data{i}.z)+data{i}.TractionZ).^2);
% plot(data{i}.z,nablaS,'k-','DisplayName','$|\nabla\V{\sigma}_{zz}+\V{t}|$')
% plot(data{i}.z,norm(data{i}.Gravity)*data{i}.Density,'r:','LineWidth',2,'DisplayName','$g \rho$')
% xlabel('z')
% legend('show')
% set(legend,'Location','South','Box','off')
% axis tight; 
% plot(data{i}.Base*[1 1],ylim,'k:')
% plot(data{i}.Surface*[1 1],ylim,'k:')
% xlim([data{i}.Base data{i}.Surface]+0.2*data{i}.FlowHeight*[-1 1])
% 
% 
% figure(4); clf; hold on
% set(gcf,'Position',P)
% set(gcf,'FileName','L0Traction')
% 
% data=dataSteadyL0;
% i=4;
% nablaS=sqrt(...
%   (deriv(data{i}.StressXZ,data{i}.z)+data{i}.TractionX).^2+...
%   (deriv(data{i}.StressYZ,data{i}.z)+data{i}.TractionY).^2+...
%   (deriv(data{i}.StressZZ,data{i}.z)+data{i}.TractionZ).^2);
% plot(data{i}.z,nablaS,'k-','DisplayName','$|\nabla\V{\sigma}_{zz}+\V{t}|$')
% plot(data{i}.z,norm(data{i}.Gravity)*data{i}.Density,'r:','LineWidth',2,'DisplayName','$g \rho$')
% xlabel('z')
% legend('show')
% set(legend,'Location','South','Box','off')
% axis tight; 
% plot(data{i}.Base*[1 1],ylim,'k:')
% plot(data{i}.Surface*[1 1],ylim,'k:')
% xlim([data{i}.Base data{i}.Surface]+0.2*data{i}.FlowHeight*[-1 1])


% figure(5); clf; hold on
% set(gcf,'Position',P)
% set(gcf,'FileName','L2StressXX')
% 
% data=dataSteadyL2;
% i=1;
% plot(data{i}.z,data{i}.StressXX,'k:','DisplayName','$\sigma_{xx}^0$')
% i=2;
% plot(data{i}.z,data{i}.StressXX,'k--','DisplayName','$\sigma_{xx}$')
% i=4;
% plot(data{i}.z,data{i}.StressXX,'k','DisplayName','$\sigma_{xx}`$')
% xlabel('z')
% legend('show')
% axis tight; 
% plot(data{i}.Base*[1 1],ylim,'k:')
% plot(data{i}.Surface*[1 1],ylim,'k:')
% xlim([data{i}.Base data{i}.Surface]+0.2*data{i}.FlowHeight*[-1 1])



print_figures();
