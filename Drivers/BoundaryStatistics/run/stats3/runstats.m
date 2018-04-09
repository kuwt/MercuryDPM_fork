if true
  addpath('../../matlab')
  addpath('~/code/MD/matlab/thomas')
  dataSteadyL0=loadstatistics('H10A22L0M0.5B0.5W0.2Stress1.stat');
  dataSteadyL2=loadstatistics('H10A26L2M0.5B0.5W0.2Stress1.stat');
end

data={dataSteadyL2};
color=lines(7);

figure(1); clf; hold on
for i=1:length(data)
  plot(data{i}.z,deriv(data{i}.StressZZ,data{i}.z),'Color',color(i,:),'DisplayName',data{i}.name)
  plot(data{i}.z,data{i}.Gravity(3)*data{i}.Density,'k:')
end
ylabel('(\nabla\sigma)_{z}')
legend('show')
axis tight; %ylim([min(ylim) 0])

figure(2); clf; hold on
for i=1:length(data)
  plot(data{i}.z,deriv(data{i}.StressZZ,data{i}.z)+data{i}.TractionZ,'Color',color(i,:),'DisplayName',data{i}.name)
  plot(data{i}.z,data{i}.Gravity(3)*data{i}.Density,'k:')
end
ylabel('(\nabla\sigma)_{z}+t_z')
legend('show')
axis tight;

figure(3); clf; hold on
for i=1:length(data)
  plot(data{i}.z,data{i}.TractionZ,'Color',color(i,:),'DisplayName',data{i}.name)
end
ylabel('t_z')
legend('show')
axis tight;

figure(4); clf; hold on
for i=1:length(data)
  plot(data{i}.z,data{i}.StressZZ,'Color',color(i,:),'DisplayName',data{i}.name)
  plot(data{i}.z,data{i}.Gravity(3)*(cumsum(data{i}.Density)-sum(data{i}.Density))*diff(data{i}.z(1:2)),'k:')
%   sum(data{i}.StressZZ)/sum(-data{i}.Gravity(3)*data{i}.Density)*2
  [sum(data{i}.Density.*data{i}.StressZZ)*diff(data{i}.z(1:2)) -data{i}.Gravity(3)/2*(sum(data{i}.Density)*diff(data{i}.z(1:2))).^2 ]
end
ylabel('(\nabla\sigma)_{z}')
legend('show')
axis tight; %ylim([min(ylim) 0])

