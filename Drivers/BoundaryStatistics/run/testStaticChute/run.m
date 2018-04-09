addpath('../../matlab');
allData=loadstatistics('*/H*.stat');
figure(1);clf;
set(gcf,'Position',[561 30 1116 917])
set(gcf,'FileName',['static2d' num2str(0:3,'_%d')])
x=[-2 5];
% allData{4}.TractionX=zeros(size(allData{3}.TractionZ));
% allData{4}.TractionZ=zeros(size(allData{3}.TractionZ));

for StatType=0:3;
data=allData{StatType+1};
data.NablaStressX=deriv(data.StressXZ,data.z);
data.NablaStressZ=deriv(data.StressZZ,data.z);

subplot(4,3,1+3*StatType)
% plot(data.z,[data.StressXX data.StressXZ data.StressZX data.StressZZ]);
plot(data.z,[data.Density*data.Gravity(1) data.NablaStressX +data.TractionX])
legend('show')
title('$(\rho g, \nabla \sigma+t)_x$')
xlim(x)

subplot(4,3,2+3*StatType)
plot(data.z,[data.Density*data.Gravity(3) data.NablaStressZ +data.TractionZ])
title('$(\rho g, \nabla \sigma+t)_y$')
xlim(x)

%plot remainder
subplot(4,3,3+3*StatType)
plot(data.z,[...
  data.Density*data.Gravity(1)-data.TractionX-data.NablaStressX,...
  data.Density*data.Gravity(3)-data.TractionZ-data.NablaStressZ])
title('$(\rho g - t - \nabla \sigma)$')
xlim(x)

end
% print_figures();
