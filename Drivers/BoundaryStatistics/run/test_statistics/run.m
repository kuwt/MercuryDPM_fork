addpath('../../matlab')
allData=loadstatistics('*/static2d.stat');

for StatType=0:3;
data=allData{StatType+1};

figure(StatType+1);clf;
set(gcf,'Position',[1117-600+StatType*200 30 560 917])
set(gcf,'FileName',['static2d' num2str(StatType,'_%d')])

subplot(6,2,2)
[C,h]=contourf(data.x,data.y,data.Density,100,'LineStyle','none');
colorbar
axis equal
title('$\rho$')

subplot(6,2,3)
contourf(data.x,data.y,data.StressXX,100,'LineStyle','none')
colorbar
axis equal
title('$\sigma_{xx}$')

subplot(6,2,4)
contourf(data.x,data.y,data.StressXY,100,'LineStyle','none')
colorbar
axis equal
title('$\sigma_{xy}$')

subplot(6,2,5)
contourf(data.x,data.y,data.StressYX,100,'LineStyle','none')
colorbar
axis equal
title('$\sigma_{yx}$')

subplot(6,2,6)
contourf(data.x,data.y,data.StressYY,100,'LineStyle','none')
colorbar
axis equal
title('$\sigma_{yy}$')

subplot(6,2,7)
contourf(data.x,data.y,data.TractionX,100,'LineStyle','none')
colorbar
axis equal
title('$t_x$')

subplot(6,2,8)
contourf(data.x,data.y,data.TractionY,100,'LineStyle','none')
colorbar
axis equal
title('$t_y$')

subplot(6,2,9)
data.NablaStressX=deriv(data.StressXX,data.x)+deriv(data.StressXY,data.y);
contourf(data.x,data.y,data.NablaStressX,100,'LineStyle','none')
colorbar
axis equal
title('$(\nabla \sigma)_x$')

subplot(6,2,10)
data.NablaStressY=deriv(data.StressYX,data.x)+deriv(data.StressYY,data.y);
contourf(data.x,data.y,data.NablaStressY,100,'LineStyle','none')
colorbar
axis equal
title('$(\nabla \sigma)_y$')

subplot(6,2,11)
contourf(data.x,data.y,data.Density*data.Gravity(1)-data.TractionX-data.NablaStressX,100,'LineStyle','none')
colorbar
axis equal
title('$(\rho g - \nabla \sigma)_x$')

subplot(6,2,12)
contourf(data.x,data.y,data.Density*data.Gravity(2)-data.TractionY-data.NablaStressY,100,'LineStyle','none')
colorbar
axis equal
title('$(\rho g - \nabla \sigma)_y$')

end

%print_figures();
