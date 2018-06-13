%%
addpath('../../matlab')


for StatType=0:3;

figure(StatType+1);clf;
set(gcf,'Position',[561+200*StatType 30 1116-600 917])
set(gcf,'FileName',['static2d' num2str(0:3,'_%d')])

data=loadstatistics([num2str(StatType) '/static2d.stat']);
dataY=loadstatistics([num2str(StatType) '/static2dY.stat']);

subplot(4,2,1)
plot(mean(data.y,2),mean(data.TractionX,2));
hold on; plot(dataY.y,dataY.TractionX,'r--');
title('$t_x$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,2)
plot(mean(data.y,2),mean(data.TractionY,2));
hold on; plot(dataY.y,dataY.TractionY,'r--');
title('$t_y$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,3)
plot(mean(data.y,2),mean(data.StressXX,2));
hold on; plot(dataY.y,dataY.StressXX,'r--');
title('$\sigma_{xx}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,4)
plot(mean(data.y,2),mean(data.Density,2));
hold on; plot(dataY.y,dataY.Density,'r--');
title('$\rho$','Interpreter','tex')
%xlabel('$y$') 
% 
% subplot(4,2,5)
% plot(mean(data.y,2),mean(data.StressYX,2));
% hold on; plot(dataY.y,dataY.StressYX,'r--');
% title('$\sigma_{yx}$','Interpreter','tex')
% %xlabel('$y$') 
% 
% subplot(4,2,6)
% plot(mean(data.y,2),mean(data.StressYY,2));
% hold on; plot(dataY.y,dataY.StressYY,'r--');
% title('$\sigma_{yy}$','Interpreter','tex')
% %xlabel('$y$') 

data.NablaStressX=deriv(data.StressXX,data.x)+deriv(data.StressXY,data.y);
data.NablaStressY=deriv(data.StressYX,data.x)+deriv(data.StressYY,data.y);
dataY.NablaStressX=deriv(dataY.StressXY,dataY.y);
dataY.NablaStressY=deriv(dataY.StressYY,dataY.y);

subplot(4,2,5)
plot(mean(data.y,2),mean(data.NablaStressX,2));
hold on; plot(dataY.y,dataY.NablaStressX,'r--');
title('$\nabla\sigma_{x}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,6)
plot(mean(data.y,2),mean(data.NablaStressY,2));
hold on; plot(dataY.y,dataY.NablaStressY,'r--');
title('$\nabla\sigma_{y}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,7)
plot(mean(data.y,2),mean(data.Density*data.Gravity(1)-data.NablaStressX-data.TractionX,2));
hold on; plot(dataY.y,dataY.Density*dataY.Gravity(1)-dataY.NablaStressX-dataY.TractionX,'r--');
title('$(\rho g-\nabla\sigma-t)_x$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,8)
plot(mean(data.y,2),mean(data.Density*data.Gravity(2)-data.NablaStressY-data.TractionY,2));
hold on; plot(dataY.y,dataY.Density*dataY.Gravity(2)-dataY.NablaStressY-dataY.TractionY,'r--');
title('$(\rho g-\nabla\sigma-t)_y$','Interpreter','tex')
%xlabel('$y$') 

end
