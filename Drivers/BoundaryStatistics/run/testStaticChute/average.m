%%
addpath('../../matlab')

data=loadstatistics('3/static2d.stat');

figure(3);clf;
set(gcf,'FileName','static2d')

subplot(4,2,1)
plot(mean(data.y,2),mean(data.Density,2));
title('$\rho$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,3)
plot(mean(data.y,2),mean(data.StressXX,2));
title('$\sigma_{xx}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,4)
plot(mean(data.y,2),mean(data.StressXY,2));
title('$\sigma_{xy}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,5)
plot(mean(data.y,2),mean(data.StressYX,2));
title('$\sigma_{yx}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,6)
plot(mean(data.y,2),mean(data.StressYY,2));
title('$\sigma_{yy}$','Interpreter','tex')
%xlabel('$y$') 

data.NablaStressX=deriv(data.StressXX,data.x)+deriv(data.StressXY,data.y);
data.NablaStressY=deriv(data.StressYX,data.x)+deriv(data.StressYY,data.y);

subplot(4,2,7)
plot(mean(data.y,2),mean(data.Density*data.Gravity(1)-data.NablaStressX,2));
title('$(\rho g-\nabla\sigma)_x$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,8)
plot(mean(data.y,2),mean(data.Density*data.Gravity(2)-data.NablaStressY,2));
title('$(\rho g-\nabla\sigma)_y$','Interpreter','tex')
%xlabel('$y$') 



%%

data=loadstatistics('0/static2d.stat');

figure(4);clf;
% set(gcf,'Position',[830 570 22 15])
set(gcf,'FileName','static2d')

subplot(4,2,1)
plot(mean(data.y,2),mean(data.Density,2));
title('$\rho$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,3)
plot(mean(data.y,2),mean(data.StressXX,2));
title('$\sigma_{xx}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,4)
IntTractionX=cumsum(mean(data.TractionX,2))*diff(data.y(1:2));
IntTractionX=IntTractionX-IntTractionX(end);
plot(mean(data.y,2),mean(data.StressXY,2)+IntTractionX);
title('$\sigma_{xy}+\int t_x$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,5)
plot(mean(data.y,2),mean(data.StressYX,2));
title('$\sigma_{yx}$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,6)
IntTractionY=cumsum(mean(data.TractionY,2))*diff(data.y(1:2));
IntTractionY=IntTractionY-IntTractionY(end);
plot(mean(data.y,2),[mean(data.StressYY,2)+IntTractionY]);
title('$\sigma_{yy}+\int t_y$','Interpreter','tex')
%xlabel('$y$') 

data.NablaStressX=deriv(data.StressXX,data.x)+deriv(data.StressXY,data.y);
data.NablaStressY=deriv(data.StressYX,data.x)+deriv(data.StressYY,data.y);

subplot(4,2,7)
plot(mean(data.y,2),mean(data.Density*data.Gravity(1)-data.TractionX-data.NablaStressX));
title('$(\rho g-t-\nabla\sigma)_x$','Interpreter','tex')
%xlabel('$y$') 

subplot(4,2,8)
plot(mean(data.y,2),mean(data.Density*data.Gravity(2)-data.TractionY-data.NablaStressY));
title('$(\rho g-t-\nabla\sigma)_y$','Interpreter','tex')
%xlabel('$y$') 

% print_figures();
