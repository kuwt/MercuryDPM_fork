close all

addpath('~/code/MD/matlab/thomas')
dataXZ=loadstatistics('static2d.tra.stat');
dataZ=loadstatistics('static2d.z.stat');

dataP=importdata('static2d.data',' ',1);
dataP=dataP.data(end-9:end,1:3);

figure(1);clf;
set(gcf,'Position',[0 0 560 420])
set(gcf,'FileName','2DTraction')

numcolors =15;
contourf(dataXZ.x,dataXZ.z,sqrt(dataXZ.TractionX.^2+dataXZ.TractionY.^2+dataXZ.TractionZ.^2),numcolors,'LineStyle','none');
%contourf(dataXZ.x,dataXZ.z,sqrt(dataXZ.TractionZ.^2),numcolors,'LineStyle','none');
invgray=bone(numcolors); invgray=invgray(end:-1:1,:);
colormap(invgray)
colorbar
axis image
title('$|\V{t}|_2$')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$')
ylabel('$z$');

hold on;
x=(0:.05:2*pi)';
% for i=1:5
%   patch(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),[0 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.2)
% end
% for i=6:10
%   patch(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),[0 0 0],'FaceAlpha',0.0,'EdgeAlpha',0.2)
% end
for i=1:10
  plot(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),'Color',.8*[1 1 1])
end
for i=1:5
  plot(dataP(i,1),dataP(i,3),'x','Color',.8*[1 1 1])
end



figure(2);clf;
set(gcf,'Position',[0 0 560 420])
set(gcf,'FileName','2DStress')

StressNorm=zeros(size(dataXZ.StressXX));
for i=1:length(dataXZ.StressXX(:))
   StressNorm(i)=norm([dataXZ.StressXX(i) dataXZ.StressXY(i) dataXZ.StressXZ(i);dataXZ.StressYX(i) dataXZ.StressYY(i) dataXZ.StressYZ(i);dataXZ.StressZX(i) dataXZ.StressZY(i) dataXZ.StressZZ(i)]);
%  StressNorm(i)=norm([dataXZ.StressZZ(i)]);
end
contourf(dataXZ.x,dataXZ.z,StressNorm,numcolors,'LineStyle','none')
colormap(invgray)
colorbar
axis image
title('$|\V{\sigma}|_2$')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$')
ylabel('$z$')

hold on;
x=0:.05:2*pi;
% for i=1:5
%   patch(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),[0 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.2)
% end
% for i=6:10
%   patch(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),[0 0 0],'FaceAlpha',0.0,'EdgeAlpha',0.2)
% end
for i=1:10
  plot(sin(x)/2+dataP(i,1),cos(x)/2+dataP(i,3),'Color',.8*[1 1 1])
end
for i=1:5
  plot(dataP(i,1),dataP(i,3),'x','Color',.8*[1 1 1])
end


figure(3)
set(gcf,'Position',[0 0 .8*560 420])
set(gcf,'FileName','2D')
rgb=imread('../snapshots/2D.bmp');
rgb(rgb>0)=255;
image(rgb)

colorbar
P=get(gca,'Position');
colorbar off
set(gca,'Position',P);

axis image
title(' ')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('$x$')
ylabel('$z$')

addpath('../../matlab')
print_figures();
