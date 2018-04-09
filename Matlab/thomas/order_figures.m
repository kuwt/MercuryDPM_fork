function order_figures(n)

hfig=get(0,'Children');
ScreenSize=get(0,'ScreenSize');
numx=max(3,floor(length(hfig)^(ScreenSize(3)/(ScreenSize(3)+ScreenSize(4)))));
if exist('n','var'), numx=n; end
numy=max(2,ceil(length(hfig)/numx));
Position=get(gcf,'Position');
dx=(ScreenSize(3)-Position(3))/max(1,numx-1);
dy=(ScreenSize(4)-Position(4))/max(1,numy-1);

for i=1:length(hfig)
  figure(hfig(i));
  Position=get(hfig(i),'Position');
  %disp([hfig(i) mod((i-1),numx)*dx floor((i-1)/numx)*dy   Position(3:4)]);
  set(hfig(i),'Position',[mod((i-1),numx)*dx floor((i-1)/numx)*dy   Position(3:4)]);
end

return