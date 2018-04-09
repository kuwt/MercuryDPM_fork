function add_triangle_to_loglog(slope,lx,ly,reverse,scale)
hold on

%size of triangle   
s=diff(log10(xlim))/3;
if exist('scale','var')
   s=s*scale; 
end
%distance from corner
d=diff(log10(xlim))/10;

%
x=10^(log10(min(xlim))+lx*diff(log10(xlim)));
y=10^(log10(min(ylim))+ly*diff(log10(ylim)));
if exist('reverse','var')&&reverse==true
    h=plot(x*10.^(d+s*[0 1 1 0]),y*10.^(d+s*slope*[0 0 -1 0]),'k-');
    text(x*10^(d+s),y*10^((d-s*slope)/2),num2str(abs(slope)),'HorizontalAlignment','Right');
    text(x*10^((d+s)/2),y*10^(d*.9),num2str(1),'VerticalAlignment','top');
else
    h=plot(x*10.^(d+s*[0 1 0 0]),y*10.^(d+s*slope*[1 0 0 1]),'k-');
    text(x*10^(d),y*10^((d+s*slope)/2),num2str(abs(slope)),'HorizontalAlignment','Right');
    text(x*10^((d+s)/2),y*10^(d*.9),num2str(1),'VerticalAlignment','top');
end
set(h,'HandleVisibility','off')
return