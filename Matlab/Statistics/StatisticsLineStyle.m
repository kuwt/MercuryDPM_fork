function StatisticsLineStyle(h)
%create new color if none is given
if ~exist('color','var')||isempty(color)
  C=length(get(gca,'Children'));
  CO=get(gca,'ColorOrder');
  CC=size(CO,1);
  LSO=get(gca,'LineStyleOrder');
  set(h,'Color',CO(mod(C-1,CC)+1,:));
  set(h,'LineStyle',LSO(mod(ceil(C/CC)-1,length(LSO))+1,:));
  set(h,'LineWidth',2);
end
return