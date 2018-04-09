function h=plotPrincipalStress(data,h,color,legend_)
% this function creates a plot of the energy ratio ene_kin/ene_ela in the
% chute, including a legend. You should be able to specify figure, suplot, 
% hold outside of the function.

%create new axis if none is given
if ~exist('h','var')||isempty(h)
  axis;
  h=gca;
end

%create new legend if none is given
if ~exist('legend_','var')||isempty(legend_)
  legend_=['$h=' num2str(data.FlowHeight,2) ', a=' num2str(data.ChuteAngle,2) '$'];
end

%if you want principle stresses
data=loadPrincipalStress(data);

for i=1:3
  %plot
  subplot(3,1,i)
  %indexCenter = data.z<.66*data.FlowHeight&data.z>.33*data.FlowHeight;
  %h=plot(data.z(indexCenter),data.PrincipalAngle(i,indexCenter));
  h=plot(data.z,data.PrincipalAngle(i,:));
  StatisticsLineStyle(h)
  %label axes
  xlabel('$z$')
  ylabel(['$\phi_' num2str(i) '$'])
  axis tight
end

%add to legend
hleg = legend; legend([get(hleg,'String'), legend_],'Location','South');

return
