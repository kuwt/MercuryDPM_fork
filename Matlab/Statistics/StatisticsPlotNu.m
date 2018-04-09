function h=StatisticsPlotNu(data,h,color,legend_)
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


%plot
h=plot(h,data.z,data.Nu);

StatisticsLineStyle(h)

%label axes
xlabel('$z$')
ylabel('$\nu$')
axis tight

%add to legend
hleg = legend; legend([get(hleg,'String'), legend_],'Location','South');

return
