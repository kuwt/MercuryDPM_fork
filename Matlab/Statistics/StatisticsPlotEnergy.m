function h=StatisticsPlotEnergy(data,h,color,legend_)
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
  legend_=['$h=' num2str(data.FlowHeight,3) ', \theta=' num2str(data.ChuteAngle,2) '^\circ$'];
end


%plot
ind=data.Ene.Time>1;
z=data.Ene.Kin(ind)./data.Ene.Ela(ind);
% Speed=data.Speed(data.Speed<100);
% z=data.ParticleDensity*data.ParticleVolume*Speed*Speed'./data.Ene.Ela(ind);
h=loglog(h,data.Ene.Time(ind),(z+z([1 1:end-1])+z([1 1 1:end-2])+z([2:end end])+z([3:end end end]))/5);

%obtain standard color/linestyle scheme
StatisticsLineStyle(h)
set(h,'LineWidth',1)

%label axes
xlabel('$t$')
ylabel('$Ene_{kin}/Ene_{ela}$')

%resize axis
axis tight

%add to legend
hleg = legend; legend([get(hleg,'String'), legend_],'Location','SouthWest');

return
