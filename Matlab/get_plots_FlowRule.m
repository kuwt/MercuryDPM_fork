function data=get_plots_FlowRule(data)

close all

if (~exist('data'))
	data=loadstatistics({ ...
		'silbert.theta.22.z.5.stat' ...
		'silbert.theta.26.z.5.stat' ...
		'silbert.theta.30.z.5.stat' ...
		'silbert.theta.22.z.10.stat' ...
		'silbert.theta.26.z.10.stat' ...
		'silbert.theta.30.z.10.stat' ...
		'silbert.theta.22.z.20.stat' ...
		'silbert.theta.26.z.20.stat' ...
		'silbert.theta.30.z.20.stat' ...
		'silbert.theta.22.z.40.stat' ...
		'silbert.theta.26.z.40.stat' ...
		'silbert.theta.30.z.40.stat' ...
		});
	data=loadstatistics({ ...
		'silbert.5.stat' ...
		'silbert.6.stat' ...
		'silbert.7.stat' ...
		});

end

for i=1:length(data)
  data{i}.name = 'thetah';
  %data{i}.label = '\theta=14^\circ H=10';
  data{i}.PlotVariable = data{i}.z;
  data{i}.PlotVariableLabel = 'z';
end

color = {'k','b','r','g','c','m','y'};
ptype = {'-','--',':','-.'};
HList = {'$H=5d$','$H=10d$','$H=20d$','$H=40d$',};
thetaList = {'$\theta=22^\circ$','$\theta=26^\circ$','$\theta=30^\circ$'};

%plot pvf for H=20,
h=figure(); clf;
set(gcf,'PaperSize',[4.8 3.5],'Position',[0 0 400 350]); 
set(gcf,'paperpositionmode','auto'); 
set(gca,'PlotBoxAspectRatio',[1 .7 1]);
H=2;
for i=3*H+(1:3)
  hold on
  plot(data{i}.PlotVariable,data{i}.Nu,color{i-3*H+1},'LineWidth',2);
end
xlabel(data{1}.PlotVariableLabel);
v=axis; axis([v(1:3) .64]);
legend(thetaList,3);
% title(['Particle Volume fraction for ' HList{H}]);
filename = [data{1}.name 'Nu'];
%PrintLaTeX(gcf,filename,8); 

%plot velocity for H=20
h=figure(); clf;
set(gcf,'PaperSize',[4.8 3.5],'Position',[0 0 400 350]); 
set(gcf,'paperpositionmode','auto'); 
set(gca,'PlotBoxAspectRatio',[1 .7 1]);
for H=3
for i=3*H+(1:3)
  hold on
	set(gco,'Interpreter','none')
	plot(data{i}.PlotVariable/min(data{i}.FlowHeight,40),data{i}.VelocityX/data{i}.FlowVelocityX,[color{i-3*H+1} '--'],'LineWidth',2);%,ptype{H+1}
	M = 5/3*sum(data{i}.VelocityX(~isnan(data{i}.VelocityX)))/length(data{i}.VelocityX(~isnan(data{i}.VelocityX)));
end
end
z=0:.05:1;
u=(5/3)*(1-(1-z).^1.5);
plot(z,u,'k','LineWidth',2);
v = axis; axis([0 1 0 1.75]);
legend({thetaList{1:end},'Bagnold'},4,'Interpreter','none');
% title(['Velocity Profile for ' HList{3+1}],'Interpreter','none');
filename = [data{1}.name 'VelocityX'];
xlabel('$z/h$','Interpreter','none');
ylabel('$u/\bar{u}$','Interpreter','none');
%PrintLaTeX(gcf,filename,8); 

%plot StressZZ for theta=26
h=figure(); clf;
set(gcf,'PaperPosition',1.2*[0 0 4 3.5],'Position',[0 0 400 350]); 
set(gca,'PlotBoxAspectRatio',[1 .7 1]);
for th=1:2
for i=3*(0:3)+th
  hold on
	dz = data{i}.z(2)-data{i}.z(1);
	angle=str2num(data{i}.label(15:16));
	Gravity = cos(angle/180*pi);
	%plot(data{i}.PlotVariable/min(data{i}.FlowHeight,40),data{i}.VelocityX/d
	%ata{i}.FlowVelocityX,[color{i-3*H+1}],'LineWidth',2);%,ptype{H+1}
	mg = data{i}.ParticleDensity*Gravity*sum(data{i}.Nu)*dz;
  plot(data{i}.PlotVariable/data{i}.FlowHeight,data{i}.StressZZ/mg,[color{(i-th)/3+2} ptype{th}],'LineWidth',2);
end
end
xlabel('$z/h$','Interpreter','none');
ylabel('$\sigma_{zz}/(\int\!\!\rho g\, dz)$','Interpreter','none');
axis([0 1.05 0 1])
legend(HList,1,'Interpreter','none');
% title(['$\sigma_{zz}$ for ' thetaList{1} '(solid),' thetaList{2} '(dotted)'],'Interpreter','none');
filename = [data{1}.name 'StressZZ'];
%PrintLaTeX(gcf,filename,8); 

%plot K for theta=26 with Bagnold
h=figure(); clf;
set(gcf,'PaperSize',[4.8 3.5],'Position',[0 0 400 350]); 
set(gcf,'paperpositionmode','auto'); 
set(gca,'PlotBoxAspectRatio',[1 .7 1]);
for th=1:2
for H=2
	i=3*H+th;
	plot(data{i}.PlotVariable/data{i}.FlowHeight, ...
		[data{i}.StressXX./data{i}.Pressure;data{i}.StressYY./data{i}.Pressure;data{i}.StressZZ./data{i}.Pressure], ...
		[ptype{th}],'LineWidth',2);
	hold on
end
end
xlabel('$z/h$','Interpreter','none');
ylabel('$\sigma_{ii}/p$','Interpreter','none');
legend({'$\sigma_{xx}/p$','$\sigma_{yy}/p$','$\sigma_{zz}/p$'},-1,'Location','South','Interpreter','none');
axis tight;
% title(['Pressure components for ' thetaList{th-1} '(-),'  thetaList{th} '(:),' HList{H+1}],'Interpreter','none');
filename = [data{1}.name 'PressureUniformity'];
%PrintLaTeX(gcf,filename,8); 
 
 

%plot mu
h=figure(); clf;
set(gcf,'PaperPosition',[0 0 4 3],'Position',[0 0 400 300]);
set(gcf,'paperpositionmode','auto'); 
for th=1
for i=3*(0:3)+th
	hold on
	plot(data{i}.PlotVariable/data{i}.FlowHeight,data{i}.MacroFrictionCoefficient,color{(i-th)/3+2},'LineWidth',2);
end
end
plot([0 1],tan((18+4*th)/180*pi)*[1 1],'k','LineWidth',2);
xlabel('$z/h$','Interpreter','none');
ylabel('$\mu_{Macro}$','Interpreter','none');
v = axis; axis([0 1 0 .5])
List = HList;
List{end+1} = '$tan(\theta)$';
if th==1, legend(List,3,'Interpreter','none'); end
% title(['Macroscopic friction coefficient for ' thetaList{th}],'Interpreter','none');
filename = [data{1}.name 'Mu'];
%PrintLaTeX(gcf,filename,8); 

%plot coordination number
h=figure(); clf;
set(gcf,'PaperPosition',[0 0 4 3],'Position',[0 0 400 300]);
set(gcf,'paperpositionmode','auto'); 
for th=1:3
	for H=0:3
		i=3*H+th;
		FlowHeight(H+1) = data{i}.FlowHeight;
		CoordinationNumber(H+1) = data{i}.TotalCoordinationNumber;
	end
	hold on
	plot(FlowHeight,CoordinationNumber,[color{th} 'o:']);
end
xlabel('$h$','Interpreter','none');
ylabel('$Z$','Interpreter','none');
legend(thetaList,9,'Interpreter','none');
% title(['Coordination number for ' thetaList{th}],'Interpreter','none');
filename = [data{1}.name 'Z'];
%PrintLaTeX(gcf,filename,8); 

%plot Fabric anisotropy
h=figure(); clf;
set(gcf,'PaperPosition',[0 0 4 3],'Position',[0 0 400 300]);
set(gcf,'paperpositionmode','auto'); 
for th=1:3
	for H=0:3
		i=3*H+th;
		FlowHeight(H+1) = data{i}.FlowHeight;
		FabricAnisotropy(H+1) = abs(data{i}.TotalFabricXZ);
	end
	hold on
	plot(FlowHeight,FabricAnisotropy,[color{th} 'o:']);
end
xlabel('$h$','Interpreter','none');
ylabel('$|F_{xz}|$','Interpreter','none');
legend(thetaList,9,'Interpreter','none');
% title(['Fabric anisotropy for ' thetaList{th}],'Interpreter','none');
filename = [data{1}.name 'FabricAnisotropy'];
%PrintLaTeX(gcf,filename,8); 

%plot pressure anisotropy
h=figure(); clf;
set(gcf,'PaperPosition',[0 0 4 5],'Position',[0 0 400 500]);
set(gcf,'paperpositionmode','auto');
ylabel_ = {'xx','yy','zz'};
for th=1:3
	for H=0:3
		i=3*H+th;
		FlowHeight(H+1) = data{i}.FlowHeight;
		ind=~isnan(data{i}.Pressure);
		NormalStressOverPressure(:,H+1) = [ ...
			sum(data{i}.StressXX(ind)) sum(data{i}.StressYY(ind)) sum(data{i}.StressZZ(ind))]/ sum(data{i}.Pressure(ind));
	end
	for j=1:3
		subplot(3,1,j);
		hold on
		plot(FlowHeight,NormalStressOverPressure(j,:),[color{th} 'o:']);
		ylabel(['$|\sigma_{' ylabel_{j} '}|/P$'],'Interpreter','none');
	end
end
xlabel('$h$','Interpreter','none');
legend(thetaList,2,'Interpreter','none');
% title(['Pressure anisotropy for ' thetaList{th}],'Interpreter','none');
filename = [data{1}.name 'PressureAnisotropy'];
%PrintLaTeX(gcf,filename,8); 

return