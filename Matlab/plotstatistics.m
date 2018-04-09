function [data,info]=plotstatistics(filenames)


if (exist('filenames'))
    if (iscell(filenames))
      [data,info]=plot_statistics(filenames)
    else
        [data,info]=plot_statistics({filenames})
    end
else
if strcmp(pwd,'/storage/usr/people/weinhartt/thomas/MDCLRexamples/chute_MDCLR')
  [data, info] = plot_statistics({'chute_periodic.6.tav9.data';'tav9.data'});
elseif strcmp(pwd,'/storage/usr/people/weinhartt/code/MD/DRIVERS/Chute/run/silbert_autorun')
  [data, info] = plot_statistics({ ...
    'silbert.theta.14.z.5.stat' ...
    'silbert.theta.18.z.5.stat' ...
    'silbert.theta.22.z.5.stat' ...
    'silbert.theta.26.z.5.stat' ...
    'silbert.theta.30.z.5.stat' ...
    });
elseif strcmp(pwd,'/storage/usr/people/weinhartt/code/MD/DRIVERS/Chute/run/silbert')
  [data, info] = plot_statistics({ ...
    'silbert.13.stat', ...
    'silbert.13.2.stat', ...
    'tav9.data' ...
    });
elseif strcmp(pwd,'/storage/usr/people/weinhartt/code/MD/DRIVERS/Chute/run/chute_periodic')
  [data, info] = plot_statistics({'chute_periodic.6.stat','tav9.data'});
else
  [data, info] = plot_statistics({'segregation.1.large'});
end
end
return





function [data,info]=plot_statistics(filenames)
%[data] = make_readable(data, info);
for i = 1:length(filenames)
  %load file
  [data, info] = load_file(filenames{i});
  [data, info] = get_standard_variables(data,info);
  %[data, info] = get_I(data,info);
  %[data, info] = extract_variables(data,info,[27:32]);
  %[data, info] = get_momentum_equation(data,info);
  %[data, info] = get_energy_equation(data,info);
  
  %description{end+1} = [ label{mod(i,14)+1} ': ' filenames{i} ', ' info.description ', '];


figure(1);
if strcmp(pwd,'/storage/usr/people/weinhartt/code/MD/DRIVERS/Chute/run/silbert_autorun') 
  % Temperature on logplog
  % VelocityX on logplog
  for i = [2 5]
    subplot(4,4,i);
    set(gca,'YScale','log');
  end
end
 subplot(1,2,1);
 v = axis;
 axis([v(1:2) 0 .6]);
 subplot(1,2,2);
 v = axis;
 axis([v(1:2) .3 .6]);

htextbox = annotation('textbox',[0 .9 1 .1],'EdgeColor','none');
%set(htextbox,'String',description,'Interpreter','none')
%['Labels: black: ' description1 ', red: ' description2]

%saveas(gcf ,['figure1.pdf']);

return


end




function [data, info] = load_data_file(filename)
%load stat file data

disp(['loading data from ' filename]);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);
  
% reads variable names from line 1 and text output
% from line 2
info.variable_names = { ...
  'Nu'; 'Density'; ...
  'MomentumX'; 'MomentumY'; 'MomentumZ'; ...
  'MomentumFluxXX'; 'MomentumFluxXY'; 'MomentumFluxXZ'; 'MomentumFluxYY'; 'MomentumFluxYZ'; 'MomentumFluxZZ'; ...
  'EnergyFluxX'; 'EnergyFluxY'; 'EnergyFluxZ'; ...
  'NormalStressXX'; 'NormalStressXY'; 'NormalStressXZ'; 'NormalStressYY'; 'NormalStressYZ'; 'NormalStressZZ'; ...
  'TangentialStressXX'; 'TangentialStressXY'; 'TangentialStressXZ'; 'TangentialStressYY'; 'TangentialStressYZ'; 'TangentialStressZZ'; ...
  'FabricXX'; 'FabricXY'; 'FabricXZ'; 'FabricYY'; 'FabricYZ'; 'FabricZZ'; ...
  'CollisionalHeatFluxX'; 'CollisionalHeatFluxY'; 'CollisionalHeatFluxZ'; ...
  'Dissipation'; ...
  };
text = textscan([rawdata.textdata{8} ' ' rawdata.textdata{9}],'%s ');
info.text = [text{1}{5} num2str(str2double(text{1}{9})*100)];
rho_text = textscan(rawdata.textdata{11},'%s ');
rho = str2double(rho_text{1}{end});

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data.time = 0;
data.coordinates = rawdata.data(:,4:6);
nu = max(rawdata.data(:,10),1e-200*ones(size(rawdata.data,1),1));
data.variables = [nu ...
  nu*rho ...
  rawdata.data(:,[51:53]).*nu(:,[1 1 1])*rho ...
  rawdata.data(:,[41:43 45 46 49]) ...
  nan(size(rawdata.data,1),3) ...
  rawdata.data(:,[21:23 25 26 29]) ...
  rawdata.data(:,[31:33 35 36 39]) ...
  rawdata.data(:,[71:73 75 76 79]) ...
  nan(size(rawdata.data,1),4) ...
  ];

info.description = ['Luding ' info.text];

return



function [data, info] = get_energy_equation(data,info)

if isfield(data,'variance'), data = rmfield(data,'variance'); end

info.variable_names = {'PVF'; ...
  'VelocityX'; 'VelocityY'; 'VelocityZ'; ...
  'StressXX'; 'StressXY'; 'StressXZ'; 'StressYY'; 'StressYZ'; 'StressZZ'; ...
  'Temperature'; ...
  'Pressure'; 'Dissipation'; ...
  'HeatFluxX'; 'HeatFluxY'; 'HeatFluxZ';};
Nu = data.variables(:,1);
Velocity = data.variables(:,[3:5])./data.variables(:,[2 2 2]);
ContactStress = data.variables(:,[15:20])+data.variables(:,[21:26]);
Temperature = sum(data.variables(:,[6 9 11])./data.variables(:,[2 2 2])-Velocity.^2 ,2)/3;
Pressure = sum(ContactStress(:,[1 4 6]),2)/3;
Potential = data.variables(:,37);
Dissipation = data.variables(:,36);

Density = data.variables(:,2);
EnergyFlux = data.variables(:,[12:14]);
MomentumFlux = data.variables(:,[6:11]);
CollisionalHeatFlux = data.variables(:,[33:35]) - 0.5*( ...
    ContactStress(:,[1 2 3]).*Velocity(:,[1 1 1]) ...
  + ContactStress(:,[2 4 5]).*Velocity(:,[2 2 2]) ...
  + ContactStress(:,[3 5 6]).*Velocity(:,[3 3 3]) ...
  + Potential(:,[1 1 1]).*Velocity);
KineticHeatFlux = (EnergyFlux ...
  - .5*(sum(MomentumFlux(:,[1 4 6]),2)*ones(1,3)).*Velocity ...
  - MomentumFlux(:,[1 2 3]).*Velocity(:,[1 1 1]) ...
  - MomentumFlux(:,[2 4 5]).*Velocity(:,[2 2 2]) ...
  - MomentumFlux(:,[3 5 6]).*Velocity(:,[3 3 3]) ...
  + Density(:,[1 1 1]).*(sum(Velocity.*Velocity,2)*ones(1,3)).*Velocity) ./Density(:,[1 1 1]);
% EnergyFlux_b-MomentumFlux_aa*Velocity_b-2*MomentumFlux_ab*Velocity_a+2*Velocity_a*Velocity_a*Velocity_b
data.variables = [ Nu Velocity ContactStress Temperature Pressure Dissipation CollisionalHeatFlux ];

return

function [data, info] = get_momentum_equation(data,info)

[data, info] = get_standard_variables(data,info)

NablaVelocity = nabla(data.variables(:,2:4),data.coordinates);
VelocityDotNablaVelocity = data.variables(:,[2 2 2]).*NablaVelocity{1} ...
                         + data.variables(:,[3 3 3]).*NablaVelocity{2} ...
                         + data.variables(:,[4 4 4]).*NablaVelocity{3} ;

NablaStress = nabla(data.variables(:,5:10),data.coordinates);
NablaDotStress = NablaStress{1}(:,[1 2 3]) ...
               + NablaStress{2}(:,[2 4 5]) ...
               + NablaStress{3}(:,[3 5 6]) ;

density = 1.9;
%Gravity = [0.42261826 0 -0.90630779];
Gravity = [0 0 -1];
DensityNuGravity = density * data.variables(:,1) * Gravity;
DensityNuVelocityDotNablaVelocity = density * data.variables(:,[1 1 1]) .* VelocityDotNablaVelocity;
Remainder = DensityNuVelocityDotNablaVelocity -DensityNuGravity +NablaDotStress;

data.variables = [DensityNuVelocityDotNablaVelocity DensityNuGravity NablaDotStress Remainder];

%sum_of_all = density_nu_VelocityZ_dzVelocity-density_nu_Gravity+NablaStress;
%data.variables = [density_nu_VelocityZ_dzVelocity density_nu_Gravity NablaStress sum_of_all];
%index = data.coordinates(:,3)>.005 & data.coordinates(:,3)<.015;
d = 1;
index = data.coordinates(:,3)>5*d;

data.variables = data.variables(index,:);
data.coordinates = data.coordinates(index,:);

info.variable_names = {
  'density_nu_VelocityZ_dzVelocityX';
  'density_nu_VelocityZ_dzVelocityY';
  'density_nu_VelocityZ_dzVelocityZ';
  'density_nu_GravityX';
  'density_nu_GravityY';
  'density_nu_GravityZ';
  'NablaStressX';
  'NablaStressY';
  'NablaStressZ';
  'SumX';
  'SumY';
  'SumZ';
  };

return


function [data, info] = get_I(data,info)

[data, info] = get_standard_variables(data,info)
strainrate = ...
  (data.variables([2:end end],2)-data.variables([1 1:end-1],2)) ./ ...
  (data.coordinates([2:end end],3)-data.coordinates([1 1:end-1],3));
  %[.1; .2*ones(198,1); .1];
d = 1;
rho = 1.9;
info.variable_names = {'I';'deviator stress ratio sD1'};
data.variables = [strainrate./sqrt(data.variables(:,12)./data.variables(:,1)/rho) data.variables(:,end-3)];

return

function [data, info] = extract_variables(data,info,indexes)

info.variable_names = info.variable_names(indexes);
data.variables = data.variables(:,indexes); 

if isfield(data,'variance'), 
  data.variance = data.variance(:,indexes); 
end

return




%set h
if (isfield(info,'h')), h = info.h; else, h=1; end
disp(['plotting into figure' num2str(h) ' ...']);

if (~isfield(info,'plottype')), info.plottype = '-'; end

% checks which dimensions are really used and stores their coordinates
dim = ~[ ...
  all(data.coordinates(1,1)==data.coordinates(:,1)) ...
  all(data.coordinates(1,2)==data.coordinates(:,2)) ...
  all(data.coordinates(1,3)==data.coordinates(:,3)) ];
dim_label_full = 'xyz';
dim_label = dim_label_full(dim);
coordinates = data.coordinates(:,dim);
disp(['dim=' num2str(dim) ' ...']);

% stores the grid dimensions
for i=1:sum(dim)
  n(i) = size(coordinates,1) / sum(coordinates(1,i)==coordinates(:,i));
end

% stores coordinates in grid
for i=1:sum(dim)
  if (length(n)==1)
    X{i} = coordinates(:,i);
  else
    X{i} = zeros(n(end:-1:1));
    X{i}(:) = coordinates(:,i);
  end
end
Z = zeros(size(X{1}));
if (variance), Zvar = zeros(size(X{1})); end

%create figure
figure(h);
%clear figure if info.hold does not exist
if (~isfield(info,'hold')), clf; end
%resize figure     
set(gcf,'Position',[300 1 1379 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
%plot info text
%info_text = [info.text ' t=' mat2str(data.time)];
%global htextbox
%if ~isempty('htextbox'), htextbox = annotation('textbox',[0 0 1 1]); end
%set(htextbox,'String',info_text ,'Interpreter','none')

% defines number of variables and plot order
numvar = size(data.variables,2);
plot_i = 1:numvar;
if (numvar==1), 
  ynum = 1;
  xnum = 1;
elseif (numvar==2), 
  ynum = 1;
  xnum = 2;
elseif (numvar<=4), 
  ynum = 2;
  xnum = 2;
% elseif (numvar==6), %plotting a symmetric matrix
%   ynum = 3;
%   xnum = 3;
%   plot_i = [1:3 5 6 9];
elseif (numvar<=9), 
  ynum = 3;
  xnum = 3;
elseif (numvar<=12), 
  ynum = 4;
  xnum = 3;
elseif (numvar<=16), 
  ynum = 4;
  xnum = 4;
  plot_i = [1 5 9 13 6 7 8 11 12 16 2 3 4 10 14 15];
elseif (numvar<=24)
  ynum = 4;
  xnum = 6;
  plot_i = [1:2 4:numvar+1];
elseif (numvar<=36)
  ynum = 6;
  xnum = 6;
%  plot_i = [1:2 4:15 19:36 16:18 3];
elseif (numvar<=37)
  ynum = 6;
  xnum = 6;
  plot_i = [1 1:numvar];
end

if (sum(dim)==1)
  % for 1D output
  %set(htextbox,'String',['plotting 1-dimensional data' info_text],'Interpreter','none')
  
	for i=1:numvar
    subplot(ynum,xnum,plot_i(i)); 
		hold on
    Z(:) = data.variables(:,i);
    plot(X{1},Z,info.plottype);
    plot(sum(X{1})/length(X{1}),sum(Z)/length(Z),'rx');
    if (variance), 
      Zvar(:) = sqrt(data.variance(:,i)); 
      plot(X{1},Z-Zvar,':');
      plot(X{1},Z+Zvar,':');
    end
    title(info.variable_names{i},'Interpreter','none');
    xlabel(dim_label);
		hold off
	end

elseif (sum(dim)==2)
  % for 2D output
  
	if (~(isfield(info,'dim_avg') && dim(info.dim_avg)) )
		%plotting without averaging
		%set(htextbox,'String',['plotting 2-dimensional data' info_text],'Interpreter','none')
		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); hold on
			Z(:) = data.variables(:,i);
			binplot(X{1},X{2},Z);
			colorbar;
			shading flat;
			axis equal;
			title(info.variable_names{i},'Interpreter','none');
			xlabel(dim_label(1));
			ylabel(dim_label(2));
		end
	else
		%plotting averaged data
 
		dim_avg = info.dim_avg(1);
    dim_avg_bool = zeros(1,3);
    dim_avg_bool(dim_avg) = true;
		dim_avg_local = find(find(dim)==dim_avg);
 		%set(htextbox,'String',['plotting ' dim_label_full(dim_avg) '-averaged 2-dimensional data' info_text],'Interpreter','none')
		dim_other = find(dim-dim_avg_bool);
		dim_other_local = find(find(dim)==dim_other);
 		Xcut = squeeze(sum(X{dim_other_local},3-dim_avg_local)/n(dim_avg_local));
		for i=1:numvar
			subplot(ynum,xnum,plot_i(i));
			hold on 
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(Z,3-dim_avg_local)/n(dim_avg_local));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			xlabel(dim_label(dim_other_local));
			title(info.variable_names{i},'Interpreter','none');
			hold off
		end
	end

elseif (sum(dim)==3)
	%requires at least one averaging dimension
	if (~isfield(info,'dim_avg')), info.dim_avg = 2; end
  if (~isfield(info,'avg_not_cut')), info.avg_not_cut = true; end

	if (length(info.dim_avg)==1)
		%plotting 2D averaged data
		dim_avg = info.dim_avg(1);

% 		if (info.avg_not_cut)
% 			set(htextbox,'String',[...
% 				'plotting ' dim_label(dim_avg) '-averaged 3-dimensional data' ...
% 				info_text],'Interpreter','none')
% 		else
% 			set(htextbox,'String',[...
% 				'plotting center ' dim_label(dim_avg) '-cut of 3-dimensional data' ...
% 				info_text],'Interpreter','none')
% 		end

		dim_other = [1:1:dim_avg-1 dim_avg+1:1:3];
		Xcut{1} = squeeze(sum(X{dim_other(1)},4-dim_avg)/n(dim_avg));
		Xcut{2} = squeeze(sum(X{dim_other(2)},4-dim_avg)/n(dim_avg));

		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); 
			hold on
			Z(:) = data.variables(:,i);
			if (info.avg_not_cut)
				Zcut = squeeze(sum(Z,4-dim_avg)/n(dim_avg));
			else
				Zperm = permute(Z,4-[dim_other(end:-1:1) dim_avg]);
				Zcut = Zperm(:,:,round(end/2));
			end
			binplot(Xcut{1},Xcut{2},Zcut);
			colorbar;
			shading flat;
			axis equal;
			title(info.variable_names{i},'Interpreter','none');
			xlabel(dim_label(dim_other(1)));
			ylabel(dim_label(dim_other(2)));
			hold off
		end
	else
		%plotting 1D averaged data

% 		set(htextbox,'String',[...
% 		'plotting ' dim_label(info.dim_avg(1)) '-' dim_label(info.dim_avg(2)) '-averaged 3-dimensional data' ...
% 		info_text],'Interpreter','none');
% 		
		dim_other =  [1 2 3]~=info.dim_avg(1) & [1 2 3]~=info.dim_avg(2);
		Xcut =  squeeze(sum(sum(X{dim_other},4-info.dim_avg(2)),4-info.dim_avg(1)) ...
			/prod(n(info.dim_avg)));

		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); 
			hold on
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(sum(Z,4-info.dim_avg(2)),4-info.dim_avg(1)) ...
				/prod(n(info.dim_avg)));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			%axis tight
			xlabel(dim_label(dim_other));
			title(info.variable_names{i},'Interpreter','none');
			hold off
		end
		
	end
	
end

if (numvar==16)
  %scale
  for i=[1:16]
    subplot(4,4,i);
    axis tight
  end
  % PVF vertical range: [0:0.7]
  subplot(4,4,1);
  v = axis;
  axis([v(1:2) 0 .7]);
  % PVF vertical range: [0:0.7]
  subplot(4,4,4);
  v = axis;
  axis([v(1:3) .66]);
  % vy and vz with the same vertical axis range, e.g. [-0.02:0.02]
  same_scale([9 13],xnum,ynum);
  % fix vertical axis range for all diagonal stress: e.g. [-20:200]
  same_scale([6 11 16],xnum,ynum);
  % use the same range for Sxy and Syz: e.g. [-10:10]
  same_scale([7 12],xnum,ynum);
  % pls. zoom in [40:50] or [35:55] and plot a line at 45 degrees
  subplot(4,4,10);
  v = axis;
  hold on 
  plot(v(1:2),[45 45])
  %axis([v(1:2) 40 50]);
  subplot(4,4,14);
  v = axis;
  hold on 
  plot(v(1:2),[90 90])
  %axis([v(1:2) 80 91]);
  subplot(4,4,15);
  v = axis;
  hold on 
  plot(v(1:2),[45 45])
  %axis([v(1:2) 40 50]);
elseif (numvar==23 || numvar==24)
  same_scale([4:6],xnum,ynum);
  same_scale([7:12],xnum,ynum);
  same_scale([13:18],xnum,ynum);
  same_scale([19:24],xnum,ynum);
end

return


function binplot(X,Y,Z)
% plots 2D bin data better that surf or pcolor, plotting all data and
% plotting in the right xy position

dx = X(1,2)-X(1,1);
XX = [X([1,1:end],:)-dx/2, X([1,1:end],end)+dx/2];

dy = Y(2,1)-Y(1,1);
YY = [Y(:,[1,1:end])-dy/2; Y(end,[1,1:end])+dy/2];

ZZ = Z([1:end end],[1:end end]);

pcolor(XX,YY,ZZ);

return

function [data] = make_readable(data, info)
disp('writing variables into separate columns');

%not correct!
sorted = sortrows(data.coordinates);

% global_variables = textscan(info.text,'%s ')
% box(1) = str2double(global_variables{1}{6});
% box(2) = str2double(global_variables{1}{7});
% box(3) = str2double(global_variables{1}{8});
% boxvolume = prod(box)/length(data.coordinates);
 
%set nan values to zero
data.variables(isnan(data.variables))=0;

%find average of each variable
avg = sum(data.variables)/length(data.coordinates);

%store average
for i=1:size(data.variables,2);
  data = setfield(data, ['avg_' sscanf(info.variable_names{i},'%s ',16)], avg(i) );
end

return


function same_scale(plots,ynum,xnum)
%use same scale

subplot(ynum,xnum,plots(1));
axis tight
a = axis;
for i = plots(2:end)
  subplot(ynum,xnum,i);
  axis tight
  a_local = axis;
  if (a_local(3)<a(3)); a(3) = a_local(3); end
  if (a_local(4)>a(4)); a(4) = a_local(4); end
end
for i = plots
  subplot(ynum,xnum,i);
  axis(a);
end

return

