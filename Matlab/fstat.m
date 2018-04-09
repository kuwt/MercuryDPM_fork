function [data,info]=fstat(filename)
%plots data from stat files
disp(' ');
close all;

% find latest stat file
if (~exist('filename'))
  filenames = ls('*.stat -t');
  filename = strtok(filenames);
end

%load file
[datalist, info] = load_file(filename);

%average data
data = get_average_data(datalist,60);

%add options
info.dim_avg = [1 2];
info.avg_not_cut = true;
%info.h = 2;

%plot data
plot_stat_file(data,info);

[data] = make_readable(data, info);

%load file
[datalist, info] = load_file('tav8.data');

%add options
info.dim_avg = [1 2];
info.hold = true;
info.plottype = 'r';

%plot data
plot_stat_file(datalist{1},info);
return


function [data, info] = load_file(filename)
%load stat file data
if (strcmp(filename(end-4:end),'.data'))
	[data, info] = load_data_file(filename);
else
	[data, info] = load_stat_file(filename);
end

return


function [data, info] = load_stat_file(filename)
%load stat file data
disp(['loading data from ' filename]);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',3);
  
% reads variable names from line 1 and text output
% from line 2
variable_names = textscan(rawdata.textdata{1},'%s ');
info.variable_names = variable_names{1};
info.text = rawdata.textdata{2};

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data{1}.time = cell2mat(textscan(rawdata.textdata{3},'%d '));
data{1}.coordinates = rawdata.data(1:index_time(1)-1,1:3);
data{1}.variables = rawdata.data(1:index_time(1)-1,4:end);
for i=2:length(index_time)
  data{i}.time = rawdata.data(index_time(i-1),1);
  data{i}.coordinates = rawdata.data(index_time(i-1)+1:index_time(i)-1,1:3);
  data{i}.variables = rawdata.data(index_time(i-1)+1:index_time(i)-1,4:end);
end

return

function [data, info] = load_data_file(filename)
%load stat file data
disp(['loading data from ' filename]);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);
  
% reads variable names from line 1 and text output
% from line 2
info.variable_names = {'PVF'; 'Velocity_x'; 'Velocity_y'; 'Velocity_z'; 'StressXX'; 'StressXY'; 'StressXZ'; 'StressYY'; 'StressYZ'; 'StressZZ'; 'Temperature'; 'Pressure'};
info.text = ' ';

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data{1}.time = 0;
data{1}.coordinates = rawdata.data(:,4:6);
data{1}.variables = [ ...
  rawdata.data(:,[12 51 52 53]), ...
  (rawdata.data(:,[21 22 23 25 26 29]) + rawdata.data(:,[31 32 33 35 36 39])), ...
  1/3*(sum(rawdata.data(:,[54 55 56]),2)-sum(rawdata.data(:,[51 52 53]).^2,2)) ...
  rawdata.data(:,[19]), ...
  ];

return

function plot_stat_file(data,info)
% plots data of one timestep

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

%create figure
figure(h);
%clear figure if info.hold does not exist
if (~isfield(info,'hold')), clf; end
%resize figure     
set(gcf,'Position',[300 1 1379 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
%plot info text
info_text = [info.text ' t=' mat2str(data.time)];
htextbox = annotation('textbox',[0 0 1 1]);
set(htextbox,'String',info_text ,'Interpreter','none')

% defines number of variables and plot order
numvar = size(data.variables,2);
subplotlength = ceil(sqrt(numvar));
plot_i = 1:numvar;
if (numvar<=4), 
elseif (numvar<=9), 
  plot_i = [1 4 7 5 6 9 2 3 8];
elseif (numvar>=11), 
  plot_i = [1 5 9 13 6 7 8 11 12 16 2 3 4];
end

if (sum(dim)==1)
  % for 1D output
  set(htextbox,'String',['plotting 1-dimensional data' info_text],'Interpreter','none')
  
	for i=1:numvar
    subplot(subplotlength,subplotlength,plot_i(i)); 
		hold on
    Z(:) = data.variables(:,i);
    plot(X{1},Z,info.plottype);
    plot(sum(X{1})/length(X{1}),sum(Z)/length(Z),'rx');
    title(info.variable_names{i});
    xlabel(dim_label);
		hold off
	end

elseif (sum(dim)==2)
  % for 2D output
  
	if (~(isfield(info,'dim_avg')))% && dim(info.dim_avg)) )
		%plotting without averaging
		set(htextbox,'String',['plotting 2-dimensional data' info_text],'Interpreter','none')
		for i=1:numvar
			subplot(subplotlength,subplotlength,plot_i(i)); hold on
			Z(:) = data.variables(:,i);
			binplot(X{1},X{2},Z);
			colorbar;
			shading flat;
			axis equal;
			title(info.variable_names{i});
			xlabel(dim_label(1));
			ylabel(dim_label(2));
		end
	else
		%plotting averaged data
 
		dim_avg = info.dim_avg(1);
    dim_avg_bool = zeros(1,3);
    dim_avg_bool(dim_avg) = true;
		dim_avg_local = find(find(dim)==dim_avg);
 		set(htextbox,'String',['plotting ' dim_label_full(dim_avg) '-averaged 2-dimensional data' info_text],'Interpreter','none')
		dim_other = find(dim-dim_avg_bool);
		dim_other_local = find(find(dim)==dim_other);
 		Xcut = squeeze(sum(X{dim_other_local},3-dim_avg_local)/n(dim_avg_local));
		for i=1:numvar
			subplot(subplotlength,subplotlength,plot_i(i));
			hold on 
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(Z,3-dim_avg_local)/n(dim_avg_local));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			xlabel(dim_label(dim_other_local));
			title(info.variable_names{i});
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

		if (info.avg_not_cut)
			set(htextbox,'String',[...
				'plotting ' dim_label(dim_avg) '-averaged 3-dimensional data' ...
				info_text],'Interpreter','none')
		else
			set(htextbox,'String',[...
				'plotting center ' dim_label(dim_avg) '-cut of 3-dimensional data' ...
				info_text],'Interpreter','none')
		end

		dim_other = [1:1:dim_avg-1 dim_avg+1:1:3];
		Xcut{1} = squeeze(sum(X{dim_other(1)},4-dim_avg)/n(dim_avg));
		Xcut{2} = squeeze(sum(X{dim_other(2)},4-dim_avg)/n(dim_avg));

		for i=1:numvar
			subplot(subplotlength,subplotlength,plot_i(i)); 
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
			title(info.variable_names{i});
			xlabel(dim_label(dim_other(1)));
			ylabel(dim_label(dim_other(2)));
			hold off
		end
	else
		%plotting 1D averaged data

		set(htextbox,'String',[...
		'plotting ' dim_label(info.dim_avg(1)) '-' dim_label(info.dim_avg(2)) '-averaged 3-dimensional data' ...
		info_text],'Interpreter','none');
		
		dim_other =  [1 2 3]~=info.dim_avg(1) & [1 2 3]~=info.dim_avg(2);
		Xcut =  squeeze(sum(sum(X{dim_other},4-info.dim_avg(2)),4-info.dim_avg(1)) ...
			/prod(n(info.dim_avg)));

		for i=1:numvar
			subplot(subplotlength,subplotlength,plot_i(i)); 
			hold on
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(sum(Z,4-info.dim_avg(2)),4-info.dim_avg(1)) ...
				/prod(n(info.dim_avg)));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			%axis tight
			xlabel(dim_label(dim_other));
			title(info.variable_names{i});
			hold off
		end
		
	end
	
end

saveas(h ,['figure' num2str(h ) '.pdf']);
return

function avg = get_average_data(data,minimum)

minimum = min(minimum,length(data));
avg.time = [data{minimum}.time data{end}.time];
avg.coordinates = data{minimum}.coordinates;
avg.variables = data{minimum}.variables;
for i=minimum+1:length(data), 
  avg.variables = avg.variables + data{i}.variables; 
end
avg.variables = avg.variables/(length(data)-minimum+1);

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
 
numvar = size(data.variables,2);
avg = sum(data.variables)/length(data.coordinates);
for i=1:numvar
  %data = setfield(data, info.variable_names{i}, data.variables(:,i) );
  data = setfield(data, ['avg_' info.variable_names{i}], avg(i) );
end
return