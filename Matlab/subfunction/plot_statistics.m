function [data,info]=plot_statistics(filenames)
%plots data from stat files

%clean up
disp(' '); close all; 
description=[];
plottype = {'y:','k','r','g','b','c','m','y','k:','r:','g:','b:','c:','m:'};
label = {'yellow...','black','red','green','blue','cyan','magenta', ...
  'yellow','black...','red...','green...','blue...','cyan...','magenta...'};

for i = 1:length(filenames)
  %load file
  [data, info] = load_file(filenames{i});
  [data, info] = get_standard_variables(data,info);
  %[data, info] = get_I(data,info);
  %[data, info] = extract_variables(data,info,[27:32]);
  %[data, info] = get_momentum_equation(data,info);
  
  description{end+1} = [ label{mod(i,14)+1} ': ' filenames{i} ', ' info.description ', '];

  %plot data
  info.plottype = plottype{mod(i,14)+1};
  plotdata(data,info);
end