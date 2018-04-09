function [data, info] = load_file(filename)
%load stat file data

if (strcmp(filename(end-4:end),'.data'))
	[data, info] = load_data_file(filename);
else
	[data, info] = load_stat_file(filename,1);
end

%add options
info.dim_avg = [1 2];
info.avg_not_cut = true;
info.hold = true;

%data.variables = data.variables(data.coordinates(:,3)<0.0188,:);
%data.coordinates = data.coordinates(data.coordinates(:,3)<0.0188,:);

return