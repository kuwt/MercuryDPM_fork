function [data, info] = load_stat_file(filename,num)
%load stat file data




disp(['loading data from ' filename]);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',3);

  
% reads variable names from line 1 and text output
% from line 2

variable_names = textscan(rawdata.textdata{1},'%s ');
info.variable_names = variable_names{1};



if (size(rawdata.data,2)==26)
  disp('WARNING:outdated stat file')
  rawdata.data = [rawdata.data(:,1:14) zeros(size(rawdata.data,1),3) rawdata.data(:,15:end) zeros(size(rawdata.data,1),10)];
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
end

text = textscan(rawdata.textdata{2},'%s ');
info.text = text{1};

% writes the raw data into bits of data for each timestep
index_time = [0; find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];


data.time = cell2mat(textscan(rawdata.textdata{3},'%d '));

data.coordinates = rawdata.data(index_time(num)+1:index_time(num+1)-1,1:3);
data.variables = rawdata.data(index_time(num)+1:index_time(num+1)-1,[4:end]);
if (length(index_time)>1)
  data.variance = rawdata.data(index_time(num)+1:index_time(num+1)-1,[4:end]);
end

info.description = ['Goldhirsch w/d=' num2str(str2double(info.text{2})/1e-3)];

return