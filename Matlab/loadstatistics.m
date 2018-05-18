% loadstatistics (version 1.0) by Thomas Weinhart
%
% Use this function to load a MercuryDPM .stat file into a MATLAB struct.
% See http://mercurydpm.org/documentation/MercuryCG for more information.
%
% Usage: 
%   data=loadstatistics('chute.stat');
%       loads a single stat file 
%   data=loadstatistics({'chute.stat','hopper.stat'})
%       loads several stat files 
%   data=loadstatistics('*.stat');
%       loads all stat files in this directory
%   data=loadstatistics('chute.data');
%       loads a MDCLR data file (Stefan Luding's particle simulation code) 
%   data=loadstatistics(..,'basic');
%       Does not compute secondary cg fields. 

% Calls loadFile for each stat file specified. 
% If multiple stat files are specified, the individual outputs are grouped in a cell array. 
function data=loadstatistics(filenames,opt)
% creates an empty structure opt, if no optional arguments are given
if ~exist('opt','var'); opt=struct(); end

% Checks if no filename is given
if ~exist('filenames','var')
   error('Function requires at least stat file to be specified, e.g. loadstatistics(''chute.stat'')')
end

% Checks if multiple files are specified by a cell array of strings, 
% e.g. loadstatistics({'chute.stat','hopper.stat'})
if iscell(filenames) 
    % loadFile is used to load each individual file. 
    % The the multiple outouts are then grouped in a cell array.
	data=cell(size(filenames));
	for i=1:length(filenames)
		data{i} = loadFile(filenames{i},opt);
	end
else 
    % Checks if multiple files are specified by a basic regular expression, 
    % e.g. loadstatistics('*.stat') or loadstatistics('chute?.stat')
	if (~isempty(strfind(filenames,'?'))||~isempty(strfind(filenames,'*'))) 
        % If the filename contains * or ?, ls is used to extract a cell 
        % array of filenames.
        % Then loadstatistics is called with the cell array as argument 
		data = loadstatistics(strread(ls(filenames),'%s '),opt); 
    else
        % If argument is a single filename, loadFile is called to load this file. 
		data = loadFile(filenames,opt);
	end
end

return

% Loads a single .stat or .data file
function data = loadFile(filename,opt)
disp(['loading data from ' filename]);

% Loads the primary coarse-grained variables from the stat file.
% Distinguishes if the file to be loaded is a MercuryDPM stat file
% or a MDCLR data file
if (strcmp(filename(end-4:end),'.data'))
	data = loadDataFile(filename);
else
	data = loadStatFile(filename);
end

% makeReadable renames and reshapes coordinates and variables
% from simple 1D-arrays to meshgrid-like tensors (easier for plotting)
% getStandardVariables is called to add secondary coarse-grained
% fields; this can be disabled by the optional variable 'basic'
% If the output is a cell array (i.e. a sequence of time steps), these
% steps are done for each time step separately.
if (~iscell(data))
    data = makeReadable(data);
    if ~isfield(opt,'basic')
        data = getStandardVariables(data,opt);
    end
else
    for i=1:length(data)    
        data{i} = makeReadable(data{i});
        if ~isfield(opt,'basic')
            data{i} = getStandardVariables(data{i},opt);
        end
    end
end
return

% Loads the primary variables stored in a stat file
function data = loadStatFile(filename)

% Extracts the base string (the part before the .stat extension) from the 
% file name and stores it in the output structure.
data.name = filename(1:strfind(filename,'.stat')-1);

% Load the stat file using importdata, with the first three lines as header files
% (such that matlab recognizes the array size)
rawdata = importdata(filename,' ',3);

% \todo Not sure why this is here
if iscell(rawdata), return; end

% Checks if the stat file is from a very early version of MercuryDPM
if (size(rawdata.data,2)==26)
  % Loads old .stat files
  disp('WARNING:outdated stat file')
  rawdata.data = [ ...
    rawdata.data(:,1:14) ...
    zeros(size(rawdata.data,1),3) ...
    rawdata.data(:,15:end) ...
    zeros(size(rawdata.data,1),10)];
  data.variable_names = { ...
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
else
  % Loads a current MercuryDPM stat file
  % Reads variable names from line 1 and text output from line 2
  variable_names = textscan(rawdata.textdata{1},'%s ');
  data.variable_names = variable_names{1};
end

% Reads additional fields from the restart and ene files if they exist
% \todo make restart optional
data = readRestartFile(filename,data);
data = readEneFile(filename,data);

% Reads header of stat file and stores it in the output structure
text = textscan(rawdata.textdata{2},'%s ');
data.text = text{1};

% Checks how many data sets are in the stat file 
% (data sets begin with a time interval followed by coordinates and data values) 
% and stores an array of time intervals in the output structure 
% Writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,3))); size(rawdata.data,1)+1];
data.time = cell2mat(textscan(rawdata.textdata{3},'%f '));

% Reads the cg length scale from the header information
data.w = cell2mat(textscan(rawdata.textdata{2}(3:end),'%f '));

% Reads from the header information if gradient/variance need to be computed
% (for simplicity, it is assumed that you never compute both variance and gradient)
doGradient=any(cellfun(@(d)strcmp(d,'doGradient'),data.text));
doVariance=any(cellfun(@(d)strcmp(d,'doVariance'),data.text));

% Reads the coordinates from the data of the first time step 
% (assumes coordinates don't change with time)
data.coordinates = rawdata.data(1:index_time(1)-1,1:3);

% checks if there is only one data set 
if (length(index_time)==1)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
% checks if there is only one time step, plus variance data 
elseif (doVariance&&length(index_time)==2)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
    data.variance = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
% checks if there is only one time step, plus gradient data 
elseif (doGradient&&length(index_time)==4)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
    data.gradx = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
    data.grady = rawdata.data(index_time(2)+1:index_time(3)-1,4:end);
    data.gradz = rawdata.data(index_time(3)+1:index_time(4)-1,4:end);
% The following cases have all multiple time step, thus the data output is a cell array.
% First case: multiple time steps, no gradient data
elseif ~doGradient
    dataTemplate = data;
    data = cell(1,length(index_time)-1);
    data{1}=dataTemplate;
    data{1}.variables = rawdata.data(1:index_time(1)-1,4:end);
    for i=2:length(index_time)-1
        data{i}=dataTemplate;
        data{i}.time = rawdata.data(index_time(i-1),1:2)';
        data{i}.variables = rawdata.data(index_time(i-1)+1:index_time(i)-1,4:end);
    end
    disp(['multiple time steps (' num2str(length(index_time)-1) '); creating cell output'])
% Second case: multiple time steps with gradient data
else
    dataTemplate = data;
    data = cell(1,length(index_time)/4);
    data{1}=dataTemplate;
    data{1}.variables = rawdata.data(1:index_time(1)-1,4:end);
    data{1}.gradx = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
    data{1}.grady = rawdata.data(index_time(2)+1:index_time(3)-1,4:end);
    data{1}.gradz = rawdata.data(index_time(3)+1:index_time(4)-1,4:end);
    for i=8:4:length(index_time)
        data{i/4}=dataTemplate;
        data{i/4}.time = rawdata.data(index_time(i-4),1:2)';
        data{i/4}.variables = rawdata.data(index_time(i-4)+1:index_time(i-3)-1,4:end);
        data{i/4}.gradx = rawdata.data(index_time(i-3)+1:index_time(i-2)-1,4:end);
        data{i/4}.grady = rawdata.data(index_time(i-2)+1:index_time(i-1)-1,4:end);
        data{i/4}.gradz = rawdata.data(index_time(i-1)+1:index_time(i  )-1,4:end);
    end
    disp(['multiple time steps (' num2str(length(index_time)/4) '); creating cell output'])
% \todo the case multiple time steps with variance data is never considered
end
return

% Loads the primary variables stored in a data file (i.e. MDCLR output)
function data = loadDataFile(filename)
data.name = filename(1:end-5);

% Load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);

% reads variable names from line 1 and text output
% from line 2
data.variable_names = { ...
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
data.text = [text{1}{5} num2str(str2double(text{1}{9})*100)];
rho_text = textscan(rawdata.textdata{11},'%s ');
rho = str2double(rho_text{1}{end});

% writes the raw data into bits of data for each timestep
%index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data.time = 0;
data.coordinates = rawdata.data(:,4:6);
nu = max(rawdata.data(:,10),1e-200*ones(size(rawdata.data,1),1));
data.variables = [nu ...
	nu*rho ...
	rawdata.data(:,51:53).*nu(:,[1 1 1])*rho ...
	rawdata.data(:,[41:43 45 46 49]) ...
	nan(size(rawdata.data,1),3) ...
	rawdata.data(:,[21:23 25 26 29]) ...
	rawdata.data(:,[31:33 35 36 39]) ...
	rawdata.data(:,[71:73 75 76 79]) ...
	nan(size(rawdata.data,1),4) ...
	];

%data.description = ['Luding ' data.text];

return

% Renames coordinates to x,y,z and variables to Density, Momentum, ...
% Then reshapes coordinates and variables from simple 1D-arrays to 
% meshgrid-like tensors (easier for plotting)
function [data] = makeReadable(data)

% determines what shape to use
data.nx=length(data.coordinates(:,1))/sum(data.coordinates(:,1)==data.coordinates(1,1));
data.ny=length(data.coordinates(:,2))/sum(data.coordinates(:,2)==data.coordinates(1,2));
data.nz=length(data.coordinates(:,3))/sum(data.coordinates(:,3)==data.coordinates(1,3));
shape=size(squeeze(zeros(data.nz,data.ny,data.nx)));
if (prod(shape(:))~=size(data.coordinates,1))
    disp('Warning: cannot put xyz values on mesh')
    shape=[size(data.coordinates,1) 1];
end

% Renames and reshapes the coordinate values
data.x = reshape(data.coordinates(:,1),shape);
data.y = reshape(data.coordinates(:,2),shape);
data.z = reshape(data.coordinates(:,3),shape);
data = rmfield(data,{'coordinates'});

% Renames and reshapes the variables
% The array 'variable' is expanded into smaller arrays, each containing one
% variable only (f.e. Density), with the field names defined in
% variable_names. 
% If available, variance, gradx, grady, gradz are also expanded the same way
if isfield(data,'variance')
	for i=1:length(data.variable_names)
  	data.([sscanf(data.variable_names{i},'%s ',16) '_var']) ...
      =reshape(data.variance(:,i),shape);
	end
	data = rmfield(data,{'variance'});
end

if isfield(data,'gradx')
	for i=1:length(data.variable_names)
  	data.([sscanf(data.variable_names{i},'%s ',16) '_dx']) ...
      =reshape(data.gradx(:,i),shape);
  	data.([sscanf(data.variable_names{i},'%s ',16) '_dy']) ...
      =reshape(data.grady(:,i),shape);
  	data.([sscanf(data.variable_names{i},'%s ',16) '_dz']) ...
      =reshape(data.gradz(:,i),shape);
	end
	data = rmfield(data,{'gradx','grady','gradz'});
end

for i=1:length(data.variable_names)
	data.(sscanf(data.variable_names{i},'%s ',16)) ...
    = reshape(data.variables(:,i),shape);
end
data = rmfield(data,{'variables','variable_names'});
return

% Extracts secondary cg fields (fields computed from the primary fields 
% in the stat file), such as kinetic stress, velocity, granular temperature, 
% friction.
% If a division (e.g. for  Velocity=Momentum/Density) causes nan values, 
% these values are set to zero
function data = getStandardVariables(data,opt)

data.MomentumX(isnan(data.MomentumX)) = 0;
data.MomentumY(isnan(data.MomentumY)) = 0;
data.MomentumZ(isnan(data.MomentumZ)) = 0;

data.VelocityX = data.MomentumX./data.Density;
data.VelocityY = data.MomentumY./data.Density;
data.VelocityZ = data.MomentumZ./data.Density;

data.VelocityX(isnan(data.VelocityX)) = 0;
data.VelocityY(isnan(data.VelocityY)) = 0;
data.VelocityZ(isnan(data.VelocityZ)) = 0;

data.Temperature = (...
	(data.MomentumFluxXX + data.MomentumFluxYY + data.MomentumFluxZZ)./data.Density ...
	- (data.VelocityX.^2 + data.VelocityY.^2 + data.VelocityZ.^2) )/3;
data.Temperature(isnan(data.Temperature)) = 0;

data.CoordinationNumber=(data.FabricXX+data.FabricYY+data.FabricZZ)./data.Nu;
data.CoordinationNumber(isnan(data.CoordinationNumber)) = 0;

if isfield(opt,'basic'); return; end

% \todo rename traction to interaction force density or drag
data.TractionX = data.NormalTractionX + data.TangentialTractionX;
data.TractionY = data.NormalTractionY + data.TangentialTractionY;
data.TractionZ = data.NormalTractionZ + data.TangentialTractionZ;
data.KineticStressXX = data.MomentumFluxXX - data.Density.*data.VelocityX.*data.VelocityX;
data.KineticStressXY = data.MomentumFluxXY - data.Density.*data.VelocityX.*data.VelocityY;
data.KineticStressXZ = data.MomentumFluxXZ - data.Density.*data.VelocityX.*data.VelocityZ;
data.KineticStressYY = data.MomentumFluxYY - data.Density.*data.VelocityY.*data.VelocityY;
data.KineticStressYZ = data.MomentumFluxYZ - data.Density.*data.VelocityY.*data.VelocityZ;
data.KineticStressZZ = data.MomentumFluxZZ - data.Density.*data.VelocityZ.*data.VelocityZ;
data.StressXX = data.NormalStressXX + data.TangentialStressXX + data.KineticStressXX;
data.StressXY = data.NormalStressXY + data.TangentialStressXY + data.KineticStressXY; 
data.StressXZ = data.NormalStressXZ + data.TangentialStressXZ + data.KineticStressXZ;
data.StressYX = data.NormalStressYX + data.TangentialStressYX + data.KineticStressXY;
data.StressYY = data.NormalStressYY + data.TangentialStressYY + data.KineticStressYY;
data.StressYZ = data.NormalStressYZ + data.TangentialStressYZ + data.KineticStressYZ;
data.StressZX = data.NormalStressZX + data.TangentialStressZX + data.KineticStressXZ;
data.StressZY = data.NormalStressZY + data.TangentialStressZY + data.KineticStressYZ;
data.StressZZ = data.NormalStressZZ + data.TangentialStressZZ + data.KineticStressZZ;
data.Pressure =  (data.StressXX + data.StressYY + data.StressZZ)/3;

%macroscopic friction coefficient, or deviator-stress ratio sD1:=(S1-S3)/(2*p)  or sD2 (to be discussed)
data.MacroFrictionCoefficient = zeros(size(data.x));
data.Sd = zeros(size(data.x));
data.SdStar = zeros(size(data.x));
for j=1:length(data.x(:))
	Stress = [ ...
		data.StressXX(j) data.StressXY(j) data.StressXZ(j); ...
		data.StressXY(j) data.StressYY(j) data.StressYZ(j); ...
		data.StressXZ(j) data.StressYZ(j) data.StressZZ(j) ];
	if (isequal(Stress,zeros(3))|| sum(sum(isnan(Stress))) )
		data.MacroFrictionCoefficient(j) = NaN;
		data.Sd(j) = nan;
		data.SdStar(j) = nan;
	else
		[v1,d1]=eig(Stress);
		[v,d]=dsort_ev(v1,d1);
        d = diag(d);
        p=sum(d)/3;
		data.MacroFrictionCoefficient(j) = -data.StressXZ(j)/data.StressZZ(j);
		data.Sd(j) = (d(1)-d(3))/(2*p);
		data.SdStar(j) = sqrt(((d(1)-d(3))^2+(d(2)-d(3))^2+(d(1)-d(2))^2) / (6*p^2));
	end
end

% computes gradients of the secondary fields if doGradient is true
if (isfield(data,'Nu_dz'))
  for i='XYZ'
  for k='xyz'
  vector=[i '_d' k];
  data.(['Velocity' vector]) = ...
      (data.(['Momentum' vector]).*data.Density-data.(['Momentum' i]).*data.(['Density_d' k]))./data.Density.^2;
  data.(['Velocity' vector])(isnan(data.(['Velocity' vector]))) = 0;
  end
  end
  
  for i='XYZ'
  for j='XYZ'
  for k='xyz'
  tensor=[i j '_d' k];
  if i<j, symtensor=tensor; else symtensor=tensor([2 1 3:end]); end
  data.(['KineticStress' tensor]) = data.(['MomentumFlux' symtensor]) ...
      - data.(['Momentum' i '_d' k]).*data.(['Velocity' j]) ...
      - data.(['Momentum' j '_d' k]).*data.(['Velocity' i]);
  %data.KineticStressXZ_dz = data.MomentumFluxXZ_dz - data.Density_dz.*data.VelocityX.*data.VelocityZ - data.Density.*data.VelocityX_dz.*data.VelocityZ - data.Density.*data.VelocityX.*data.VelocityZ_dz;
  data.(['Stress' tensor]) = data.(['NormalStress' tensor]) + data.(['TangentialStress' tensor]) + data.(['KineticStress' tensor]);
  end
  end
  end
end

% Uncomment these lines to compute 
% - extra variables that only exist for stattype Z (this was implemented for the chute flow)
% - the remainder of the steady-state momentum equation
% if (data.nz>1&&data.nx==1&&data.ny==1), 
%     data = getDepthVariables(data); 
%     [data.MomentumEquationsRemainder, data.MomentumEquationsMaximum]=getMomentumEquation(data);
% end
return

% extracts extra variables that only exist for stattype Z
function data = getDepthVariables(data)

dz=diff(data.z(1:2,1));
data.ShearRate = diff(data.VelocityX)./dz;
data.ShearRate(end+1) = data.ShearRate(end);
data.InertialParameter=data.ShearRate*data.d./sqrt(data.StressZZ/data.ParticleDensity);

data.FlowNu = sum(data.Nu,1)*dz/diff(data.Domain([5 6]));
data.FlowDensity = sum(data.Density,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumX = sum(data.MomentumX,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumY = sum(data.MomentumY,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumZ = sum(data.MomentumZ,1)*dz/diff(data.Domain([5 6]));
data.FlowVelocityX = data.FlowMomentumX./data.FlowDensity;
data.FlowVelocityY = data.FlowMomentumY./data.FlowDensity;
data.FlowVelocityZ = data.FlowMomentumZ./data.FlowDensity;

data.FlowMomentumX = sum(data.MomentumX,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumY = sum(data.MomentumY,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumZ = sum(data.MomentumZ,1)*dz/diff(data.Domain([5 6]));
data.FlowVelocityX = data.FlowMomentumX./data.FlowDensity;
data.FlowVelocityY = data.FlowMomentumY./data.FlowDensity;
data.FlowVelocityZ = data.FlowMomentumZ./data.FlowDensity;
FlowMomentumFluxXX = sum(data.MomentumFluxXX,1)*dz/diff(data.Domain([5 6]));
FlowMomentumFluxYY = sum(data.MomentumFluxYY,1)*dz/diff(data.Domain([5 6]));
FlowMomentumFluxZZ = sum(data.MomentumFluxZZ,1)*dz/diff(data.Domain([5 6]));

data.BottomFrictionCoefficient = -data.StressXZ(1,:)/data.StressZZ(1,:);

% height and average density of the flow
if ~(data.nx==1&&data.ny==1), 
  ContactStressZZ=data.NormalStressZZ+data.TangentialStressZZ;
  data.FlowHeight=max(data.z.*(ContactStressZZ~=0),[],1);
else
  kappa=0.02;
  IntDensity = cumsum(data.Density);
  Base = min(data.z(IntDensity>=kappa*max(IntDensity)));
  Surface = max(data.z(IntDensity<=(1-kappa)*max(IntDensity)));
  data.FlowHeight=(Surface-Base)/(1-2*kappa);
  data.Base=Base-kappa*data.FlowHeight;
  data.Surface=Surface+kappa*data.FlowHeight;

  FlowHeightIndex = data.z>data.Base & data.z<data.Base+data.FlowHeight;
  if ~isempty(FlowHeightIndex)
    data.Froude=norm([data.FlowVelocityX data.FlowVelocityY data.FlowVelocityZ])/sqrt(data.FlowHeight*(-data.Gravity(3)));
  end
end

% % height and average density of the flow
% FlowHeightIndex = min([length(data.Nu),round(length(data.Nu)/2)-1+find(data.Nu(round(end/2):end)<0.3,1)]);
% if ~isempty(FlowHeightIndex)
% 	AveragingInterval = ceil(FlowHeightIndex/3):ceil(FlowHeightIndex*2/3);
% 	data.FlowDensity = sum(data.Density(AveragingInterval))/length(AveragingInterval);
% 	data.FlowNu = sum(data.Nu(AveragingInterval))/length(AveragingInterval);
% 	data.FlowHeight = data.z(1) + sum(data.Density)/length(data.z)/data.FlowDensity * (data.z(end)-data.z(1));
% 
%   index = data.z>0&data.z<data.FlowHeight;
% 	data.FlowVelocityX = average(data.VelocityX(index));
% 	data.FlowVelocityY = average(data.VelocityY(index));
% 	data.FlowVelocityZ = average(data.VelocityZ(index));
%   data.Froude=norm([data.FlowVelocityX data.FlowVelocityY data.FlowVelocityZ])/sqrt(data.FlowHeight*(-data.Gravity(3)));
% end
return

function [Vs, Ds] = dsort_ev(V, D)
    Ds=D;
    Vs=V;
%     tmpvec=[0;0;0];
    if((Ds(1,1))<(Ds(2,2)))
        tmp=Ds(1,1);
        Ds(1,1)=Ds(2,2);
        Ds(2,2)=tmp;
        tmpvec=Vs(:,1);
        Vs(:,1)=Vs(:,2);
        Vs(:,2)=tmpvec;
    end
    if((Ds(2,2))<(Ds(3,3)))
        tmp=Ds(2,2);
        Ds(2,2)=Ds(3,3);
        Ds(3,3)=tmp;
        tmpvec=Vs(:,2);
        Vs(:,2)=Vs(:,3);
        Vs(:,3)=tmpvec;
    end
    if((Ds(1,1))<(Ds(2,2)))
        tmp=Ds(1,1);
        Ds(1,1)=Ds(2,2);
        Ds(2,2)=tmp;
        tmpvec=Vs(:,1);
        Vs(:,1)=Vs(:,2);
        Vs(:,2)=tmpvec;
    end
    if((Ds(2,2))<(Ds(3,3)))
        tmp=Ds(2,2);
        Ds(2,2)=Ds(3,3);
        Ds(3,3)=tmp;
        tmpvec=Vs(:,2);
        Vs(:,2)=Vs(:,3);
        Vs(:,3)=tmpvec;
    end

return

function [remainder,maximum] = getMomentumEquation(data)
if isfield(data,'Nu_dx'),
  NablaMomentumX = data.ParticleDensity*[data.MomentumX_dx data.MomentumX_dy data.MomentumX_dz];
  NablaMomentumY = data.ParticleDensity*[data.MomentumY_dx data.MomentumY_dy data.MomentumY_dz];
  NablaMomentumZ = data.ParticleDensity*[data.MomentumZ_dx data.MomentumZ_dy data.MomentumZ_dz];
  NablaDotStress = [...
      data.StressXX_dx+data.StressXY_dy+data.StressXZ_dz ...
      data.StressYX_dx+data.StressYY_dy+data.StressYZ_dz ...
      data.StressZX_dx+data.StressZY_dy+data.StressZZ_dz ];
else
  disp('estimating gradient!')
  NablaMomentumX = nabla(data.ParticleDensity*data.Nu.*data.VelocityX,data.x,data.y,data.z);
  NablaMomentumY = nabla(data.ParticleDensity*data.Nu.*data.VelocityY,data.x,data.y,data.z);
  NablaMomentumZ = nabla(data.ParticleDensity*data.Nu.*data.VelocityZ,data.x,data.y,data.z);

  NablaDotStress = nabla([data.StressXX data.StressXY data.StressXZ data.StressYY data.StressYZ data.StressZZ],data.x,data.y,data.z);
end
    
VelocityDotNablaMomentum = [...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumX,2) ...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumY,2) ...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumZ,2) ];

DensityGravity = data.ParticleDensity * data.Nu * data.Gravity;

Traction = [data.TractionX data.TractionY data.TractionZ ];

remainder = VelocityDotNablaMomentum -DensityGravity +NablaDotStress +Traction;

maximum = max([max(VelocityDotNablaMomentum) max(DensityGravity) max(NablaDotStress)]);
return

% Adds the ene file information to the Matlab output 
% (only called if ene file is found)
function data = readEneFile(statname,data)

%load ene data
filename = [statname(1:end-5) '.ene'];
if ~exist(filename,'file'), 
  %if unsuccessful, load ene data with last part of filename removed
  dots=strfind(statname,'.');
  if (length(dots)>1), filename = [statname(1:min(dots(end-1:end))) 'ene']; end
  if ~exist(filename,'file')
    warning([filename ' not found; the output will not contain the variables stored in the ene file']);
  else
    %if unsuccessful, give up
    disp([statname(1:end-5) '.ene not found, using ' filename ' instead']); 
  end
end

if exist(filename,'file'),
  % load raw data from file, with the first three lines as header files
  % (such that matlab recognizes the array size)
  rawdata = importdata(filename,' ',1);
  % if tabs are used as delimiters
  if (~isfield(rawdata,'data')||size(rawdata.data,2)~=8)
    rawdata = importdata(filename,'\t',1);
  end
  if (~isfield(rawdata,'data')||size(rawdata.data,2)~=8)
      disp([filename ' not readable']);
      return;
  end

  data.Ene.Time = rawdata.data(:,1);
  data.Ene.Gra = rawdata.data(:,2);
  data.Ene.Kin = rawdata.data(:,3);
  data.Ene.Rot = rawdata.data(:,4);
  data.Ene.Ela = rawdata.data(:,5);
  data.Ene.ComX = rawdata.data(:,6);
  data.Ene.ComY = rawdata.data(:,7);
  data.Ene.ComZ = rawdata.data(:,8);
else
  disp([filename ' not found'])
end

return

% Adds the restart file information to the Matlab output
function data = readRestartFile(statname,data)

%load restart data
filename = [statname(1:end-4) 'restart'];
if ~exist(filename,'file'), 
  %if unsuccessful, load ene data with last part of filename removed
  dots=strfind(statname,'.');
  if (length(dots)>1), filename = [statname(1:min(dots(end-1:end))) 'restart']; end
  if ~exist(filename,'file')
    warning([filename ' not found; the output will not contain variables such as particle size and gravity']);
    return
  else
    %if unsuccessful, give up
    disp([statname(1:end-4) 'restart not found, using ' filename ' instead']); 
  end
end

%load restart data
fid=fopen(filename);
tdata = textscan(fid,'%s');
tdata = tdata{1};
fclose(fid);

%this assumes monodispersed particles
if (length(tdata{1})==1) 
	%old restart files
	data.ParticleDensity = str2double(tdata(29));
	data.d = 2*str2double(tdata(28));
	data.Gravity = str2double(tdata(2:4))';
	data.N = str2double(tdata(end-3));
	data.Domain=str2double(tdata(5:10))';
else
	%new restart version
	i = find(strcmp(tdata,'density'));
    if isempty(i), i = find(strcmp(tdata,'rho')); end
  %if length(i)>1, disp('warning: two species'); end
	data.ParticleDensity = str2double(tdata(i+1));
	i = find(strcmp(tdata,'slidingFrictionCoefficient'));
    if isempty(i), i = find(strcmp(tdata,'mu')); end
	data.Mu = str2double(tdata(i+1));
	i = find(strcmp(tdata,'gravity'));
	data.Gravity = str2double(tdata(i+(1:3)))';
	i = find(strcmp(tdata,'Particles'));
	data.N = str2double(tdata{i+1});
	i = find(strcmp(tdata,'xMin'));
	if isempty(i), i = find(strcmp(tdata,'xmin')); end
	data.Domain = str2double(tdata(i+(1:2:11)))';
	i = find(strcmp(tdata,'FixedParticleRadius'));
	if ~isempty(i), data.FixedParticleRadius = str2double(tdata{i+1}); end
	i = find(strcmp(tdata,'MinInflowParticleRadius'));
	if ~isempty(i), 
		data.InflowParticleRadius = str2double(tdata(i+[1 3]))'; 
		data.d = sum(data.InflowParticleRadius);
	else
		i = find(strcmp(tdata,'Particles'));
    if (str2double(tdata{i+1})==0)
      data.d=nan;
    else
		j = min(find(strcmp(tdata(i+2:end),'radius')));
        data.d = 2*str2double(tdata{i+j+2});
    end
	end
end
data.ParticleVolume = pi/6*data.d^3;
data.DomainVolume = prod(data.Domain([2 4 6])-data.Domain([1 3 5]));
data.ChuteAngle = round(atand(-data.Gravity(1)/data.Gravity(3))*400)/400;
return

function NablaVariable = nabla(variable,x,y,z)

% stores the grid dimensions
n(1) = length(x)/sum(x(1)==x);
n(2) = length(y)/sum(y(1)==y);
n(3) = length(z)/sum(z(1)==z);

% stores coordinates in grid
X = zeros(n(end:-1:1));
X(:) = x;
Y = zeros(n(end:-1:1));
Y(:) = y;
Z = zeros(n(end:-1:1));
Z(:) = z;
V = zeros(n(end:-1:1));

% stores differentials (small number if X/Y/Z is constant)
dX = max(1e-60,(X(:,:,[2:end end])-X(:,:,:)));
dY = max(1e-60,(Y(:,[2:end end],:)-Y(:,:,:)));
dZ = max(1e-60,(Z([2:end end],:,:)-Z(:,:,:)));

if size(variable,2)==1 %scalar variable
  V(:) = variable;
  dXV = (V(:,:,[2:end end])-V(:,:,:)) ./dX;
  dYV = (V(:,[2:end end],:)-V(:,:,:)) ./dY;
  dZV = (V([2:end end],:,:)-V(:,:,:)) ./dZ;
  
  NablaVariable = [dXV(:) dYV(:) dZV(:)];
elseif size(variable,2)==3 %vector
  V(:) = variable(1,:);
  dXV = (V(:,:,[2:end end])-V(:,:,:)) ./dX;
  V(:) = variable(2,:);
  dYV = (V(:,[2:end end],:)-V(:,:,:)) ./dY;
  V(:) = variable(3,:);
  dZV = (V([2:end end],:,:)-V(:,:,:)) ./dZ;

  NablaVariable = dXV(:) + dYV(:) + dZV(:);
elseif size(variable,2)==6 %symmetric matrix
  V(:) = variable(:,1);
  dXVXX = (V(:,:,[2:end end])-V(:,:,:)) ./dX;
  V(:) = variable(:,2);
  dXVYX = (V(:,:,[2:end end])-V(:,:,:)) ./dX;
  V(:) = variable(:,3);
  dXVZX = (V(:,:,[2:end end])-V(:,:,:)) ./dX;
  
  V(:) = variable(:,2);
  dYVXY = (V(:,[2:end end],:)-V(:,:,:)) ./dY;
  V(:) = variable(:,4);
  dYVYY = (V(:,[2:end end],:)-V(:,:,:)) ./dY;
  V(:) = variable(:,5);
  dYVZY = (V(:,[2:end end],:)-V(:,:,:)) ./dY;
  
  V(:) = variable(:,3);
  dZVXZ = (V([2:end end],:,:)-V(:,:,:)) ./dZ;
  V(:) = variable(:,5);
  dZVYZ = (V([2:end end],:,:)-V(:,:,:)) ./dZ;
  V(:) = variable(:,6);
  dZVZZ = (V([2:end end],:,:)-V(:,:,:)) ./dZ;

  NablaVariable = [dXVXX(:)+dYVXY(:)+dZVXZ(:) dXVYX(:)+dYVYY(:)+dZVYZ(:) dXVZX(:)+dYVZY(:)+dZVZZ(:)];
else 
  disp('Error');
end
return

function NablaVariable = nablaCenter(variable,x,y,z)

% stores the grid dimensions
n(1) = length(x)/sum(x(1)==x);
n(2) = length(y)/sum(y(1)==y);
n(3) = length(z)/sum(z(1)==z);

% stores coordinates in grid
X = zeros(n(end:-1:1));
X(:) = x;
Y = zeros(n(end:-1:1));
Y(:) = y;
Z = zeros(n(end:-1:1));
Z(:) = z;
V = zeros(n(end:-1:1));

% stores differentials (small number if X/Y/Z is constant)
dX = max(1e-60,(X(:,:,[2:end end])-X(:,:,[1 1:end-1])));
dY = max(1e-60,(Y(:,[2:end end],:)-Y(:,[1 1:end-1],:)));
dZ = max(1e-60,(Z([2:end end],:,:)-Z([1 1:end-1],:,:)));

if size(variable,2)==1 %scalar variable
  V(:) = variable;
  dXV = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  dYV = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  dZV = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  
  NablaVariable = [dXV(:) dYV(:) dZV(:)];
elseif size(variable,2)==3 %vector
  V(:) = variable(1,:);
  dXV = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(2,:);
  dYV = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(3,:);
  dZV = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;

  NablaVariable = dXV(:) + dYV(:) + dZV(:);
elseif size(variable,2)==6 %symmetric matrix
  V(:) = variable(:,1);
  dXVXX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(:,2);
  dXVYX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(:,3);
  dXVZX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  
  V(:) = variable(:,2);
  dYVXY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(:,4);
  dYVYY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(:,5);
  dYVZY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  
  V(:) = variable(:,3);
  dZVXZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  V(:) = variable(:,5);
  dZVYZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  V(:) = variable(:,6);
  dZVZZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;

  NablaVariable = [dXVXX(:)+dYVXY(:)+dZVXZ(:) dXVYX(:)+dYVYY(:)+dZVYZ(:) dXVZX(:)+dYVZY(:)+dZVZZ(:)];
else 
  disp('Error');
end
return
