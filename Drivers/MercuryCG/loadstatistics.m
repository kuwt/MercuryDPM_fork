% loadstatistics (version 1.0) by Thomas Weinhart
%
% loads a MercuryDPM stat (or MDCLR data) file into a matlab struct
%
% it accepts a single filename, a cell of filenames, or a string that will
% be inserted into ls
%
% Usages: 
%   data=loadstatistics('chuteDemo.stat');
%   data=loadstatistics('*.stat');
%   data=loadstatistics({'chuteDemo.stat','hopperDemo.stat'})
% 
function data=loadstatistics(filenames,opt)
if ~exist('opt','var'); opt=struct(); end

if iscell(filenames) 
    % if argument is a cell of filenames, load_file is used to load each of
    % the files
	data=cell(size(filenames));
	for i=1:length(filenames)
		data{i} = load_file(filenames{i},opt);
	end
else 
	if (~contains(filenames,'?')&&~contains(filenames,'*')) 
        % if argument is a single filename, load_file is used for loading
		data = load_file(filenames,opt);
    else
        % if argument contains * or ?, run it through ls to procuce a cell
        % of filenames 
		data = loadstatistics(strread(ls(filenames),'%s '),opt); 
	end
end

end

% load_file loads a single .stat or .data file
function data = load_file(filename,opt)
disp(['loading data from ' filename]);

% distinguishes between stat and data files
if (strcmp(filename(end-4:end),'.data'))
	data = load_data_file(filename);
else
	data = load_stat_file(filename);
end

% \todo{why is this needed?}
%data.label = data.name;
if (~iscell(data))
    data = make_readable(data);
    if ~isfield(opt,'basic')
        data = get_standard_variables(data,opt);
    end
else
    for i=1:length(data)    
        data{i} = make_readable(data{i});
        if ~isfield(opt,'basic')
            data{i} = get_standard_variables(data{i},opt);
        end
    end
end
end

%loads all variables stored in a stat file
function data = load_stat_file(filename)
data.name = filename(1:strfind(filename,'.stat')-1);

% load raw data from file, with the first three lines as header files
% (such that matlab recognizes the array size)
rawdata = importdata(filename,' ',3);

if iscell(rawdata), return; end

if (size(rawdata.data,2)==26)
  % for old .stat files
  disp('WARNING:outdated stat file')
  rawdata.data = [ ...
    rawdata.data(:,1:14) ...
    zeros(size(rawdata.data,1),3) ...
    rawdata.data(:,15:end) ...
    zeros(size(rawdata.data,1),10)];
  data.variable_names = { ...
    'VolumeFraction'; 'Density'; ...
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
  % reads variable names from line 1 and text output from line 2
  variable_names = textscan(rawdata.textdata{1},'%s ');
  data.variable_names = variable_names{1};
end

%also allows input from the restart and ene files if they exist
data = read_restart(filename,data);
data = read_ene(filename,data);

%reads header file of stat file
text = textscan(rawdata.textdata{2},'%s ');
data.text = text{1};

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,3))); size(rawdata.data,1)+1];
data.time = cell2mat(textscan(rawdata.textdata{3},'%f '));
data.w = cell2mat(textscan(rawdata.textdata{2}(3:end),'%f '));

data.coordinates = rawdata.data(1:index_time(1)-1,1:3);
doGradient=any(cellfun(@(d)strcmp(d,'doGradient'),data.text));
doVariance=any(cellfun(@(d)strcmp(d,'doVariance'),data.text));

if (length(index_time)==1)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
elseif (doVariance&&length(index_time)==2)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
    data.variance = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
elseif (doGradient&&length(index_time)==4)
    data.variables = rawdata.data(1:index_time(1)-1,4:end);
    data.gradx = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
    data.grady = rawdata.data(index_time(2)+1:index_time(3)-1,4:end);
    data.gradz = rawdata.data(index_time(3)+1:index_time(4)-1,4:end);
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
end

% \todo{why is this needed?}
%data.description = ['Goldhirsch w/d=' num2str(str2double(data.text{2})/1e-3)];

end

%loads from data file the same variables that are present in the
%stat file (for comparison to Stefan's code)
function data = load_data_file(filename)
data.name = filename(1:end-5);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);

% reads variable names from line 1 and text output
% from line 2
data.variable_names = { ...
	'VolumeFraction'; 'Density'; ...
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
VolumeFraction = max(rawdata.data(:,10),1e-200*ones(size(rawdata.data,1),1));
data.variables = [VolumeFraction ...
	VolumeFraction*rho ...
	rawdata.data(:,51:53).*VolumeFraction(:,[1 1 1])*rho ...
	rawdata.data(:,[41:43 45 46 49]) ...
	nan(size(rawdata.data,1),3) ...
	rawdata.data(:,[21:23 25 26 29]) ...
	rawdata.data(:,[31:33 35 36 39]) ...
	rawdata.data(:,[71:73 75 76 79]) ...
	nan(size(rawdata.data,1),4) ...
	];

%data.description = ['Luding ' data.text];

end

% rewrites coordinates and variables into n dimensional shape where n is
% the VolumeFractionmber of dimensions in stattype (easy for plotting) 
function [data] = make_readable(data)

data.nx=length(data.coordinates(:,1))/sum(data.coordinates(:,1)==data.coordinates(1,1));
data.ny=length(data.coordinates(:,2))/sum(data.coordinates(:,2)==data.coordinates(1,2));
data.nz=length(data.coordinates(:,3))/sum(data.coordinates(:,3)==data.coordinates(1,3));
shape=size(squeeze(zeros(data.nz,data.ny,data.nx)));
if (prod(shape(:))~=size(data.coordinates,1))
    disp('Warning: cannot put xyz values on mesh')
    shape=[size(data.coordinates,1) 1];
end
data.x = reshape(data.coordinates(:,1),shape);
data.y = reshape(data.coordinates(:,2),shape);
data.z = reshape(data.coordinates(:,3),shape);
data = rmfield(data,{'coordinates'});

%the array 'variable' is expanded into smaller arrays, each containing one
%variable only (f.e. Density), with the field names defined in
%variable_names. 
%If available, variance, gradx, grady, gradz are also expanded the same way

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

if isfield(data,'Nu')
   data.VolumeFraction = data.Nu;
   data = rmfield(data,'Nu');
end
end

% extracts a load of extra variables (such as stress) from the basic
% microscopic fields 
function data = get_standard_variables(data,opt)

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

data.CoordinationVolumeFractionmber=(data.FabricXX+data.FabricYY+data.FabricZZ)./data.VolumeFraction;
data.CoordinationVolumeFractionmber(isnan(data.CoordinationVolumeFractionmber)) = 0;

if isfield(opt,'basic'); return; end

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

if (isfield(data,'VolumeFraction_dz'))
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

if (data.nz>1&&data.nx==1&&data.ny==1), 
    data = get_depth_variables(data); 
    if isfield(data,'ParticleDensity')
        [data.MomentumEquationsRemainder, data.MomentumEquationsMaximum]=get_momentum_equation(data);
    end
end

end

% extracts extra variables that only exist for stattype Z
function data = get_depth_variables(data)

dz=diff(data.z(1:2,1));
nz=size(data.z,1);
data.ShearRate = diff(data.VelocityX)./dz;
data.ShearRate(end+1) = data.ShearRate(end);
if isfield(data,'d')
    data.InertialParameter=data.ShearRate*data.d./sqrt(data.StressZZ/data.ParticleDensity);
end

data.FlowVolumeFraction = mean(data.VolumeFraction);
data.FlowDensity = mean(data.Density,1);
data.FlowMomentumX = mean(data.MomentumX,1);
data.FlowMomentumY = mean(data.MomentumY,1);
data.FlowMomentumZ = mean(data.MomentumZ,1);
data.FlowVelocityX = data.FlowMomentumX./data.FlowDensity;
data.FlowVelocityY = data.FlowMomentumY./data.FlowDensity;
data.FlowVelocityZ = data.FlowMomentumZ./data.FlowDensity;

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
  if ~isempty(FlowHeightIndex)&&isfield(data,'Gravity')
    data.Froude=norm([data.FlowVelocityX data.FlowVelocityY data.FlowVelocityZ])/sqrt(data.FlowHeight*(-data.Gravity(3)));
  end
end

end

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

end

function [remainder,maximum] = get_momentum_equation(data)
if isfield(data,'VolumeFraction_dx'),
  NablaMomentumX = data.ParticleDensity*[data.MomentumX_dx data.MomentumX_dy data.MomentumX_dz];
  NablaMomentumY = data.ParticleDensity*[data.MomentumY_dx data.MomentumY_dy data.MomentumY_dz];
  NablaMomentumZ = data.ParticleDensity*[data.MomentumZ_dx data.MomentumZ_dy data.MomentumZ_dz];
  NablaDotStress = [...
      data.StressXX_dx+data.StressXY_dy+data.StressXZ_dz ...
      data.StressYX_dx+data.StressYY_dy+data.StressYZ_dz ...
      data.StressZX_dx+data.StressZY_dy+data.StressZZ_dz ];
else
  disp('estimating gradient!')
  NablaMomentumX = nabla(data.ParticleDensity*data.VolumeFraction.*data.VelocityX,data.x,data.y,data.z);
  NablaMomentumY = nabla(data.ParticleDensity*data.VolumeFraction.*data.VelocityY,data.x,data.y,data.z);
  NablaMomentumZ = nabla(data.ParticleDensity*data.VolumeFraction.*data.VelocityZ,data.x,data.y,data.z);

  NablaDotStress = nabla([data.StressXX data.StressXY data.StressXZ data.StressYY data.StressYZ data.StressZZ],data.x,data.y,data.z);
end
    
VelocityDotNablaMomentum = [...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumX,2) ...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumY,2) ...
  sum([data.VelocityX data.VelocityX data.VelocityZ].*NablaMomentumZ,2) ];

DensityGravity = data.ParticleDensity * data.VolumeFraction * data.Gravity;

Traction = [data.TractionX data.TractionY data.TractionZ ];

remainder = VelocityDotNablaMomentum -DensityGravity +NablaDotStress +Traction;

maximum = max([max(VelocityDotNablaMomentum) max(DensityGravity) max(NablaDotStress)]);
end

function data = read_ene(statname,data)

%load ene data
filename = [statname(1:end-5) '.ene'];
if ~exist(filename,'file'), 
  %if unsuccessful, load ene data with last part of filename removed
  dots=strfind(statname,'.');
  if (length(dots)>1), filename = [statname(1:min(dots(end-1:end))) 'ene']; end
  if ~exist(filename,'file')
    disp([filename ' not found']); 
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
end

end

function data = read_restart(statname,data)

%load restart data
filename = strrep(statname,'.stat','.restart');
fid=fopen(filename);
%if name.restart could not be opened, it maybe because name has an appendix 
%(e.g. problem.X.stat tries to load problem.X.restart). 
%Thus, try loading restart data with last part of name removed 
%(e.g. load problem.restart)
if fid==-1 
  oldfilename=filename;
  dots=strfind(statname,'.');
  if ~isscalar(dots)
    filename = [statname(1:min(dots(end-1:end))) 'restart'];
  end
  fid=fopen(filename);
  if fid~=-1
    disp([oldfilename ' not found, using ' filename ' instead']); 
  end
end
if fid==-1
    disp([filename ' not found']);
    tdata = { ...
    'restart_version','1' ...
    };
else
    tdata = textscan(fid,'%s');
    tdata = tdata{1};
    fclose(fid);
end

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
    if length(i)>1, disp('multiple species detected'); end
	if ~isempty(i), data.ParticleDensity = str2double(tdata(i+1)); end
    
	i = find(strcmp(tdata,'slidingFrictionCoefficient'));
    if isempty(i), i = find(strcmp(tdata,'mu')); end
	if ~isempty(i), data.Mu = str2double(tdata(i+1)); end
	
    i = find(strcmp(tdata,'gravity'));
	if ~isempty(i), data.Gravity = str2double(tdata(i+(1:3)))'; end
	
    i = find(strcmp(tdata,'Particles'));
	if ~isempty(i), data.N = str2double(tdata{i+1}); end
	
    i = find(strcmp(tdata,'xMin'));
	if isempty(i), i = find(strcmp(tdata,'xmin')); end
	if ~isempty(i), data.Domain = str2double(tdata(i+(1:2:11)))'; end
	
    i = find(strcmp(tdata,'FixedParticleRadius'));
	if ~isempty(i), data.FixedParticleRadius = str2double(tdata{i+1}); end
	
    i = find(strcmp(tdata,'MinInflowParticleRadius'));
	if ~isempty(i), 
		data.InflowParticleRadius = str2double(tdata(i+[1 3]))'; 
		data.d = sum(data.InflowParticleRadius);
	else
		i = find(strcmp(tdata,'Particles'));
        if ~isempty(i), 
            if (str2double(tdata{i+1})==0)
                data.d=nan;
            else
                j = min(find(strcmp(tdata(i+2:end),'radius')));
                data.d = 2*str2double(tdata{i+j+2});
            end
        end
	end
end
if isfield(data,'d'), data.ParticleVolume = pi/6*data.d^3; end
if isfield(data,'Domain'), data.DomainVolume = prod(data.Domain([2 4 6])-data.Domain([1 3 5])); end
if isfield(data,'Gravity'), data.ChuteAngle = round(atand(-data.Gravity(1)/data.Gravity(3))*400)/400; end
end

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

% stores differentials (small VolumeFractionmber if X/Y/Z is constant)
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
end

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

% stores differentials (small VolumeFractionmber if X/Y/Z is constant)
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
end
